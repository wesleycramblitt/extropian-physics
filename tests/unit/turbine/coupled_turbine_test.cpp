// coupled_turbine_test.cpp
// Unit tests for the actuator-disk-style coupled turbine-in-grid driver:
// spin-up soak + energy accounting, agreement with the reduced-order
// MomentumBalance reference, determinism, validation guards, and field
// writer stamping.

#include <exd/physics/turbine/coupled_turbine.hpp>

#include <exd/physics/turbine/turbine_simulator.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/physics/fluid/forces/force_evaluator.hpp>
#include <exd/physics/io/field_writer.hpp>
#include <exd/physics/io/output_policy.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

using namespace exd::physics;
using namespace exd::physics::turbine;

namespace
{

using exd::physics::fluid::fdm::FDMBoundaryType;
using exd::physics::fluid::fdm3::BoundaryFace;
using exd::physics::fluid::fdm3::FDM3Config;

constexpr double kPi = 3.14159265358979323846;

// ── Fixtures ────────────────────────────────────────────────────────

// Small axial rotor: r_hub = 0.1, r_tip = 0.4, chord 0.2, 3 blades.
exd::geometry::TurbineDefinition make_turbine()
{
    exd::geometry::TurbineDefinition t;
    t.flow_path.hub_points = {{0.0f, 0.1f}, {1.0f, 0.1f}};
    t.flow_path.shroud_points = {{0.0f, 0.4f}, {1.0f, 0.4f}};

    exd::geometry::BladeRow row;
    row.type = exd::geometry::BladeRowType::Rotor;
    row.blade_count.value = 3.0f;
    row.leading_edge_hub = {0.25f, 0.1f};
    row.leading_edge_shroud = {0.25f, 0.4f};
    row.trailing_edge_hub = {0.45f, 0.1f};
    row.trailing_edge_shroud = {0.45f, 0.4f};

    exd::geometry::BladeSection s0;
    s0.span = 0.0f;
    s0.chord.value = 0.2f;
    s0.stagger.value = -2.0f;
    exd::geometry::BladeSection s1 = s0;
    s1.span = 1.0f;
    row.sections = {s0, s1};

    t.blade_rows = {row};
    return t;
}

// 20³ box sized 6×R_tip lateral / 8×R_tip axial for R_tip = 0.4 m.
FDM3Config make_grid(double v_inf)
{
    auto g = default_grid_config(v_inf, 20);
    g.lx = 2.4;
    g.ly = 2.4;
    g.lz = 3.2;
    g.rho = 1.225;
    g.mu = 1.81e-5;
    g.dt = 0.02;
    g.time_integration = exd::physics::fluid::fdm::TimeIntegration::RK4;
    g.advection_scheme = exd::physics::fluid::fdm::AdvectionScheme::Hybrid;
    g.pressure_max_iterations = 300;
    g.pressure_tolerance = 1e-6;
    g.sor_omega = 1.5;
    return g;
}

CoupledTurbineConfig base_config()
{
    CoupledTurbineConfig cfg;
    cfg.grid = make_grid(1.0);
    cfg.turbine = make_turbine();
    cfg.element_count = 16;
    cfg.rotor_inertia = 0.02;
    cfg.rotor_origin = {1.2, 1.2, 0.45 * cfg.grid.lz};
    cfg.fluid_steps_per_exchange = 10;
    cfg.force_relaxation = 0.4;
    cfg.ramp_time_s = 2.5;      // ≥ 10·window = 10·(10·0.02) = 2.0
    cfg.smear_cells = 2.5;
    cfg.dt = 0.02;
    cfg.max_steps = 1500;
    cfg.history_interval = 1;
    // Opposing torque curve sized so Q(ω) ≈ T_load around ω ≈ 8-12 rad/s.
    cfg.generator.omega_pts = {0.0, 5.0, 10.0, 15.0, 20.0};
    cfg.generator.torque_pts = {0.0, 0.004, 0.009, 0.012, 0.012};
    return cfg;
}

// One-shot reduced-order (MomentumBalance) cp at a prescribed ω against a
// uniform -Z freestream — the standalone reference for the coupled disk.
double reduced_order_cp(const exd::geometry::TurbineDefinition& turbine,
                        double omega, double v_inf, double rho, double mu,
                        int element_count)
{
    TurbineConfig tc;
    tc.element_count = element_count;
    tc.default_airfoil = "naca0012";
    tc.force.type = exd::physics::fluid::forces::ForceEvaluatorType::MomentumBalance;

    coupling::UniformFieldConfig uf;
    uf.velocity = {0.0, 0.0, -v_inf};
    uf.rho = rho;
    uf.mu = mu;
    uf.p_ref = 0.0;
    auto field = coupling::make_uniform_field(uf);

    std::vector<std::string> gwarnings;
    const auto blade = make_blade_geometry(turbine, element_count, "naca0012", gwarnings);
    const double R = blade.r_tip;
    const double area = kPi * R * R;
    const double pden = 0.5 * rho * area * v_inf * v_inf * v_inf;

    mechanics::RotationalState state{omega, 0.0};
    auto step = step_turbine(state, turbine, *field, tc);
    if (!step.ok || pden <= 0.0) return 0.0;
    return step.power / pden;
}

std::string temp_dir()
{
    auto base = std::filesystem::temp_directory_path();
    auto d = base / ("exd_coupled_turbine_test_" + std::to_string(::getpid()));
    std::error_code ec;
    std::filesystem::remove_all(d, ec);
    std::filesystem::create_directories(d, ec);
    return d.string();
}

} // anonymous namespace

// ── Local blade-element evaluator ───────────────────────────────────

TEST_CASE("blade_element: aggregates station velocities, drives rotor")
{
    auto turbine = make_turbine();
    std::vector<std::string> gwarnings;
    const auto blade = make_blade_geometry(turbine, 16, "naca0012", gwarnings);
    REQUIRE(blade.stations.size() == 16U);

    // 8 sampled velocities per station, all uniform -Z inflow.
    using exd::physics::fluid::forces::SurfaceFlow;
    SurfaceFlow flow;
    constexpr int n_azi = 8;
    const std::size_t n = blade.stations.size() * static_cast<std::size_t>(n_azi);
    flow.points.resize(n);
    flow.normals.resize(n);
    flow.velocity.resize(n);
    flow.shear_traction.resize(n);
    flow.pressure.resize(n);
    flow.area.resize(n);
    flow.element_index.resize(n);
    flow.density = 1.225;
    flow.viscosity = 1.81e-5;
    flow.p_ref = 0.0;
    for (std::size_t e = 0; e < blade.stations.size(); ++e)
        for (int k = 0; k < n_azi; ++k)
        {
            const std::size_t s = e * static_cast<std::size_t>(n_azi) + static_cast<std::size_t>(k);
            flow.element_index[s] = static_cast<int32_t>(e);
            flow.velocity[s] = {0.0, 0.0, -1.0};
        }

    exd::physics::fluid::forces::ForceEvaluatorParams params;
    params.type = exd::physics::fluid::forces::ForceEvaluatorType::BladeElement;
    ModelStatus status;
    auto evaluator = exd::physics::fluid::forces::make_force_evaluator(params, status);
    REQUIRE(evaluator);

    mechanics::RotationAxis axis{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    std::vector<mechanics::ElementForce3D> per_element;
    evaluator->compute(blade, flow, 10.0, axis, per_element, status);
    REQUIRE(status.ok);

    // Same output shape as the other coefficient models: stations · blades.
    CHECK(per_element.size() == blade.stations.size() * 3U);
    for (const auto& f : per_element)
        for (double c : f.force)
            CHECK(std::isfinite(c));

    // Thrust into the inflow (-e_z), torque in the direction of omega (+e_z).
    const mechanics::MomentResult m = mechanics::integrate_moment(per_element, axis);
    REQUIRE(m.valid);
    CHECK(m.axial_force < 0.0);
    CHECK(m.torque > 0.0);

    // The mean is taken over ALL points of a station: zero out one of the 8
    // azimuthal samples and the result changes smoothly (no velocity[0] use).
    flow.velocity[0] = {0.0, 0.0, 0.0};
    std::vector<mechanics::ElementForce3D> per2;
    evaluator->compute(blade, flow, 10.0, axis, per2, status);
    REQUIRE(status.ok);
    CHECK(std::abs(per2[0].force[2]) < std::abs(per_element[0].force[2]));
}

// ── Spin-up soak ────────────────────────────────────────────────────

TEST_CASE("coupled turbine: spin-up soak settles with consistent energy balance")
{
    auto cfg = base_config();
    ModelStatus status;
    auto result = run_coupled_turbine(cfg, status);

    REQUIRE(result.valid);
    REQUIRE(status.ok);
    REQUIRE(result.exchanges > 10);
    REQUIRE(result.history.size() >= 2U);

    // The rotor spins up and settles: last-quarter omega spread is small.
    CHECK(result.final_omega > 0.0);
    const auto& h = result.history;
    const std::size_t q0 = (3 * h.size()) / 4;
    double lo = h[q0].omega, hi = h[q0].omega, sum = 0.0;
    for (std::size_t i = q0; i < h.size(); ++i)
    {
        lo = std::min(lo, h[i].omega);
        hi = std::max(hi, h[i].omega);
        sum += h[i].omega;
    }
    const double mean_q = sum / static_cast<double>(h.size() - q0);
    REQUIRE(mean_q > 0.0);
    CHECK((hi - lo) / mean_q < 0.25);

    // Dimensionless performance in the physical range.
    CHECK(result.final_cp > 0.0);
    CHECK(result.final_cp <= 16.0 / 27.0 + 1e-9); // Betz limit
    CHECK(result.final_tsr > 0.0);

    // Sign table check: the disk decelerates the -Z inflow (inflow enters at
    // the ZMax inlet and travels toward smaller z), so the mean axial speed
    // DOWNSTREAM of the rotor plane is smaller than upstream (wake deficit),
    // measured over the rotor disc area.
    {
        const auto& field = result.fluid.field;
        const double ox = cfg.rotor_origin[0];
        const double oy = cfg.rotor_origin[1];
        const double oz = cfg.rotor_origin[2];
        const double r_tip = 0.4;
        const double r_tip2 = r_tip * r_tip;
        double w_up = 0.0, w_down = 0.0;
        int n_up = 0, n_down = 0;
        for (int k = 0; k < cfg.grid.nz; ++k)
            for (int j = 0; j < cfg.grid.ny; ++j)
                for (int i = 0; i < cfg.grid.nx; ++i)
                {
                    const double x = (static_cast<double>(i) + 0.5) * cfg.grid.dx();
                    const double y = (static_cast<double>(j) + 0.5) * cfg.grid.dy();
                    const double z = (static_cast<double>(k) + 0.5) * cfg.grid.dz();
                    const double dxp = x - ox;
                    const double dyp = y - oy;
                    if (dxp * dxp + dyp * dyp > r_tip2) continue;
                    const std::size_t id = field.index(i, j, k);
                    if (z > oz)          // between the inlet (ZMax) and the disk
                    {
                        w_up += std::fabs(field.w[id]);
                        ++n_up;
                    }
                    else if (z < oz)     // between the disk and the outlet (ZMin)
                    {
                        w_down += std::fabs(field.w[id]);
                        ++n_down;
                    }
                }
        REQUIRE(n_up > 0);
        REQUIRE(n_down > 0);
        INFO("upstream |w|=", w_up / n_up, " downstream |w|=", w_down / n_down);
        CHECK(w_down / n_down < w_up / n_up);
    }

    for (const auto& s : h)
    {
        CHECK(std::isfinite(s.omega));
        CHECK(std::isfinite(s.torque));
        CHECK(std::isfinite(s.power));
    }

    // Energy balance: aero_work covers rotor KE + load work + fluid KE change
    // (with a fudge for the ramp transient).
    const double sink = result.rotor_ke_change + result.load_work +
                        result.fluid_ke_change -
                        std::max(0.05 * result.aero_work, 1e-6);
    INFO("aero_work=", result.aero_work, " rotor_ke=", result.rotor_ke_change,
         " load_work=", result.load_work, " fluid_ke=", result.fluid_ke_change,
         " sink=", sink);
    CHECK(result.aero_work >= sink);

    CHECK(std::isfinite(result.aero_work));
    CHECK(std::isfinite(result.fluid_ke_change));
}

// ── Agreement with the reduced-order reference ──────────────────────

TEST_CASE("coupled turbine: disk cp agrees with MomentumBalance reference")
{
    auto cfg = base_config();
    ModelStatus status;
    auto result = run_coupled_turbine(cfg, status);
    REQUIRE(result.valid);

    const double v_inf = 1.0;
    const double cp_ref = reduced_order_cp(cfg.turbine, result.final_omega,
                                           v_inf, cfg.grid.rho, cfg.grid.mu,
                                           cfg.element_count);
    REQUIRE(cp_ref > 0.0);
    const double ratio = result.final_cp / cp_ref;
    INFO("coupled cp=", result.final_cp, " reference cp=", cp_ref,
         " ratio=", ratio, " omega=", result.final_omega);
    CHECK(ratio >= 0.4);
    CHECK(ratio <= 2.0);
}

// ── Determinism ─────────────────────────────────────────────────────

TEST_CASE("coupled turbine: deterministic given the same config")
{
    auto cfg = base_config();
    ModelStatus s1, s2;
    auto r1 = run_coupled_turbine(cfg, s1);
    auto r2 = run_coupled_turbine(cfg, s2);
    REQUIRE(r1.valid);
    REQUIRE(r2.valid);
    CHECK(std::fabs(r1.final_omega - r2.final_omega) < 1e-9);
    CHECK(r1.exchanges == r2.exchanges);
    CHECK(r1.history.size() == r2.history.size());
}

// ── Validation guards ───────────────────────────────────────────────

TEST_CASE("coupled turbine: grid without an inlet is rejected")
{
    auto cfg = base_config();
    cfg.grid.boundary_conditions.clear();
    cfg.grid.boundary_conditions = {
        {BoundaryFace::ZMax, FDMBoundaryType::Outlet},
        {BoundaryFace::ZMin, FDMBoundaryType::Outlet},
        {BoundaryFace::XMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::XMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
    };
    ModelStatus status;
    auto result = run_coupled_turbine(cfg, status);
    CHECK_FALSE(result.valid);
    CHECK_FALSE(result.error.empty());
}

TEST_CASE("coupled turbine: dt above the smear-layer stability guard is rejected")
{
    auto cfg = base_config();
    cfg.dt = 0.1;             // >= smear_cells*min(dx,dy,dz)/(3*V_inf) = 0.1
    cfg.ramp_time_s = 12.0;   // keep the ramp guard from firing first: window=1, 10*1=10
    ModelStatus status;
    auto result = run_coupled_turbine(cfg, status);
    CHECK_FALSE(result.valid);
    CHECK_FALSE(result.error.empty());
}

TEST_CASE("coupled turbine: short ramp below 10 exchange windows is rejected")
{
    auto cfg = base_config();
    cfg.ramp_time_s = 0.5;    // < 10·window = 10·(10·0.02) = 2.0
    ModelStatus status;
    auto result = run_coupled_turbine(cfg, status);
    CHECK_FALSE(result.valid);
    CHECK_FALSE(result.error.empty());
}

// ── Field writer stamps ─────────────────────────────────────────────

TEST_CASE("coupled turbine: field writer emits stamps through the scheduler")
{
    auto cfg = base_config();
    cfg.max_steps = 220;
    const std::string dir = temp_dir();

    ModelStatus status;
    auto writer = io::make_fld_writer({dir}, status);
    REQUIRE(status.ok);
    REQUIRE(writer);

    io::OutputScheduler sched(io::OutputPolicy{50, 0.0});
    cfg.field_writer = writer.get();
    cfg.output_scheduler = &sched;

    auto result = run_coupled_turbine(cfg, status);
    REQUIRE(result.valid);

    // timeline.txt exists and has at least one stamp line.
    const std::string tl = dir + "/timeline.txt";
    {
        std::ifstream f(tl);
        REQUIRE(f.good());
        int lines = 0;
        std::string line;
        while (std::getline(f, line)) ++lines;
        CHECK(lines >= 1);
    }
    // The first emphasized stamp (step 0 → step_00000001.fld) exists.
    CHECK(std::filesystem::exists(dir + "/step_00000001.fld"));
}
// ── Rotor machine-state CSV (Wave 7) ───────────────────────────────

TEST_CASE("coupled: rotor CSV streams one row per fluid step")
{
    auto cfg = base_config();
    cfg.max_steps = 220;
    cfg.csv_path = std::filesystem::temp_directory_path().string()
                   + "/coupled_rotor_" + std::to_string(::getpid()) + ".csv";
    ModelStatus status;
    auto r = run_coupled_turbine(cfg, status);
    REQUIRE(r.valid);

    std::ifstream f(cfg.csv_path);
    REQUIRE(f.good());
    std::string line;
    int rows = 0;
    bool header_ok = false;
    while (std::getline(f, line))
    {
        if (rows == 0) header_ok = line.rfind("time,omega_rad_s", 0) == 0;
        ++rows;
    }
    CHECK(header_ok);
    CHECK(rows == static_cast<int>(cfg.max_steps) + 1); // header + per-step rows
}

// ── Water turbine capability (Wave 8) ──────────────────────────────

TEST_CASE("coupled: water turbine (hydro) runs with seawater properties")
{
    // Same capability as the wind case — only the working fluid changes:
    // seawater rho/mu, current speed. The full stack (3D FDM field, BCs,
    // blade-element coupling, rotor dynamics) is fluid-agnostic.
    auto cfg = base_config();
    const double v_cur = 2.0;
    cfg.grid = default_grid_config(v_cur, 16, 3.0, 4.0);
    cfg.grid.rho = 1025.0;       // seawater (kg/m³)
    cfg.grid.mu = 1.08e-3;       // seawater dynamic viscosity (Pa·s)
    cfg.rotor_origin = {1.2, 1.2, 1.2};
    cfg.dt = 0.0;                // grid dt (CFL-based)
    cfg.fluid_steps_per_exchange = 10;
    cfg.force_relaxation = 0.4;
    cfg.ramp_time_s = 4.0; // >= 10·window (window = 10 steps × 0.019 s)
    cfg.max_steps = 600;
    // Aero torque scales with ρ: the load curve AND the rotor inertia must
    // scale with the fluid density to land at the same ω transient as air.
    // (Physical: a water turbine's rotor inertia and braking torque are
    //  ~ρ_w/ρ_a ≈ 837× those of an equal-size air rotor.)
    const double rho_ratio = 1025.0 / 1.225;
    cfg.rotor_inertia = 0.02 * rho_ratio;
    cfg.generator.omega_pts = {0.0, 5.0, 10.0, 15.0, 20.0};
    cfg.generator.torque_pts = {0.0, 0.004 * rho_ratio, 0.009 * rho_ratio,
                                0.012 * rho_ratio, 0.012 * rho_ratio};

    ModelStatus status;
    auto r = run_coupled_turbine(cfg, status);
    REQUIRE(r.valid);

    CHECK(r.final_omega > 0.0);
    CHECK(r.final_tsr > 0.0);
    CHECK(r.final_cp > 0.0);
    CHECK(r.final_cp < 16.0 / 27.0 + 1e-6); // Betz bound
    CHECK(r.exchanges > 5U);
    for (const auto& h : r.history)
        CHECK(std::isfinite(h.omega));
}
