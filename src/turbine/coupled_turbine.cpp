// coupled_turbine.cpp
// Coupled turbine-in-grid driver (Wave 4): actuator-disk-style turbine <-> fdm3.
//
// Loop structure (staggered, classic body-force coupling):
//
//   for each fluid step:
//     if (step % fluid_steps_per_exchange == 0)  EXCHANGE:
//       1. sample the disk plane (8 azimuths per station) from the fdm3 field
//       2. evaluate local blade-element forces over a grouped SurfaceFlow
//       3. per-station ring torque Q_e / thrust T_e via integrate_moment
//       4. build the Gaussian-smeared body force (fluid receives the NEGATED
//          blade force), under-relaxed and ramped
//       5. advance the rotor by the exchange window against the generator load
//     advance the fdm3 solver one step under the current body force
//
// Sign table (matches the standalone UniformDisk rules):
//   blade thrust magnitude   T_e > 0  acts on the rotor in -e_z (downwind)
//   fluid body force along e_z        = +T_e·η/(ρ·2π·r_e·dr_e)   (decelerates -Z inflow)
//   blade driving torque     Q_e > 0  rotates the rotor to +ω
//   fluid tangential force f_θ        = -Q_e·η/(ρ·2π·r_e²·dr_e)  (counter-swirl)
//   f_x = f_θ·(-sinφ), f_y = f_θ·cosφ   (φ = atan2(y-oy, x-ox))

#include <exd/physics/turbine/coupled_turbine.hpp>

#include <exd/physics/io/series_writer.hpp>

#include <exd/physics/turbine/turbine_simulator.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/physics/fluid/forces/force_evaluator.hpp>
#include <exd/physics/fluid/forces/flow_types.hpp>
#include <exd/physics/mechanics/dynamics.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace exd::physics::turbine
{

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr char kDefaultAirfoil[] = "naca0012";
constexpr int kAzimuthSamples = 8;

using exd::physics::fluid::fdm::FDMBoundaryType;
using exd::physics::fluid::fdm3::BoundaryFace;

/// Find the axial inlet BC, derive V_inf (magnitude of its w_value) and the
/// face it sits on.
bool find_inlet(const fdm3::FDM3Config& grid, const fdm3::FDM3BoundaryCondition*& inlet,
                double& v_inf)
{
    inlet = nullptr;
    for (const auto& bc : grid.boundary_conditions)
    {
        if (bc.type == FDMBoundaryType::Inlet)
        {
            inlet = &bc;
            break;
        }
    }
    if (!inlet) return false;
    if (inlet->face == BoundaryFace::ZMin && inlet->w_value > 0.0)
        v_inf = inlet->w_value;   // upflow inlet
    else
        v_inf = std::fabs(inlet->w_value);
    return true;
}

void merge_warnings(CoupledTurbineResult& result, ModelStatus& status,
                    const std::vector<std::string>& warns)
{
    for (const auto& w : warns)
    {
        if (std::find(result.warnings.begin(), result.warnings.end(), w) ==
            result.warnings.end())
        {
            result.warnings.push_back(w);
        }
        if (std::find(status.warnings.begin(), status.warnings.end(), w) ==
            status.warnings.end())
        {
            status.warnings.push_back(w);
        }
    }
}

bool fail(ModelStatus& status, CoupledTurbineResult& result, const std::string& msg)
{
    status.ok = false;
    status.error = msg;
    result.valid = false;
    result.error = msg;
    return false;
}

} // anonymous namespace

// ────────────────────────────────────────────────────────────────────
// Default grid layout (documented convention in the header)
// ────────────────────────────────────────────────────────────────────

fdm3::FDM3Config default_grid_config(double v_inf, int n_per_axis,
                                     double domain_radius_mult,
                                     double domain_length_mult)
{
    namespace fdm = exd::physics::fluid::fdm;
    fdm3::FDM3Config cfg;
    cfg.nx = n_per_axis;
    cfg.ny = n_per_axis;
    cfg.nz = n_per_axis;
    cfg.lx = domain_radius_mult;   // nominal unit-radius disk × mult
    cfg.ly = domain_radius_mult;
    cfg.lz = domain_length_mult;
    cfg.rho = 1.225;
    cfg.mu = 1.81e-5;
    cfg.dt = 0.25 * std::min({cfg.dx(), cfg.dy(), cfg.dz()}) / std::max(v_inf, 1e-9);
    cfg.time_integration = fdm::TimeIntegration::ForwardEuler;
    cfg.advection_scheme = fdm::AdvectionScheme::Hybrid;
    cfg.pressure_max_iterations = 300;
    cfg.pressure_tolerance = 1e-6;
    cfg.sor_omega = 1.5;
    // Full projection (velocity under-relaxation 1.0): the coupled force field
    // is dynamic, so each step should fully enforce incompressibility.
    cfg.velocity_under_relaxation = 1.0;
    cfg.pressure_under_relaxation = 0.3;
    cfg.initial_u = 0.0;
    cfg.initial_v = 0.0;
    cfg.initial_w = -v_inf;
    cfg.initial_p = 0.0;
    cfg.boundary_conditions = {
        {BoundaryFace::ZMax, FDMBoundaryType::Inlet, 0.0, 0.0, -v_inf},
        {BoundaryFace::ZMin, FDMBoundaryType::Outlet},
        {BoundaryFace::XMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::XMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
    };
    return cfg;
}

// ────────────────────────────────────────────────────────────────────
// Coupled run
// ────────────────────────────────────────────────────────────────────

CoupledTurbineResult run_coupled_turbine(const CoupledTurbineConfig& cfg,
                                         ModelStatus& status)
{
    using exd::physics::fluid::forces::BladeGeometry;
    using exd::physics::fluid::forces::BladeStation;
    using exd::physics::fluid::forces::ForceEvaluatorType;
    using exd::physics::fluid::forces::SurfaceFlow;
    using exd::physics::mechanics::MomentResult;
    using exd::physics::mechanics::RotationalState;

    CoupledTurbineResult result;
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    // ── 1. Grid validation (all fdm3 checks) ──
    {
        std::string err;
        std::vector<std::string> warns;
        if (!cfg.grid.validate(err, warns))
            return fail(status, result, err), result;
        merge_warnings(result, status, warns);
    }

    // ── 2. Turbine geometry (warnings merged) ──
    const double dt = (cfg.dt > 0.0) ? cfg.dt : cfg.grid.dt;
    BladeGeometry blade;
    {
        std::vector<std::string> geo_warnings;
        blade = make_blade_geometry(cfg.turbine, cfg.element_count, kDefaultAirfoil,
                                    geo_warnings);
        merge_warnings(result, status, geo_warnings);
    }
    if (blade.stations.empty() || !(blade.r_tip > blade.r_hub) || blade.r_hub <= 0.0)
    {
        const std::string msg = "invalid blade geometry" +
            (result.warnings.empty() ? std::string() : ": " + result.warnings.back());
        return fail(status, result, msg), result;
    }
    const std::size_t n_stations = blade.stations.size();
    const int blade_count = blade.blade_count;

    // ── 3. Coupling-specific checks ──
    const double window = static_cast<double>(cfg.fluid_steps_per_exchange) * dt;
    if (cfg.element_count < 4)
        return fail(status, result, "coupled turbine: element_count must be >= 4"), result;
    if (cfg.rotor_inertia <= 0.0)
        return fail(status, result, "coupled turbine: rotor_inertia must be positive"), result;
    if (cfg.force_relaxation <= 0.0 || cfg.force_relaxation > 1.0)
        return fail(status, result,
                    "coupled turbine: force_relaxation must be in (0, 1]"), result;
    if (cfg.smear_cells < 1.0)
        return fail(status, result, "coupled turbine: smear_cells must be >= 1"), result;
    if (cfg.fluid_steps_per_exchange < 1)
        return fail(status, result,
                    "coupled turbine: fluid_steps_per_exchange must be >= 1"), result;
    if (cfg.max_steps < 1)
        return fail(status, result, "coupled turbine: max_steps must be >= 1"), result;
    if (cfg.history_interval < 1)
        return fail(status, result,
                    "coupled turbine: history_interval must be >= 1"), result;
    if (cfg.dt < 0.0)
        return fail(status, result,
                    "coupled turbine: dt must be >= 0 (0 = use grid.dt)"), result;
    if (dt <= 0.0)
        return fail(status, result,
                    "coupled turbine: time step must be positive"), result;
    if (cfg.ramp_time_s < 10.0 * window)
    {
        return fail(status, result,
            "coupled turbine: ramp_time_s must be >= 10 * fluid_steps_per_exchange * dt "
            "(forcing-discontinuity guard)"), result;
    }

    // ── 4. Inlet / V_inf ──
    const fdm3::FDM3BoundaryCondition* inlet = nullptr;
    double v_inf = 0.0;
    if (!find_inlet(cfg.grid, inlet, v_inf))
        return fail(status, result,
                    "coupled turbine needs an Inlet BC (w_value = -V_inf)"), result;
    if (v_inf <= 0.0)
        return fail(status, result,
                    "coupled turbine: inlet axial velocity must be non-zero"), result;
    {
        const bool ok_face = inlet->face == BoundaryFace::ZMax ||
                             (inlet->face == BoundaryFace::ZMin && inlet->w_value > 0.0);
        if (!ok_face)
        {
            const std::string w = "coupled turbine: inlet is not on the expected "
                "axial face (ZMax with w_value = -V_inf, or ZMin with w_value = +V_inf); "
                "results follow the documented -Z inflow convention";
            merge_warnings(result, status, {w});
        }
    }

    // ── 5. Stability guard: the smeared force layer must be resolved ──
    {
        const double hmin = std::min({cfg.grid.dx(), cfg.grid.dy(), cfg.grid.dz()});
        if (dt >= cfg.smear_cells * hmin / (3.0 * v_inf))
        {
            return fail(status, result,
                "coupled turbine: dt too large; must be < "
                "smear_cells*min(dx,dy,dz)/(3*V_inf) so the smeared force "
                "layer is resolved"), result;
        }
    }

    // ── 6. Rotor-disc containment (warnings) ──
    {
        const double ox = cfg.rotor_origin[0];
        const double oy = cfg.rotor_origin[1];
        const double oz = cfg.rotor_origin[2];
        const double r_tip = blade.r_tip;
        if (ox - r_tip < 0.0 || ox + r_tip > cfg.grid.lx ||
            oy - r_tip < 0.0 || oy + r_tip > cfg.grid.ly)
        {
            merge_warnings(result, status, {
                "coupled turbine: rotor disc not fully inside the lateral grid "
                "bounds; disk samples near the wall read zero velocity (stall "
                "guard acts)"});
        }
        if (oz < 0.1 * cfg.grid.lz || oz > 0.9 * cfg.grid.lz)
        {
            merge_warnings(result, status, {
                "coupled turbine: rotor plane is poorly centered along lz; "
                "keep it ~0.1..0.9 of the axial domain"});
        }
    }

    // ── Setup: fluid solver ──
    fdm3::FDM3Solver solver;
    {
        ModelStatus fstatus;
        if (!solver.initialize(cfg.grid, fstatus))
        {
            result.valid = false;
            result.error = fstatus.error;
            status.ok = false;
            status.error = fstatus.error;
            merge_warnings(result, status, fstatus.warnings);
            return result;
        }
        merge_warnings(result, status, fstatus.warnings);
    }

    // ── Rotor machine-state CSV (optional; one flushed row per step) ──
    std::unique_ptr<io::CsvSeriesWriter> csv;
    if (!cfg.csv_path.empty())
    {
        csv = std::make_unique<io::CsvSeriesWriter>(
            cfg.csv_path,
            std::vector<std::string>{"omega_rad_s", "angle_rad", "torque_Nm",
                                     "axial_force_N", "power_W", "exchange"},
            true /* flush per row: real-time visible */, &status);
        if (!status.ok)
        {
            result.valid = false;
            result.error = status.error;
            return result;
        }
    }

    // Axis: rotor plane maps to grid z = rotor_origin[2]; the blade geometry's
    // own z_rotor axial offset is folded into the axis origin so the per-blade
    // reference points land on the plane.
    const double ox = cfg.rotor_origin[0];
    const double oy = cfg.rotor_origin[1];
    const double oz = cfg.rotor_origin[2];
    const mechanics::RotationAxis axis{{ox, oy, oz - blade.z_rotor}, {0.0, 0.0, 1.0}};

    // Force evaluator: local blade element, created ONCE and reused across
    // exchanges (no per-exchange evaluator allocation).
    fluid::forces::PolarDatabase polars;
    polars.add_builtin_polars();
    fluid::forces::ForceEvaluatorParams fparams;
    fparams.type = ForceEvaluatorType::BladeElement;
    fparams.polars = &polars;
    std::unique_ptr<fluid::forces::IForceEvaluator> evaluator =
        fluid::forces::make_force_evaluator(fparams, status);

    // Generator load + rotational dynamics.
    std::unique_ptr<mechanics::IMomentModel> load;
    if (!cfg.generator.omega_pts.empty())
        load = mechanics::make_curve_moment(cfg.generator);
    std::unique_ptr<mechanics::IRotationalDynamics> dynamics =
        mechanics::make_rigid_rotor_dynamics(
            {cfg.rotor_inertia, mechanics::RotationalIntegration::Heun});

    RotationalState state{0.0, 0.0};

    // Initial fluid kinetic energy (analytic, uniform initial field).
    const double dV = cfg.grid.dx() * cfg.grid.dy() * cfg.grid.dz();
    const double vol_box = static_cast<double>(cfg.grid.nx) * cfg.grid.ny *
                           cfg.grid.nz * dV;
    const double ke_start_volume =
        0.5 * cfg.grid.rho * (cfg.grid.initial_u * cfg.grid.initial_u +
                              cfg.grid.initial_v * cfg.grid.initial_v +
                              cfg.grid.initial_w * cfg.grid.initial_w) * vol_box;

    // Body-force buffers (interior cells, i + nx*(j + ny*k)); reused.
    const std::size_t cells = static_cast<std::size_t>(cfg.grid.nx) *
                              cfg.grid.ny * cfg.grid.nz;
    std::vector<double> fx(cells, 0.0), fy(cells, 0.0), fz(cells, 0.0);
    std::vector<double> fx_old(cells, 0.0), fy_old(cells, 0.0), fz_old(cells, 0.0);

    // Sampling geometry + SurfaceFlow buffers (reused across exchanges; the
    // vectors are sized once and only the VALUES are refreshed).
    const std::size_t n_samp = n_stations * static_cast<std::size_t>(kAzimuthSamples);
    SurfaceFlow flow;
    flow.points.resize(n_samp, {0.0, 0.0, 0.0});
    flow.normals.resize(n_samp, {0.0, 0.0, 0.0});
    flow.velocity.resize(n_samp, {0.0, 0.0, 0.0});
    flow.shear_traction.resize(n_samp, {0.0, 0.0, 0.0});
    flow.pressure.resize(n_samp, 0.0);
    flow.area.resize(n_samp, 0.0);
    flow.element_index.resize(n_samp, 0);
    flow.density = cfg.grid.rho;
    flow.viscosity = cfg.grid.mu;
    flow.p_ref = 0.0;
    for (std::size_t e = 0; e < n_stations; ++e)
        for (int k = 0; k < kAzimuthSamples; ++k)
        {
            const std::size_t sidx = e * static_cast<std::size_t>(kAzimuthSamples) +
                                     static_cast<std::size_t>(k);
            flow.element_index[sidx] = static_cast<int32_t>(e);
        }

    // Work / exchange accounting.
    double aero_work = 0.0;
    double load_work = 0.0;
    std::vector<double> power_samples;
    power_samples.reserve(cfg.max_steps / cfg.fluid_steps_per_exchange + 1);
    double total_torque = 0.0;
    double total_axial = 0.0;
    uint64_t exchanges = 0;

    const double hmin = std::min({cfg.grid.dx(), cfg.grid.dy(), cfg.grid.dz()});
    const double eps = cfg.smear_cells * hmin;
    const double eta_norm = 1.0 / (eps * std::sqrt(2.0 * kPi));
    const double rho = cfg.grid.rho;
    const double beta = cfg.force_relaxation;
    const double tau_ramp = cfg.ramp_time_s;

    std::vector<mechanics::ElementForce3D> per_element;
    per_element.reserve(n_stations * static_cast<std::size_t>(blade_count));

    // ── Main loop ─────────────────────────────────────────────────────
    for (uint64_t step = 0; step < cfg.max_steps; ++step)
    {
        const double step_t = static_cast<double>(step) * dt;

        if (step % cfg.fluid_steps_per_exchange == 0)
        {
            // 1. Sample disk plane (8 azimuths per station) + refresh the
            //    SurfaceFlow velocity/pressure arrays through the adapter.
            for (std::size_t e = 0; e < n_stations; ++e)
            {
                const double r_e = std::max(blade.stations[e].r, 1e-9);
                for (int k = 0; k < kAzimuthSamples; ++k)
                {
                    const double phi = state.angle_rad +
                        2.0 * kPi * static_cast<double>(k) / kAzimuthSamples;
                    const std::size_t sidx = e * static_cast<std::size_t>(kAzimuthSamples) +
                                             static_cast<std::size_t>(k);
                    flow.points[sidx] = {
                        ox + r_e * std::cos(phi),
                        oy + r_e * std::sin(phi),
                        oz,
                    };
                }
            }
            // Fresh adapter snapshot of the current solver field.
            std::unique_ptr<coupling::IFlowField3D> adapter =
                coupling::make_fdm3_field_adapter(solver);
            for (std::size_t s = 0; s < n_samp; ++s)
            {
                std::array<double, 3> vel{0.0, 0.0, 0.0};
                double p = 0.0;
                if (adapter->sample(flow.points[s], vel, p))
                {
                    flow.velocity[s] = vel;
                    flow.pressure[s] = p;
                }
                else
                {
                    // Out-of-bounds sample (disk near the lateral wall): zero
                    // velocity so the station feeds W=0 and the stall guard.
                    flow.velocity[s] = {0.0, 0.0, 0.0};
                    flow.pressure[s] = 0.0;
                }
            }

            // 2. Evaluate blade-element forces (grouped by element_index).
            per_element.clear();
            evaluator->compute(blade, flow, state.omega, axis, per_element, status);
            if (!status.ok)
            {
                result.valid = false;
                result.error = status.error;
                return result;
            }
            merge_warnings(result, status, status.warnings);

            // 3. Per-station ring torque / thrust.
            total_torque = 0.0;
            total_axial = 0.0;
            std::vector<double> station_torque(n_stations, 0.0);
            std::vector<double> station_axial(n_stations, 0.0);
            for (std::size_t e = 0; e < n_stations; ++e)
            {
                const std::size_t b0 = e * static_cast<std::size_t>(blade_count);
                const MomentResult m = mechanics::integrate_moment(
                    std::span<const mechanics::ElementForce3D>(
                        &per_element[b0], static_cast<std::size_t>(blade_count)),
                    axis);
                if (!m.valid)
                    return fail(status, result,
                                "axis degenerate in moment integration"), result;
                station_torque[e] = m.torque;
                station_axial[e] = m.axial_force;
                total_torque += m.torque;
                total_axial += m.axial_force;
            }

            // 4. Build the smeared body force (fluid receives -blade force).
            std::fill(fx.begin(), fx.end(), 0.0);
            std::fill(fy.begin(), fy.end(), 0.0);
            std::fill(fz.begin(), fz.end(), 0.0);
            for (std::size_t e = 0; e < n_stations; ++e)
            {
                const BladeStation& st = blade.stations[e];
                const double r_e = std::max(st.r, 1e-9);
                const double dr_e = std::max(st.dr, 1e-9);
                // Ring totals: T_e = thrust magnitude (blade thrust is -e_z);
                // Q_e = ring torque (driving +ω).
                const double T_e = -station_axial[e];
                const double Q_e = station_torque[e];
                const double az_denom = rho * 2.0 * kPi * r_e * dr_e;
                const double th_denom = rho * 2.0 * kPi * r_e * r_e * dr_e;
                if (az_denom <= 0.0 || th_denom <= 0.0) continue;
                const double band_half = 0.5 * dr_e;
                const double band2 = band_half * band_half;

                for (int k = 0; k < cfg.grid.nz; ++k)
                {
                    const double z = (static_cast<double>(k) + 0.5) * cfg.grid.dz();
                    const double zt = z - oz;
                    if (zt > 4.0 * eps || zt < -4.0 * eps) continue;
                    const double eta = std::exp(-zt * zt / (2.0 * eps * eps)) * eta_norm;
                    for (int j = 0; j < cfg.grid.ny; ++j)
                    {
                        const double y = (static_cast<double>(j) + 0.5) * cfg.grid.dy();
                        for (int i = 0; i < cfg.grid.nx; ++i)
                        {
                            const double x = (static_cast<double>(i) + 0.5) * cfg.grid.dx();
                            const double rho_cell = std::sqrt((x - ox) * (x - ox) +
                                                               (y - oy) * (y - oy));
                            if (rho_cell < 1e-12) continue;
                            const double drho = rho_cell - r_e;
                            if (drho * drho > band2) continue;
                            const std::size_t id = static_cast<std::size_t>(i) +
                                static_cast<std::size_t>(cfg.grid.nx) *
                                (static_cast<std::size_t>(j) +
                                 static_cast<std::size_t>(cfg.grid.ny) * k);
                            fz[id] += T_e * eta / az_denom;
                            const double f_theta = -Q_e * eta / th_denom;
                            const double phi_cell = std::atan2(y - oy, x - ox);
                            fx[id] += f_theta * (-std::sin(phi_cell));
                            fy[id] += f_theta * std::cos(phi_cell);
                        }
                    }
                }
            }

            // Under-relax + ramp (same ramp factor applied to the relaxed
            // combination each exchange).
            const double ramp = std::min(1.0, step_t / tau_ramp);
            for (std::size_t c = 0; c < cells; ++c)
            {
                fx[c] = ramp * (beta * fx[c] + (1.0 - beta) * fx_old[c]);
                fy[c] = ramp * (beta * fy[c] + (1.0 - beta) * fy_old[c]);
                fz[c] = ramp * (beta * fz[c] + (1.0 - beta) * fz_old[c]);
            }
            fx_old = fx;
            fy_old = fy;
            fz_old = fz;
            {
                ModelStatus fstatus;
                if (!solver.set_body_force(fx, fy, fz, fstatus))
                {
                    result.valid = false;
                    result.error = fstatus.error;
                    status.ok = false;
                    status.error = fstatus.error;
                    return result;
                }
            }

            // 5. Rotor advance over the frozen-force window.
            const double omega_before = state.omega;
            double external = 0.0;
            if (load)
            {
                external = load->moment(state, status);
                if (!status.ok)
                {
                    result.valid = false;
                    result.error = status.error;
                    return result;
                }
            }
            const double net = total_torque - external;
            state = dynamics->advance(window, net, state, status);
            if (!status.ok)
            {
                result.valid = false;
                result.error = status.error;
                return result;
            }

            // Work over the frozen-force window uses the trapezoidal mean of
            // the window's end omegas (a more accurate estimate than the
            // pre-advance value; keeps the energy bookkeeping consistent with
            // the rotor KE change reported by the dynamics).
            const double omega_avg = 0.5 * (omega_before + state.omega);
            aero_work += total_torque * omega_avg * window;
            load_work += external * omega_avg * window;
            power_samples.push_back(total_torque * omega_before);
            ++exchanges;
        }

        // Fluid step under the current body force.
        {
            ModelStatus fstatus;
            if (!solver.step(dt, fstatus))
            {
                result.valid = false;
                result.error = fstatus.error;
                status.ok = false;
                status.error = fstatus.error;
                merge_warnings(result, status, fstatus.warnings);
                return result;
            }
        }

        const double t_now = solver.time();

        // History.
        if (cfg.record_history && (step % cfg.history_interval == 0))
        {
            CoupledTurbineStep hs;
            hs.t = t_now;
            hs.omega = state.omega;
            hs.angle_rad = state.angle_rad;
            hs.torque = total_torque;
            hs.axial_force = total_axial;
            hs.power = total_torque * state.omega;
            hs.exchange = exchanges;
            result.history.push_back(hs);
        }

        // Rotor machine-state CSV (real-time: one flushed row per step).
        if (csv) csv->write_row(t_now, std::vector<double>{state.omega, state.angle_rad,
                                                           total_torque, total_axial,
                                                           total_torque * state.omega,
                                                           static_cast<double>(exchanges)});

        // Field stamping (drivers own writer stamps).
        if (cfg.field_writer)
        {
            bool emit = false;
            if (cfg.output_scheduler)
            {
                cfg.output_scheduler->set_now(t_now);
                emit = cfg.output_scheduler->should_emit(step);
            }
            else
            {
                emit = (step % cfg.field_stamp_interval == 0);
            }
            if (emit)
            {
                const auto& field = solver.field();
                io::FieldGeometry geo;
                geo.origin = {0.5 * cfg.grid.dx(), 0.5 * cfg.grid.dy(), 0.5 * cfg.grid.dz()};
                geo.spacing = {cfg.grid.dx(), cfg.grid.dy(), cfg.grid.dz()};
                geo.dims = {static_cast<uint32_t>(cfg.grid.nx),
                            static_cast<uint32_t>(cfg.grid.ny),
                            static_cast<uint32_t>(cfg.grid.nz)};

                std::vector<float> vel(3 * cells, 0.0f);
                std::vector<float> pres(cells, 0.0f);
                for (int k = 0; k < cfg.grid.nz; ++k)
                    for (int j = 0; j < cfg.grid.ny; ++j)
                        for (int i = 0; i < cfg.grid.nx; ++i)
                        {
                            const std::size_t id = field.index(i, j, k);
                            vel[3 * id + 0] = static_cast<float>(field.u[id]);
                            vel[3 * id + 1] = static_cast<float>(field.v[id]);
                            vel[3 * id + 2] = static_cast<float>(field.w[id]);
                            pres[id] = static_cast<float>(field.p[id]);
                        }
                if (!cfg.field_writer->begin_stamp(t_now, static_cast<uint64_t>(step + 1)))
                {
                    merge_warnings(result, status,
                                   {"coupled turbine: field stamp begin failed"});
                }
                else
                {
                    cfg.field_writer->write_vector_field("velocity", geo, vel);
                    cfg.field_writer->write_scalar_field("pressure", geo, pres);
                    cfg.field_writer->end_stamp();
                }
            }
        }
    }

    if (csv) csv->close();

    // ── Final summary ────────────────────────────────────────────────
    result.fluid.valid = true;
    result.fluid.field = solver.field();
    result.fluid.steps_taken = solver.step_count();
    result.fluid.final_time = solver.time();
    result.fluid.history = {solver.last_step()};
    result.fluid.converged = false;

    result.valid = true;
    result.final_omega = state.omega;
    result.final_tsr = (v_inf > 0.0) ? state.omega * blade.r_tip / v_inf : 0.0;
    result.exchanges = exchanges;
    result.aero_work = aero_work;
    result.load_work = load_work;

    // cp from the mean aero power over the LAST 20% of exchanges.
    {
        const std::size_t n = power_samples.size();
        const std::size_t n_take = std::max<std::size_t>(1, (n * 20 + 99) / 100);
        const std::size_t start = (n > n_take) ? n - n_take : 0;
        double sum = 0.0;
        for (std::size_t i = start; i < n; ++i) sum += power_samples[i];
        const double p_mean = (n > 0) ? sum / static_cast<double>(n - start) : 0.0;
        const double area = kPi * blade.r_tip * blade.r_tip;
        const double power_denom = 0.5 * rho * area * v_inf * v_inf * v_inf;
        result.final_cp = (power_denom > 0.0) ? p_mean / power_denom : 0.0;
    }

    result.rotor_ke_change = 0.5 * cfg.rotor_inertia * (state.omega * state.omega);

    // Final fluid kinetic energy (interior cells) -> fluid_ke_change.
    {
        const auto& field = result.fluid.field;
        double ke = 0.0;
        for (int k = 0; k < cfg.grid.nz; ++k)
            for (int j = 0; j < cfg.grid.ny; ++j)
                for (int i = 0; i < cfg.grid.nx; ++i)
                {
                    const std::size_t id = field.index(i, j, k);
                    ke += field.u[id] * field.u[id] + field.v[id] * field.v[id] +
                          field.w[id] * field.w[id];
                }
        const double ke_end_volume = 0.5 * rho * ke * dV;
        result.fluid_ke_change = ke_end_volume - ke_start_volume;
    }

    return result;
}

} // namespace exd::physics::turbine