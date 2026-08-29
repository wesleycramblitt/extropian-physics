// turbine_simulator_test.cpp
// Unit tests for the turbine application assembly: blade-geometry mapping,
// generator load curve, the single coupled step, and the uniform-freestream
// time-marching simulation.

#include <exd/physics/turbine/turbine_simulator.hpp>

#include <exd/geometry/turbine.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/fluid/forces/flow_types.hpp>
#include <exd/physics/mechanics/moment_model.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>
#include <exd/physics/mechanics/status.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace exd::geometry;
using namespace exd::physics::turbine;

namespace
{

using exd::physics::fluid::forces::ForceEvaluatorType;
using exd::physics::fluid::forces::Freestream;
using exd::physics::mechanics::ModelStatus;
using exd::physics::mechanics::RotationalState;

// ── Fixtures ────────────────────────────────────────────────────────

TurbineDefinition make_turbine()
{
    TurbineDefinition t;
    t.flow_path.hub_points = {{0.0f, 0.5f}, {1.0f, 0.5f}};
    t.flow_path.shroud_points = {{0.0f, 5.0f}, {1.0f, 5.0f}};

    BladeRow row;
    row.type = BladeRowType::Rotor;
    row.blade_count.value = 3.0f;
    row.rotational_speed.value = 30.0f; // rpm; informational (simulator starts from config)
    row.leading_edge_hub = {0.3f, 0.5f};
    row.leading_edge_shroud = {0.3f, 5.0f};
    row.trailing_edge_hub = {0.8f, 0.5f};
    row.trailing_edge_shroud = {0.8f, 5.0f};

    BladeSection s0;
    s0.span = 0.0f;
    s0.chord.value = 0.4f; // real chord comes from the LE/TE distances below
    s0.stagger.value = -2.0f;
    BladeSection s1 = s0;
    s1.span = 1.0f;
    row.sections = {s0, s1};

    t.blade_rows = {row};
    return t;
}

Freestream inflow()
{
    Freestream fs;
    fs.velocity = {0.0, 0.0, -10.0}; // downwind axial inflow; machine axis = +Z
    fs.rho = 1.225;
    fs.mu = 1.81e-5;
    fs.p_ref = 101325.0;
    return fs;
}

TurbineConfig step_config()
{
    TurbineConfig cfg;
    cfg.element_count = 16;
    cfg.inertia = 500.0;
    cfg.dt = 0.05;
    cfg.force.type = ForceEvaluatorType::MomentumBalance;
    cfg.generator = {}; // empty curve = no load
    return cfg;
}

std::unique_ptr<exd::physics::coupling::IFlowField3D> uniform_field(const Freestream& fs)
{
    return exd::physics::coupling::make_uniform_field(
        exd::physics::coupling::UniformFieldConfig{fs.velocity, fs.rho, fs.mu, fs.p_ref});
}

} // anonymous namespace

// ── make_blade_geometry ─────────────────────────────────────────────

TEST_CASE("make_blade_geometry: maps the single rotor row onto blade stations")
{
    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(make_turbine(), 16, "naca0012", warnings);

    REQUIRE_FALSE(blade.stations.empty());
    CHECK(blade.stations.size() == 16U);

    // r_hub = 0.5*(0.5 + 0.5), r_tip = 0.5*(5 + 5), z_rotor = 0.5*(0.3 + 0.8).
    CHECK(blade.r_hub == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(blade.r_tip == doctest::Approx(5.0).epsilon(1e-6));
    CHECK(blade.z_rotor == doctest::Approx(0.55).epsilon(1e-6));
    CHECK(blade.blade_count == 3);

    const auto& first = blade.stations.front();
    CHECK(first.dr == doctest::Approx(4.5 / 16.0).epsilon(1e-6));
    CHECK(first.r == doctest::Approx(0.5 + 0.5 * 4.5 / 16.0).epsilon(1e-6));
    // |te - le| = sqrt(0.5^2 + 0^2) = 0.5 for both hub and shroud.
    CHECK(first.chord == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(first.twist_deg == doctest::Approx(-2.0).epsilon(1e-6));
    CHECK(first.thickness_ratio == doctest::Approx(0.1).epsilon(1e-6)); // fixture default
    CHECK(first.airfoil == "naca0012");

    double prev_r = -1.0;
    for (const auto& s : blade.stations)
    {
        CHECK(s.r > prev_r);
        CHECK(s.twist_deg == doctest::Approx(-2.0).epsilon(1e-6));
        prev_r = s.r;
    }
}

TEST_CASE("make_blade_geometry: no rotor row yields empty default and warning")
{
    auto t = make_turbine();
    t.blade_rows[0].type = BladeRowType::Stator;

    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);

    CHECK(blade.stations.empty());
    CHECK(blade.r_hub == 0.0);
    CHECK(blade.r_tip == 0.0);
    bool found = false;
    for (const auto& w : warnings)
        if (w.find("no Rotor") != std::string::npos) found = true;
    CHECK(found);
}

TEST_CASE("make_blade_geometry: multiple rotor rows yield empty default and warning")
{
    auto t = make_turbine();
    const BladeRow rotor_copy = t.blade_rows.front();
    t.blade_rows.push_back(rotor_copy);

    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);

    CHECK(blade.stations.empty());
    bool found = false;
    for (const auto& w : warnings)
        if (w.find("Rotor rows") != std::string::npos) found = true;
    CHECK(found);
}

TEST_CASE("make_blade_geometry: non-rotor rows are ignored with a warning")
{
    auto t = make_turbine();
    BladeRow stator;
    stator.type = BladeRowType::Stator;
    t.blade_rows = {stator, t.blade_rows.front()};

    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);

    REQUIRE_FALSE(blade.stations.empty());
    bool found = false;
    for (const auto& w : warnings)
        if (w.find("non-rotor") != std::string::npos) found = true;
    CHECK(found);
}

TEST_CASE("make_blade_geometry: inverted or non-positive radii are fatal")
{
    auto t = make_turbine();
    t.blade_rows[0].trailing_edge_shroud = {0.8f, -5.0f}; // r_tip = 0
    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);
    CHECK(blade.stations.empty());
    CHECK(blade.r_tip == 0.0);
}

TEST_CASE("make_blade_geometry: empty sections give zero stagger with warning")
{
    auto t = make_turbine();
    t.blade_rows[0].sections.clear();

    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);

    REQUIRE(blade.stations.size() == 8U);
    CHECK(blade.stations.front().twist_deg == 0.0);
    CHECK(blade.stations.front().thickness_ratio == doctest::Approx(0.12));
    bool found = false;
    for (const auto& w : warnings)
        if (w.find("sections empty") != std::string::npos) found = true;
    CHECK(found);
}

TEST_CASE("make_blade_geometry: stagger clamps flat with extrapolation warnings")
{
    auto t = make_turbine();
    t.blade_rows[0].sections.resize(1);
    t.blade_rows[0].sections[0].span = 0.5f;
    t.blade_rows[0].sections[0].stagger.value = -4.0f;

    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(t, 8, "naca0012", warnings);

    REQUIRE(blade.stations.size() == 8U);
    CHECK(blade.stations.front().twist_deg == doctest::Approx(-4.0).epsilon(1e-6));
    int count = 0;
    for (const auto& w : warnings)
        if (w.find("extrapolated") != std::string::npos) ++count;
    CHECK(count >= 1);
}

TEST_CASE("make_blade_geometry: element_count < 4 is rejected")
{
    std::vector<std::string> warnings;
    auto blade = make_blade_geometry(make_turbine(), 2, "naca0012", warnings);
    CHECK(blade.stations.empty());
}

// ── make_generator_curve ────────────────────────────────────────────

TEST_CASE("make_generator_curve: samples the P/(eta*omega) hyperbola")
{
    auto curve = make_generator_curve(1000.0, 0.9, 1.0, 64);
    REQUIRE(curve.omega_pts.size() == 64U);

    CHECK(curve.omega_pts.front() == doctest::Approx(1.0).epsilon(1e-12));
    // The first sample sits exactly on the flat-floor torque value.
    CHECK(curve.torque_pts.front() == doctest::Approx(1000.0 / 0.9).epsilon(1e-9));

    // Every sample reproduces the analytic hyperbola at its own omega.
    for (std::size_t i = 0; i < curve.omega_pts.size(); ++i)
        CHECK(curve.torque_pts[i] ==
              doctest::Approx(1000.0 / (0.9 * curve.omega_pts[i])).epsilon(1e-6));

    // Strict monotonicity of both columns.
    for (std::size_t i = 1; i < curve.omega_pts.size(); ++i)
    {
        CHECK(curve.omega_pts[i] > curve.omega_pts[i - 1]);
        CHECK(curve.torque_pts[i] < curve.torque_pts[i - 1]);
    }

    // High end spans 4x the reference omega.
    const double w_ref = std::max(1.0, 1000.0 / (0.9 * 1.0));
    CHECK(curve.omega_pts.back() == doctest::Approx(4.0 * w_ref).epsilon(1e-6));
}

TEST_CASE("make_generator_curve: moment model clamps flat below min_omega")
{
    auto curve = make_generator_curve(1000.0, 0.9, 1.0, 32);
    auto model = exd::physics::mechanics::make_curve_moment(curve);

    RotationalState below{0.5, 0.0};
    ModelStatus status;
    const double t_low = model->moment(below, status);
    CHECK(status.ok);
    CHECK(t_low == doctest::Approx(1000.0 / 0.9).epsilon(1e-9));
}

TEST_CASE("make_generator_curve: invalid parameters return the empty no-load curve")
{
    CHECK(make_generator_curve(0.0, 0.9, 1.0).omega_pts.empty());
    CHECK(make_generator_curve(1000.0, 0.0, 1.0).omega_pts.empty());
    CHECK(make_generator_curve(1000.0, 0.9, 0.0).omega_pts.empty());
}

// ── step_turbine ────────────────────────────────────────────────────

TEST_CASE("step_turbine: accelerates a rotor from rest in a uniform inflow")
{
    auto turbine = make_turbine();
    auto cfg = step_config();
    auto field = uniform_field(inflow());

    RotationalState state{0.0, 0.0};
    auto r1 = step_turbine(state, turbine, *field, cfg);
    REQUIRE(r1.ok);
    CHECK(state.omega > 0.0);
    CHECK(r1.state.omega > 0.0);
    CHECK(r1.aero.valid);
    CHECK(r1.aero.torque > 0.0);
    CHECK(r1.aero.axial_force < 0.0); // thrust acts downwind (-e_z)
    CHECK(r1.per_element.size() ==
          static_cast<std::size_t>(cfg.element_count) * 3U);
    for (const auto& f : r1.per_element)
        for (double c : f.force)
            CHECK(std::isfinite(c));
    // Power from rest is zero (omega_before == 0).
    CHECK(r1.power == doctest::Approx(0.0).epsilon(1e-9));

    // A second step keeps accelerating; power is now positive.
    const double w1 = state.omega;
    auto r2 = step_turbine(state, turbine, *field, cfg);
    REQUIRE(r2.ok);
    CHECK(state.omega > w1);
    CHECK(r2.power > 0.0);
}

TEST_CASE("step_turbine: invalid config fails without touching state")
{
    auto cfg = step_config();
    cfg.element_count = 2; // < 4
    auto field = uniform_field(inflow());

    RotationalState state{3.0, 1.0};
    auto r = step_turbine(state, make_turbine(), *field, cfg);
    CHECK_FALSE(r.ok);
    CHECK_FALSE(r.status.error.empty());
    CHECK(state.omega == doctest::Approx(3.0));
    CHECK(state.angle_rad == doctest::Approx(1.0));
}

TEST_CASE("step_turbine: invalid blade geometry fails and leaves state untouched")
{
    auto t = make_turbine();
    t.blade_rows[0].type = BladeRowType::Stator; // no rotor row
    auto cfg = step_config();
    auto field = uniform_field(inflow());

    RotationalState state{2.0, 0.5};
    auto r = step_turbine(state, t, *field, cfg);
    CHECK_FALSE(r.ok);
    CHECK_FALSE(r.status.error.empty());
    CHECK(state.omega == doctest::Approx(2.0));
    CHECK(state.angle_rad == doctest::Approx(0.5));
}

// ── simulate_turbine ────────────────────────────────────────────────

TEST_CASE("simulate_turbine: spins up and settles against a matched generator")
{
    auto turbine = make_turbine();

    TurbineConfig cfg;
    cfg.element_count = 16;
    cfg.dt = 0.02;
    cfg.max_steps = 1500;
    cfg.inertia = 100.0; // low inertia so the rotor reaches equilibrium inside the window
    cfg.initial_omega = 0.0;
    cfg.record_history = true;
    cfg.history_interval = 1;
    cfg.generator = make_generator_curve(5000.0, 0.9, 20.0, 64);

    auto result = simulate_turbine(turbine, inflow(), cfg);

    REQUIRE(result.valid);
    CHECK(result.final_step.ok);
    CHECK(result.final_step.state.omega > 0.0);
    CHECK(result.tsr > 0.0);
    CHECK(result.cp > 0.0);
    CHECK(result.cp <= 16.0 / 27.0 + 1e-9); // Betz limit in settled operation
    CHECK(result.total_energy > 0.0);
    CHECK(result.history.size() == static_cast<std::size_t>(cfg.max_steps));

    CHECK(std::isfinite(result.final_step.power));
    CHECK(std::isfinite(result.final_step.state.angle_rad));

    // Settled: the last recorded omega is within 1% of the previous one.
    const auto& h = result.history;
    REQUIRE(h.size() >= 2U);
    const double w_last = h.back().state.omega;
    const double w_prev = h[h.size() - 2].state.omega;
    REQUIRE(w_prev > 0.0);
    CHECK(std::fabs(w_last - w_prev) / w_prev < 0.01);
}

TEST_CASE("simulate_turbine: zero freestream velocity is invalid")
{
    auto fs = inflow();
    fs.velocity = {0.0, 0.0, 0.0};

    auto result = simulate_turbine(make_turbine(), fs, TurbineConfig{});
    CHECK_FALSE(result.valid);
    CHECK_FALSE(result.error.empty());
}

TEST_CASE("simulate_turbine: a single step records exactly one history entry")
{
    TurbineConfig cfg;
    cfg.element_count = 8;
    cfg.dt = 0.02;
    cfg.max_steps = 1;
    cfg.inertia = 100.0;
    cfg.generator = make_generator_curve(5000.0, 0.9, 20.0, 16);

    auto result = simulate_turbine(make_turbine(), inflow(), cfg);
    REQUIRE(result.valid);
    CHECK(result.history.size() == 1U);
}