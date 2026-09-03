// gas_force_test.cpp
// Cycle pressure/temperature model: polytropic compression,
// Wiebe heat release, steam admission/expansion, continuity.

#include <doctest/doctest.h>

#include <cmath>

#include "../../../../src/presets/engine/engine_internal.hpp"

using namespace exd::engine::presets::engine;

namespace
{
constexpr double PI = 3.14159265358979323846;
constexpr double DEG = PI / 180.0;

EngineGeometryConfig geo_fixture()
{
    EngineGeometryConfig g;
    g.crank_radius = 0.05;
    g.rod_length = 0.20;
    g.bore = 0.086;
    g.clearance_volume = 1.0e-4;
    g.piston_mass = 0.5;
    g.flywheel_inertia = 0.02;
    return g;
}

EngineThermoConfig otto_fixture()
{
    EngineThermoConfig t;
    t.cycle = EngineCycleType::Otto;
    t.q_in_cycle = 0.0; // motored unless a test raises it
    return t;
}
} // anonymous namespace

TEST_CASE("otto: compression phase follows the polytrope pV^gc exactly")
{
    auto g = geo_fixture();
    auto t = otto_fixture();
    const double A = piston_area(g);
    const double v_bdc = g.clearance_volume + A * 2.0 * g.crank_radius;

    for (double deg = 540.0; deg < 720.0; deg += 30.0)
    {
        const double theta = deg * DEG;
        const double v = cylinder_volume(theta, g);
        double p = 0.0, T = 0.0;
        cylinder_state(theta, g, t, p, T);
        const double expected = t.p_intake * std::pow(v_bdc / v, t.gamma_compression);
        CHECK(p == doctest::Approx(expected).epsilon(1e-9));
    }
}

TEST_CASE("otto: TDC pressure is p_intake·r_c^gc (Wiebe heat zero at ignition start)")
{
    auto g = geo_fixture();
    auto t = otto_fixture();
    const double A = piston_area(g);
    const double v_tdc = g.clearance_volume;
    const double v_bdc = g.clearance_volume + A * 2.0 * g.crank_radius;
    const double r_c = v_bdc / v_tdc;
    const double p_expected = t.p_intake * std::pow(r_c, t.gamma_compression);

    double p = 0.0, T = 0.0;
    cylinder_state(0.0, g, t, p, T);
    CHECK(p == doctest::Approx(p_expected).epsilon(1e-9));
    // Heat release raises pressure mid-power-stroke and completes by 60°.
    t.q_in_cycle = 1500.0;
    cylinder_state(0.0, g, t, p, T);
    CHECK(p == doctest::Approx(p_expected).epsilon(1e-9)); // xb(0°) = 0
    double p_hot = 0.0, T_hot = 0.0;
    cylinder_state(30.0 * DEG, g, t, p_hot, T_hot);
    CHECK(p_hot > p);
    double p_end = 0.0, T_end = 0.0;
    cylinder_state(120.0 * DEG, g, t, p_end, T_end); // after burn: expansion only
    double p_poly = 0.0, T_poly = 0.0;
    t.q_in_cycle = 0.0;
    cylinder_state(120.0 * DEG, g, t, p_poly, T_poly);
    CHECK(p_end > p_poly); // heat stays in the expansion
}

TEST_CASE("otto: pressure is continuous across phase boundaries")
{
    auto g = geo_fixture();
    auto t = otto_fixture();
    t.q_in_cycle = 1500.0;

    auto p_at = [&](double deg)
    {
        double p = 0.0, T = 0.0;
        cylinder_state(deg * DEG, g, t, p, T);
        return p;
    };
    // Exhaust ramp joins expansion end to p_exhaust.
    const double p_180 = p_at(180.0);
    CHECK(std::fabs(p_180 - p_at(180.0 - 1.0)) < 0.02 * p_180);
    // Intake ramp joins p_exhaust to p_intake.
    CHECK(p_at(360.0) == doctest::Approx(t.p_exhaust).epsilon(1e-6));
    CHECK(p_at(540.0) == doctest::Approx(t.p_intake).epsilon(1e-6));
    // Compression meets the power stroke at TDC.
    CHECK(p_at(720.0 - 1.0) == doctest::Approx(p_at(0.0)).epsilon(1e-3));
}

TEST_CASE("otto: temperature via trapped mass pV = mRT")
{
    auto g = geo_fixture();
    auto t = otto_fixture();
    t.q_in_cycle = 1500.0;
    double p = 0.0, T = 0.0;
    cylinder_state(45.0 * DEG, g, t, p, T);
    const double A = piston_area(g);
    const double v_bdc = g.clearance_volume + A * 2.0 * g.crank_radius;
    const double m = t.p_intake * v_bdc / (t.r_gas * t.T_intake);
    const double expected_T = p * cylinder_volume(45.0 * DEG, g) / (m * t.r_gas);
    CHECK(T == doctest::Approx(expected_T).epsilon(1e-9));
    CHECK(T > t.T_intake); // compressed + heated
}

TEST_CASE("steam: admission, cutoff expansion, exhaust")
{
    auto g = geo_fixture();
    EngineThermoConfig t;
    t.cycle = EngineCycleType::Steam;
    t.p_boiler = 800000.0;
    t.p_condenser = 15000.0;
    t.steam_cutoff_deg = 40.0;
    t.steam_gamma = 1.13;

    double p = 0.0, T = 0.0;
    // Admission: full boiler pressure well inside the window.
    cylinder_state(20.0 * DEG, g, t, p, T);
    CHECK(p == doctest::Approx(t.p_boiler).epsilon(1e-6));
    // Expansion follows (V_cut/V)^n from the cutoff volume.
    const double v_cut = cylinder_volume(40.0 * DEG, g);
    double p_exp = 0.0, T_exp = 0.0;
    cylinder_state(90.0 * DEG, g, t, p_exp, T_exp);
    const double v90 = cylinder_volume(90.0 * DEG, g);
    CHECK(p_exp == doctest::Approx(t.p_boiler * std::pow(v_cut / v90, t.steam_gamma)).epsilon(1e-9));
    // Exhaust: condenser pressure late in the stroke.
    double p_ex = 0.0, T_ex = 0.0;
    cylinder_state(300.0 * DEG, g, t, p_ex, T_ex);
    CHECK(p_ex == doctest::Approx(t.p_condenser).epsilon(1e-6));
    // Cycle repeat: theta + 2π is identical.
    double p2 = 0.0, T2 = 0.0;
    cylinder_state(300.0 * DEG + 2.0 * PI, g, t, p2, T2);
    CHECK(p2 == doctest::Approx(p_ex).epsilon(1e-12));
}

TEST_CASE("load moment: constant + viscous + generator curve")
{
    EngineLoadConfig load;
    load.friction_constant = 0.5;
    load.friction_viscous = 0.1;
    CHECK(load_moment(load, 10.0) == doctest::Approx(1.5).epsilon(1e-12));

    load.generator_enabled = true;
    load.generator_omega_pts = {0.0, 100.0, 200.0};
    load.generator_torque_pts = {0.0, 20.0, 30.0};
    CHECK(load_moment(load, 50.0) == doctest::Approx(10.0 + 0.5 + 0.1 * 50.0).epsilon(1e-12));
    CHECK(load_moment(load, 250.0) == doctest::Approx(30.0 + 0.5 + 0.1 * 250.0).epsilon(1e-12));
    CHECK(load_moment(load, 0.0) == doctest::Approx(0.5).epsilon(1e-12));
}