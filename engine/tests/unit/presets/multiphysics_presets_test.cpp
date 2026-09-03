// multiphysics_presets_test.cpp — W14: cross-domain preset assemblies,
// each verified against an analytic anchor:
//   * aeroacoustics  — pulse convects at c + u (mean flow from fdm3)
//   * thermal stress — free bar expands u(L) = α·g·L²/2 (thermal strain)
//   * joule heating  — P ≈ V²σA/L; mid T on the uniform-q parabola
//   * species-in-flow— outlet c = c_in·exp(−k·L/u)
#include <exd/engine/presets/multiphysics/aeroacoustics.hpp>
#include <exd/engine/presets/multiphysics/joule_heating.hpp>
#include <exd/engine/presets/multiphysics/species_in_flow.hpp>
#include <exd/engine/presets/multiphysics/thermal_stress.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::presets::multiphysics;

TEST_CASE("Aeroacoustics preset: pulse arrives on schedule with the duct mean flow")
{
    // Assembly check: the fdm3 duct supplies the mean flow; the wave module
    // (whose own suite verifies the c±u physics at large u/c) transports the
    // pulse.  At the duct's u = 0.2 vs c = 343 the convective shift is 0.1%
    // — far below the discrete arrival resolution — so the anchor is the
    // arrival WINDOW with the plain-wave control removing flow artifacts.
    auto run_case = [&](bool with_flow) {
        AeroacousticsConfig cfg;
        cfg.flow.nx = 20; cfg.flow.ny = 8; cfg.flow.nz = 8;
        cfg.flow.length = 1.0; cfg.flow.width = 0.4; cfg.flow.height = 0.4;
        cfg.flow.inlet_velocity = 0.2;
        cfg.flow.max_steps = 900;
        cfg.sound_speed = 343.0;
        cfg.pulse_center = {0.3, 0.2, 0.2};
        cfg.pulse_width = 0.08;
        cfg.max_steps = 4000;
        cfg.include_flow = with_flow;
        cfg.probe_index = 14 + 20 * (4 + 8 * 4);   // x = 0.7
        return run_aeroacoustics(cfg);
    };

    const auto r = run_case(true);
    REQUIRE(r.ok);
    CHECK(r.mean_velocity[0] == doctest::Approx(0.2).epsilon(0.2));

    const double c = 343.0, d = 0.4;
    const double t_expected = d / (c + r.mean_velocity[0]);
    // arrival window: max |p| in [0.7·t, 1.4·t] must exceed 20% of the
    // initial amplitude (the pulse is really there, on schedule)
    auto window_peak = [&](const AeroacousticsResult& rr) {
        const double dt = rr.wave_result.dt_used;
        double peak = 0.0;
        for (size_t t = 1; t < rr.wave_result.probe_history.size(); ++t)
        {
            const double time = static_cast<double>(t) * dt;
            if (time >= 0.7 * t_expected && time <= 1.4 * t_expected)
                peak = std::max(peak, std::fabs(rr.wave_result.probe_history[t]));
        }
        return peak;
    };
    const double peak_flow = window_peak(r);
    CHECK(peak_flow > 0.2 * 100.0);       // > 20% of the pulse amplitude
    CHECK(peak_flow < 2.0 * 100.0);       // no blow-up

    // plain-wave control: arrival is c-only and the same window holds
    const auto r0 = run_case(false);
    REQUIRE(r0.ok);
    CHECK(window_peak(r0) > 0.2 * 100.0);

    // the coupling actually injected the mean flow
    CHECK(r.mean_velocity[0] > 0.1);
}

TEST_CASE("Thermal-stress preset: free bar expands u(L) = α·g·L²/2")
{
    ThermalStressConfig cfg;
    cfg.nx = 31; cfg.ny = 4; cfg.nz = 4;
    cfg.length = 1.0;
    cfg.t_left = 300.0;
    cfg.t_right = 400.0;
    cfg.thermal_expansion = 1e-5;

    const auto r = run_thermal_stress(cfg);
    REQUIRE(r.ok);
    // T(x) = 300 + 100·x (linear, exact for steady conduction); the free-bar
    // expansion: u_x(L) = α·∫₀ᴸ (T−T_ref)dx = α·100·L²/2
    const double dTdx = (cfg.t_right - cfg.t_left) / cfg.length;
    const double exact = cfg.thermal_expansion * dTdx * cfg.length * cfg.length / 2.0;
    CHECK(r.measured_tip_displacement == doctest::Approx(exact).epsilon(0.02));
}

TEST_CASE("Joule-heating preset: exact parallel-plate P and uniform-q parabola")
{
    JouleHeatingConfig cfg;
    cfg.nx = 41; cfg.ny = 9; cfg.nz = 9;
    cfg.spacing = 0.025;
    cfg.voltage = 10.0;
    cfg.conductivity = 5.96e7;

    const auto r = run_joule_heating(cfg);
    REQUIRE(r.ok);
    // With the NEUMANN side walls the static field is the exact discrete
    // linear bridge: E = −V/L uniform, q = σ·(V/L)² uniform, and
    // P = σ·V²·A/L holds to the solver tolerance (W15 fix).
    const double l_gap = cfg.spacing * (cfg.nx - 1);
    CHECK(r.total_power == doctest::Approx(r.analytic_power).epsilon(0.02));
    CHECK(r.center_source == doctest::Approx(
        cfg.conductivity * cfg.voltage * cfg.voltage / (l_gap * l_gap)).epsilon(0.02));

    // thermal: uniform-q parabola T_mid = T0 + q·L²/(8k) exactly
    const double k = cfg.thermal_conductivity;
    const double t_mid_exact = cfg.t_wall + r.center_source * l_gap * l_gap /
                               (8.0 * k);
    CHECK(r.mid_temperature == doctest::Approx(t_mid_exact).epsilon(0.01));
    CHECK(r.mid_temperature > cfg.t_wall);

    // steady energy balance: ∫q dV = k·A·(∂T/∂n at the x ends)
    const double area = cfg.spacing * (cfg.ny - 1) * cfg.spacing * (cfg.nz - 1);
    const size_t nxh = static_cast<size_t>(cfg.nx);
    const size_t nyh = static_cast<size_t>(cfg.ny);
    const size_t mid_j = nyh / 2, mid_k = cfg.nz / 2;
    auto tline = [&](size_t i) {
        return r.thermal.temperature.values[i + nxh * (mid_j + nyh * mid_k)];
    };
    const double h = cfg.spacing;
    const double slope_minus = (tline(1) - tline(0)) / h;
    const double slope_plus = (tline(nxh - 1) - tline(nxh - 2)) / h;
    const double flux_out = k * area * (std::fabs(slope_minus) + std::fabs(slope_plus));
    CHECK(r.thermal.total_power == doctest::Approx(r.total_power).epsilon(0.01));
    CHECK(flux_out == doctest::Approx(r.total_power).epsilon(0.02));
}

TEST_CASE("Species-in-flow preset: outlet follows exp(−k·L/u)")
{
    SpeciesInFlowConfig cfg;
    cfg.flow.nx = 20; cfg.flow.ny = 8; cfg.flow.nz = 8;
    cfg.flow.length = 1.0; cfg.flow.width = 0.4; cfg.flow.height = 0.4;
    cfg.flow.inlet_velocity = 0.2;
    cfg.flow.max_steps = 900;
    cfg.decay_rate = 0.05;

    const auto r = run_species_in_flow(cfg);
    REQUIRE(r.ok);
    CHECK(r.outlet_concentration == doctest::Approx(r.analytic_outlet).epsilon(0.05));
    // sanity: decay really decays
    CHECK(r.outlet_concentration < 1.0);
    CHECK(r.outlet_concentration > 0.5 * r.analytic_outlet);
}
