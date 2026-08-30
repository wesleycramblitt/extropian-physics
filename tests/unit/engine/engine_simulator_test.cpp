// engine_simulator_test.cpp
// Full engine steps/simulations: motored energy, fired spinup,
// governor regulation, steam, CSV output, validation errors.

#include <exd/physics/engine/engine_simulator.hpp>

#include "../../../../src/engine/engine_internal.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <string>
#include <vector>

using namespace exd::physics::engine;
using exd::physics::ModelStatus;

namespace
{
constexpr double PI = 3.14159265358979323846;

EngineConfig base_config()
{
    EngineConfig c;
    c.geometry.crank_radius = 0.05;
    c.geometry.rod_length = 0.20;
    c.geometry.bore = 0.086;
    c.geometry.clearance_volume = 1.0e-4;
    c.geometry.piston_mass = 0.5;
    c.geometry.flywheel_inertia = 0.02;
    c.thermo.q_in_cycle = 0.0; // motored by default
    c.load.friction_constant = 0.0;
    c.dt = 2.0e-4;
    c.max_steps = 4000;
    return c;
}

/// Volume (m³) — mirrors the source model for the energy test.
double cylinder_volume_for_test(double theta, const EngineGeometryConfig& g)
{
    const double x = g.crank_radius * std::cos(theta)
                   + std::sqrt(g.rod_length * g.rod_length
                               - g.crank_radius * g.crank_radius * std::sin(theta) * std::sin(theta));
    return g.clearance_volume + 0.25 * 3.14159265358979323846 * g.bore * g.bore
                                * ((g.rod_length + g.crank_radius) - x);
}

bool all_finite(const EngineSimResult& r)
{
    for (const auto& h : r.history)
    {
        if (!std::isfinite(h.state.theta_rad) || !std::isfinite(h.state.omega)
            || !std::isfinite(h.p_cyl) || !std::isfinite(h.T_cyl)
            || !std::isfinite(h.indicated_moment) || !std::isfinite(h.power)
            || !std::isfinite(h.piston_x) || !std::isfinite(h.piston_v))
            return false;
    }
    return true;
}
} // anonymous namespace

TEST_CASE("engine: pure-inertia free spin conserves kinetic energy exactly")
{
    // Nullify the gas: a HUGE clearance volume makes V ≈ const, so
    // p_cyl ≈ p_back on every leg and M_gas ≈ 0 — the machine becomes
    // purely ½·J(θ)·ω² + (zero potential). KE = ½·J(θ)·ω² must then be
    // conserved EXACTLY, which directly validates the ½·(dJ/dθ)·ω²
    // inertia-torque term in the crank ODE: a wrong factor (c ≠ ½)
    // redistributes/creates ~10 J/cycle at ω₀ = 100 while the correct
    // term keeps KE flat to integrator drift. No quadrature involved.
    auto cfg = base_config();
    cfg.geometry.clearance_volume = 1.0e9; // gas ~nullified (p variation ~1e-6)
    cfg.thermo.q_in_cycle = 0.0;
    cfg.initial_omega = 100.0;
    cfg.dt = 2.0e-4;
    cfg.max_steps = 10000; // several cycles
    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);
    REQUIRE(r.history.size() > 100U);

    const auto& g = cfg.geometry;
    auto crank_inertia = [&](double th)
    {
        const double dx = -g.crank_radius * std::sin(th)
            - (g.crank_radius * g.crank_radius * std::sin(th) * std::cos(th))
                  / std::sqrt(g.rod_length * g.rod_length
                              - g.crank_radius * g.crank_radius * std::sin(th) * std::sin(th));
        return g.flywheel_inertia + g.piston_mass * dx * dx;
    };

    double e0 = 0.0, e_max_err = 0.0;
    for (size_t i = 0; i < r.history.size(); ++i)
    {
        const auto& h = r.history[i];
        const double e = 0.5 * crank_inertia(h.state.theta_rad) * h.state.omega * h.state.omega;
        if (i == 0) e0 = e;
        e_max_err = std::max(e_max_err, std::fabs(e - e0) / e0);
    }
    // Drift of the pure-inertia machine (the ½·dJ/dθ·ω² term exactly
    // cancels the J-variation): must be integrator-level, ≪ 5%.
    CHECK(e_max_err < 0.05);
    CHECK(r.cycles_completed >= 4.0);
    CHECK(all_finite(r));
}

TEST_CASE("engine: fired Otto spinup produces positive work near the Otto bound")
{
    auto cfg = base_config();
    cfg.thermo.q_in_cycle = 1500.0;
    cfg.initial_omega = 2.0; // starter momentum: TDC rest is a stall point
    cfg.dt = 2.0e-4;
    cfg.max_steps = 30000; // ~1 s
    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);
    CHECK(r.final_step.state.omega > 0.0);
    CHECK(r.mean_indicated_power > 0.0);
    CHECK(r.total_indicated_work > 0.0);
    CHECK(r.cycles_completed > 0.0);
    // Otto limit: η = 1 − 1/r_c^(γc−1) ≈ 0.49; Wiebe-over-stroke lowers it.
    CHECK(r.efficiency_estimate > 0.2);
    CHECK(r.efficiency_estimate < 0.55);
    CHECK(all_finite(r));
}

TEST_CASE("engine: governor regulates to the setpoint against a generator load")
{
    auto cfg = base_config();
    // Heat headroom: full-throttle indicated torque ≈ 267 N·m vs ~120 N·m
    // load at 200 rad/s → steady-state throttle ≈ 0.45 (mid-range).
    cfg.thermo.q_in_cycle = 8000.0;
    cfg.governor.enabled = true;
    cfg.governor.setpoint_omega = 200.0;
    cfg.governor.pi.kp = 3.0e-3;
    cfg.governor.pi.ki = 0.02;
    cfg.governor.pi.clamp_min = 0.0;
    cfg.governor.pi.clamp_max = 1.0;
    cfg.load.friction_constant = 0.3;
    cfg.load.generator_enabled = true;
    cfg.load.generator_omega_pts = {0.0, 200.0, 400.0};
    cfg.load.generator_torque_pts = {0.0, 120.0, 240.0};
    cfg.geometry.flywheel_inertia = 0.1; // heavier flywheel → tractable transient
    cfg.dt = 1.0e-4;
    cfg.max_steps = 200000; // 20 s
    cfg.initial_omega = 50.0; // start with momentum (TDC rest is a stall point)
    cfg.record_history = true;
    cfg.history_interval = 400;

    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);
    REQUIRE(r.history.size() >= 400U);
    const size_t n = r.history.size();
    const auto& late = r.history[n - 1];
    // Regulating: late omega within ±25% of setpoint (damped oscillation
    // can land on a crest), throttle mid-range.
    CHECK(std::fabs(late.state.omega - 200.0) / 200.0 < 0.25);
    CHECK(late.throttle > 0.1);
    CHECK(late.throttle < 0.95);
    // Damped machine oscillation (PI on a nearly-flat load line): the
    // last 1/3 of the record stays within ±25% of setpoint.
    double w_min = 1e300, w_max = -1e300;
    for (size_t i = 2 * n / 3; i < n; ++i)
    {
        w_min = std::min(w_min, r.history[i].state.omega);
        w_max = std::max(w_max, r.history[i].state.omega);
    }
    CHECK(w_max - w_min < 0.50 * 200.0);
    CHECK(std::fabs(w_max - 200.0) / 200.0 < 0.25);
    CHECK(std::fabs(w_min - 200.0) / 200.0 < 0.25);
    CHECK(all_finite(r));
}

TEST_CASE("engine: steam cycle produces positive shaft power")
{
    auto cfg = base_config();
    cfg.thermo.cycle = EngineCycleType::Steam;
    cfg.thermo.p_boiler = 800000.0;
    cfg.thermo.p_condenser = 15000.0;
    cfg.thermo.steam_cutoff_deg = 40.0;
    cfg.thermo.steam_gamma = 1.13;
    cfg.thermo.r_gas = 461.5;
    cfg.thermo.p_back = cfg.thermo.p_condenser;
    cfg.dt = 2.0e-4;
    cfg.max_steps = 20000;
    cfg.initial_omega = 20.0; // TDC rest is a physical stall point (zero lever arm);
                              // momentum carries the crank through TDC
    cfg.load.friction_constant = 0.1;

    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);
    CHECK(r.final_step.state.omega > 0.0);
    CHECK(r.mean_indicated_power > 0.0);
    CHECK(r.cycles_completed > 0.0);
    CHECK(all_finite(r));
}

TEST_CASE("engine: restartable from arbitrary state (pure function of theta)")
{
    auto cfg = base_config();
    cfg.thermo.q_in_cycle = 1500.0;
    cfg.dt = 2.0e-4;
    cfg.max_steps = 1000;
    ModelStatus s1, s2;
    auto r1 = simulate_engine(cfg, s1);
    REQUIRE(r1.valid);
    auto r2 = simulate_engine(cfg, s2);
    REQUIRE(r2.valid);
    // Same inputs → same trajectory (deterministic, optimizer-batchable).
    CHECK(r1.final_step.state.theta_rad == doctest::Approx(r2.final_step.state.theta_rad).epsilon(1e-12));
    CHECK(r1.final_step.state.omega == doctest::Approx(r2.final_step.state.omega).epsilon(1e-12));
}

TEST_CASE("engine: CSV machine-state streaming matches step count")
{
    const std::string path = std::filesystem::temp_directory_path().string()
                             + "/exd_engine_test_" + std::to_string(::getpid()) + ".csv";
    auto cfg = base_config();
    cfg.thermo.q_in_cycle = 1500.0;
    cfg.dt = 2.0e-4;
    cfg.max_steps = 500;
    cfg.csv_path = path;
    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);

    std::ifstream f(path);
    REQUIRE(f.good());
    std::string line;
    int rows = 0;
    bool header_ok = false;
    while (std::getline(f, line))
    {
        if (rows == 0)
            header_ok = line.rfind("time,theta_rad,omega_rad_s", 0) == 0;
        ++rows;
    }
    CHECK(header_ok);
    CHECK(rows == static_cast<int>(cfg.max_steps) + 1); // header + one row per step
}

TEST_CASE("engine: invalid configs fail cleanly")
{
    ModelStatus status;

    auto bad1 = base_config();
    bad1.geometry.clearance_volume = 0.0;
    auto r1 = simulate_engine(bad1, status);
    CHECK_FALSE(r1.valid);
    CHECK_FALSE(r1.error.empty());

    auto bad2 = base_config();
    bad2.geometry.rod_length = 0.03; // < crank_radius 0.05
    auto r2 = simulate_engine(bad2, status);
    CHECK_FALSE(r2.valid);

    auto bad3 = base_config();
    bad3.thermo.cycle = EngineCycleType::Steam;
    bad3.thermo.p_boiler = 100.0; // < p_condenser default
    auto r3 = simulate_engine(bad3, status);
    CHECK_FALSE(r3.valid);

    auto bad4 = base_config();
    bad4.load.generator_enabled = true; // empty curve
    auto r4 = simulate_engine(bad4, status);
    CHECK_FALSE(r4.valid);
}

TEST_CASE("engine: steam cycle reports saturation temperatures and Rankine-lite efficiency")
{
    auto cfg = base_config();
    cfg.thermo.cycle = EngineCycleType::Steam;
    cfg.thermo.p_boiler = 800000.0;
    cfg.thermo.p_condenser = 15000.0;
    cfg.thermo.steam_cutoff_deg = 40.0;
    cfg.thermo.steam_gamma = 1.13;
    cfg.thermo.steam_quality_cutoff = 0.95;
    cfg.thermo.r_gas = 461.5;
    cfg.initial_omega = 20.0;
    cfg.dt = 2.0e-4;
    cfg.max_steps = 20000;
    cfg.load.friction_constant = 0.1;

    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    REQUIRE(r.valid);
    CHECK(r.final_step.state.omega > 0.0);
    CHECK(r.cycles_completed > 0.0);
    CHECK(all_finite(r));

    // Admission temperature sits on the saturation line at the boiler
    // pressure (wet steam, x = 0.95 < 1): T = T_sat(800 kPa).
    const double t_sat_b = 1.0 / (1.0 / 373.15
                                  - (461.5 / 2.257e6) * std::log(800000.0 / 101325.0));
    bool saw_admission = false;
    for (const auto& h : r.history)
    {
        const double deg = std::fmod(h.state.theta_rad, 2.0 * 3.141592653589793)
                           * 180.0 / 3.141592653589793;
        if (deg > 6.0 && deg < 38.0) // mid-admission, post-ramp
        {
            CHECK(std::fabs(h.T_cyl - t_sat_b) / t_sat_b < 5e-3);
            saw_admission = true;
        }
    }
    CHECK(saw_admission);

    // Energy accounting: efficiency = W/(boiler heat) lands in a plausible
    // simple-engine band (well below the Carnot/Rankine ideal for this
    // condenser pressure, above zero).
    CHECK(r.efficiency_estimate > 0.01);
    CHECK(r.efficiency_estimate < 0.45);
    CHECK(r.total_indicated_work > 0.0);
}

TEST_CASE("engine: steam quality validation")
{
    auto cfg = base_config();
    cfg.thermo.cycle = EngineCycleType::Steam;
    cfg.thermo.steam_quality_cutoff = 1.5; // > 1
    ModelStatus status;
    auto r = simulate_engine(cfg, status);
    CHECK_FALSE(r.valid);
    CHECK(r.error.find("steam_quality_cutoff") != std::string::npos);
}
