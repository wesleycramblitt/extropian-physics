#pragma once

// ─────────────────────────────────────────────────────
// Engine configuration (single-cylinder slider-crank).
//
// Cycle models (first pass, engineering-grade):
//   Otto   — 4-stroke: polytropic compression/expansion,
//            Wiebe heat release on the power stroke,
//            smooth sin² valve-window ramps.
//   Steam  — single-acting placeholder: constant-pressure
//            admission to cutoff, polytropic expansion
//            (n = steam_gamma, 1.13 = saturated-steam
//            placeholder — NOT a real steam table),
//            exhaust to condenser.
// Loads (friction + generator curve) are pure T(ω)
// functions of the config; no allocation in the hot loop.
// ─────────────────────────────────────────────────────

#include <exd/engine/physics/control/controller.hpp>
#include <exd/engine/numerics/integrators.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::presets::engine {

enum class EngineCycleType : uint8_t
{
    Otto,   // 4-stroke, Wiebe heat release
    Steam,  // single-acting, admission + polytropic expansion
};

struct EngineGeometryConfig
{
    double crank_radius = 0.05;    // m, > 0
    double rod_length = 0.20;      // m, > crank_radius
    double bore = 0.086;           // m, > 0
    double clearance_volume = 1.0e-4; // m³, > 0 (TDC volume; keeps pV^γ finite)
    double piston_mass = 0.5;      // kg, ≥ 0
    double flywheel_inertia = 0.02;// kg·m², > 0
};

struct EngineThermoConfig
{
    EngineCycleType cycle = EngineCycleType::Otto;

    // Shared
    double r_gas = 287.0;          // J/(kg·K) (air; steam ≈ 461.5)
    double p_intake = 101325.0;    // Pa, > 0 (Otto: intake manifold; steam: boiler ref)
    double T_intake = 300.0;       // K, > 0 (admission/trapped-gas temperature)
    double p_exhaust = 101325.0;   // Pa, > 0 (Otto back pressure)
    double p_back = 101325.0;      // Pa crankcase/condenser back pressure
                                   // (Otto: set to p_exhaust default; steam: set to
                                   //  p_condenser default — see below)

    // Otto
    double gamma_compression = 1.35; // per-phase polytrope, > 1
    double gamma_expansion = 1.35;   // > 1. NOTE: differing from
                                     // gamma_compression injects net
                                     // work per cycle with zero heat
                                     // (a heat-transfer STAND-IN; keep
                                     // equal for conservative motored
                                     // behavior — validated with warning)
    double q_in_cycle = 1500.0;      // J heat release per cycle, ≥ 0
    double wiebe_ignition_deg = 0.0; // °ATDC start of combustion
    double wiebe_burn_duration_deg = 60.0; // °, > 0
    double wiebe_a = 6.908;          // Wiebe efficiency factor
    double wiebe_m = 2.0;            // Wiebe form factor

    // Steam (placeholder model)
    double p_boiler = 800000.0;      // Pa admission pressure, > p_condenser
    double p_condenser = 15000.0;    // Pa exhaust pressure, > 0
    double steam_cutoff_deg = 40.0;  // ° admission cutoff, in (0, 180)
    double steam_gamma = 1.13;       // wet-steam polytrope exponent (saturated
                                     // steam ≈ 1.13; engineering model)
    double steam_quality_cutoff = 0.95; // dryness fraction at cutoff, in (0,1]
                                     // (0.95: slightly wet admission)
};

struct EngineLoadConfig
{
    double friction_constant = 0.3;  // N·m, ≥ 0
    double friction_viscous = 0.0;   // N·m·s/rad, ≥ 0
    bool generator_enabled = false;  // curve load below
    std::vector<double> generator_omega_pts;  // rad/s, strictly increasing
    std::vector<double> generator_torque_pts; // N·m (opposing)
};

struct EngineGovernorConfig
{
    bool enabled = false;
    double setpoint_omega = 200.0;       // rad/s
    exd::engine::physics::control::PiControllerConfig pi;      // defaults: kp/ki tuned for the simple engine
    double throttle_min = 0.0;           // heat-release fraction clamp
    double throttle_max = 1.0;
};

struct EngineConfig
{
    EngineGeometryConfig geometry;
    EngineThermoConfig thermo;
    EngineLoadConfig load;
    EngineGovernorConfig governor;

    exd::engine::numerics::IntegratorConfig integration; // default RK4

    double dt = 5.0e-4;              // s, > 0
    uint64_t max_steps = 20000;      // > 0

    double initial_theta_rad = 0.0;  // 0 = TDC start of power stroke
    double initial_omega = 0.0;      // rad/s

    bool record_history = true;
    uint64_t history_interval = 1;   // ≥ 1

    std::string csv_path;            // empty = no CSV machine-state output
};

/// Validate the config. Errors are fatal (invalid result);
/// warnings are non-fatal.
bool validate_engine_config(const EngineConfig& config,
                            std::string& error,
                            std::vector<std::string>& warnings);

} // namespace exd::engine::presets::engine
