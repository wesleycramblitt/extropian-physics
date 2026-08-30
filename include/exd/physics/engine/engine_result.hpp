#pragma once

// ─────────────────────────────────────────────────────
// Engine result types (value types, no exceptions).
// ─────────────────────────────────────────────────────

#include <exd/physics/mechanics/status.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::engine {

/// Minimal persisted machine state. Crank angle theta is
/// measured from TDC (start of the power stroke); the cycle
/// phase is derived from theta mod the cycle period, so a
/// simulation is restartable from any state.
struct EngineState
{
    double theta_rad = 0.0;  // crank angle (rad)
    double omega = 0.0;      // angular velocity (rad/s)
    uint64_t cycles = 0;     // completed power cycles
};

/// One step of the engine, evaluated AFTER the integrator
/// advanced the state. All kinematic/thermo quantities are
/// pure functions of the final crank angle.
struct EngineStepResult
{
    bool ok = false;
    mechanics::ModelStatus status;
    double t = 0.0;          // time at the END of the step
    double dt_used = 0.0;

    EngineState state;

    // Mechanism (animation-relevant)
    double piston_x = 0.0;   // piston from crank axis (m); TDC = l + r
    double piston_v = 0.0;   // piston velocity (m/s)

    // Cylinder
    double p_cyl = 0.0;      // Pa
    double T_cyl = 0.0;      // K (power/compression phases physical;
                             // valve-open phases nominal)
    double V_cyl = 0.0;      // m³

    // Torques / power
    double gas_force = 0.0;      // N on piston (+x is toward the head)
    double indicated_moment = 0.0; // N·m from gas on the crank
    double load_moment = 0.0;    // N·m opposing (friction + generator)
    double power = 0.0;          // net shaft power = (M_ind − M_load)·ω (W)

    double throttle = 1.0;       // governor heat-release fraction [0,1]
};

/// Full simulation result. History entries are recorded at
/// `history_interval` steps (like the turbine app).
struct EngineSimResult
{
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    EngineStepResult final_step;
    double total_time = 0.0;          // s
    double total_indicated_work = 0.0;// J (Σ M_ind·ω·dt)
    double mean_indicated_power = 0.0;// W
    double mean_omega = 0.0;          // rad/s
    double cycles_completed = 0.0;
    double mean_throttle = 1.0;          // time-mean governor heat fraction;
                                         // efficiency uses throttled heat
    double efficiency_estimate = 0.0;    // W_ind / heat released: Otto →
                                         // q_in rollup; steam → boiler heat
                                         // m·(h_g − h_f) per cycle (Rankine-lite)
    std::vector<EngineStepResult> history;
};

} // namespace exd::physics::engine
