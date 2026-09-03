#pragma once

// ---------------------------------------------------------------------------
// Compression-system driver: motor-driven compressor stack feeding a lumped
// plenum through a duct throttle, with optional speed governing.
//
// Full ODE state: {omega, p_plenum, mdot_duct}
//   J*domega_dt = torque_motor - torque_compressor          (shaft balance)
//   dmdot_dt    = (p0_out(stack) - p_plenum) / I            (duct inertance)
//   dp_dt       = (a^2 / V) * (mdot_duct - throttle(p))     (plenum fill)
//
// The compressor row is solved with the stage stack
// (solve_stage_stack) through the same path as the standalone plenum
// module: the compressor characteristic forwarded for the plenum rows is
//   compressor(mdot) = p0_out(omega, mdot) - p_ambient
// so the physics is identical to plenum_derivative. Torque from the stack is
// positive for work INTO the gas (a load opposing +omega on the shaft).
//
// The motor (DcMotorModel, L == 0 quasi-steady allowed) and the PI governor
// are updated exactly once per step -- controller state is mutated at step
// boundaries only, mirroring the engine convention -- and held inside the RK4
// step. The electrical drive torque is positive when assisting.
// ---------------------------------------------------------------------------

#include <exd/engine/physics/control/controller.hpp>
#include <exd/engine/physics/electromagnetics/circuit.hpp>
#include <exd/engine/physics/fluid/lumped/plenum.hpp>
#include <exd/engine/physics/fluid/turbomachinery/stage_stack.hpp>
#include <exd/engine/core/model_status.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::turbomachinery {

struct CompressionSystemConfig
{
    StageStackConfig compressor_stages; // compressor stack (work in)
    lumped::PlenumModelConfig plenum;   // lumped surge cell (lumped/plenum.hpp)
    exd::engine::physics::electromagnetics::DcMotorConfig motor;    // drive (armature; L == 0 quasi-steady)

    double shaft_inertia = 0.5;         // kg*m^2, > 0

    // Throttle: flow = gain * sqrt(max(0, p_plenum - p_ambient)) (kg/s),
    // gain in kg/(s*Pa^0.5). The default 1.5e-3 is a small opening.
    double throttle_gain = 1.5e-3;      // kg/(s*Pa^0.5)

    bool governor_enabled = false;
    exd::engine::physics::control::PiControllerConfig governor_pi; // throttle-gain modulator (speed gov.)
    double governor_setpoint_omega = 900.0;  // rad/s
    double throttle_gain_min = 0.0;          // governor clamps (kg/(s*Pa^0.5))
    double throttle_gain_max = 3.0e-3;

    double dt = 1.0e-4;                    // s, > 0
    uint64_t max_steps = 200000;           // > 0

    double initial_omega = 0.0;            // rad/s
    double initial_p_plenum = 0.0;         // Pa, 0 -> p_ambient
    double initial_mdot_duct = 0.0;        // kg/s

    bool record_history = true;
    uint64_t history_interval = 10;        // >= 1

    std::string csv_path;                  // empty = no CSV machine-state series
};

/// Validate the compression-system config. Fatal problems return false with
/// `error`; non-fatal observations live in `warnings` (mirrors the engine
/// config validator: strict on bad numbers).
bool validate_compression_system_config(const CompressionSystemConfig& config,
                                        std::string& error,
                                        std::vector<std::string>& warnings);

struct CompressionStepResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;
    double t = 0.0;               // time at the END of the step (s)
    double dt_used = 0.0;         // s
    double omega = 0.0;           // rad/s
    double p_plenum = 0.0;        // Pa
    double mdot_duct = 0.0;       // kg/s
    double pressure_ratio = 0.0;  // p0_out / p_ambient at the recorded state
    double torque_compressor = 0.0; // N*m, positive resisting
    double torque_motor = 0.0;       // N*m, positive assisting
    double power_compressor = 0.0;   // W, = torque_compressor * omega
    double power_motor = 0.0;        // W, = torque_motor * omega
    double throttle_gain = 0.0;      // kg/(s*Pa^0.5) used for this step
};

struct CompressionSystemResult
{
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    CompressionStepResult final_step;

    double total_time = 0.0;             // s
    double mean_omega = 0.0;             // rad/s
    double settle_pressure_ratio = 0.0;  // mean over the settled window
    double settle_mdot = 0.0;            // kg/s
    double total_motor_work = 0.0;       // J (motor.electrical torque * omega, rect.)
    double total_compressor_work = 0.0;  // J (compressor shaft torque * omega)
    double efficiency_estimate = 0.0;    // total_compressor_work / total_motor_work
    bool energy_balance_closed = false;  // |motor - comp - dKE|/motor_work < 5%
    std::vector<CompressionStepResult> history;
};

/// Time-march the motor-driven compression system over `max_steps`.
/// History is recorded every `history_interval` steps; when `csv_path` is
/// non-empty a per-step CSV machine-state series is streamed via
/// exd::engine::output::CsvSeriesWriter.
CompressionSystemResult simulate_compression_system(const CompressionSystemConfig& config,
                                                    exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::turbomachinery