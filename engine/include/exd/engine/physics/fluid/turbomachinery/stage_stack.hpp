#pragma once

// Stage-stack (multi-stage axial compressor/turbine) propagation.
// Stages run in series at fixed mass flow and shaft speed; total state
// (p0, T0, c_theta) is carried stage to stage. Torques sum on one shaft.

#include <exd/engine/physics/fluid/turbomachinery/stage.hpp>

#include <string>
#include <vector>

namespace exd::engine::physics::fluid::turbomachinery
{

struct StageStackConfig
{
    std::vector<StageConfig> stages; // >= 1 stage
};

/// Validate a stack config: non-empty, and every stage validates.
/// Per-stage warnings (envelope guards) are collected.
bool validate_stage_stack_config(const StageStackConfig& config,
                                 std::string& error,
                                 std::vector<std::string>& warnings);

struct StageStackResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;

    double p0_out = 0.0;     // total pressure after the last stage, Pa
    double T0_out = 0.0;     // total temperature after the last stage, K
    double c_theta_out = 0.0;// absolute swirl after the last stage, m/s
    double total_pi = 0.0;   // p0_out / p0_in over the whole stack
    double total_delta_h0 = 0.0;    // sum of stage specific work, J/kg
    double total_torque = 0.0;      // sum of stage torques (same shaft), N*m
    double total_power = 0.0;       // sum of stage power, W
    double mach_rel_max = 0.0;      // max relative LE Mach over stages
    std::vector<StageResult> per_stage;
};

/// Solve the stack at constant omega/mdot, propagating (p0, T0, c_theta).
/// Stage warnings accumulate on `status`; the first failing stage returns
/// ok = false with its error.
StageStackResult solve_stage_stack(const StageStackConfig& config,
                                   const StageInlet& inlet,
                                   double omega,
                                   double mdot,
                                   const exd::engine::physics::thermo::IEos& eos,
                                   exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::turbomachinery