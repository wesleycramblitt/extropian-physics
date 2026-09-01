// stage_stack.cpp
// Multi-stage propagation: total state (p0, T0, c_theta) is carried stage to
// stage at fixed omega/mdot; torques and specific work sum on one shaft.

#include <exd/physics/fluid/turbomachinery/stage_stack.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace exd::physics::fluid::turbomachinery
{

bool validate_stage_stack_config(const StageStackConfig& config,
                                 std::string& error,
                                 std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.stages.empty())
    {
        error = "stage_stack: at least one stage is required";
        return false;
    }
    for (size_t i = 0; i < config.stages.size(); ++i)
    {
        std::string stage_error;
        std::vector<std::string> stage_warnings;
        if (!validate_stage_config(config.stages[i], stage_error, stage_warnings))
        {
            error = "stage_stack: stage " + std::to_string(i) + ": " + stage_error;
            return false;
        }
        warnings.insert(warnings.end(), stage_warnings.begin(), stage_warnings.end());
    }
    return true;
}

StageStackResult solve_stage_stack(const StageStackConfig& config,
                                   const StageInlet& inlet,
                                   double omega,
                                   double mdot,
                                   const thermo::IEos& eos,
                                   exd::physics::ModelStatus& status)
{
    StageStackResult result;
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    std::string verror;
    std::vector<std::string> vwarnings;
    if (!validate_stage_stack_config(config, verror, vwarnings))
    {
        status.ok = false;
        status.error = verror;
        result.status = status;
        return result;
    }
    status.warnings.insert(status.warnings.end(), vwarnings.begin(), vwarnings.end());

    StageInlet current = inlet;
    result.total_delta_h0 = 0.0;
    result.total_torque = 0.0;
    result.total_power = 0.0;
    result.mach_rel_max = 0.0;
    result.total_pi = 1.0;
    result.p0_out = inlet.p0;
    result.T0_out = inlet.T0;
    result.c_theta_out = inlet.c_theta;
    result.per_stage.reserve(config.stages.size());

    for (const StageConfig& stage_config : config.stages)
    {
        ModelStatus stage_status;
        StageResult stage =
            solve_stage(stage_config, current, omega, mdot, eos, stage_status);
        result.per_stage.push_back(stage);
        status.warnings.insert(status.warnings.end(),
                               stage_status.warnings.begin(),
                               stage_status.warnings.end());
        if (!stage.ok)
        {
            status.ok = false;
            status.error = stage.status.error;
            result.status = status;
            return result;
        }

        result.total_delta_h0 += stage.delta_h0;
        result.total_torque += stage.torque;
        result.total_power += stage.power;
        result.mach_rel_max = std::max(result.mach_rel_max, stage.mach_rel_le);

        current.p0 = stage.p0_out;
        current.T0 = stage.T0_out;
        current.c_theta = stage.c_theta_out;
        result.p0_out = current.p0;
        result.T0_out = current.T0;
        result.c_theta_out = current.c_theta;
    }

    result.total_pi = result.p0_out / inlet.p0;
    result.ok = status.ok;
    result.status = status;
    return result;
}

} // namespace exd::physics::fluid::turbomachinery