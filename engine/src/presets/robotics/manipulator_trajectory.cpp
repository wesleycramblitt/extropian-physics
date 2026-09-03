#include <exd/engine/presets/robotics/manipulator_trajectory.hpp>

#include <cmath>

namespace exd::engine::presets::robotics {

ManipulatorTrajectoryResult run_manipulator_trajectory(
    const ManipulatorTrajectoryConfig& config)
{
    ManipulatorTrajectoryResult result;
    core::ModelStatus& status = result.status;
    const size_t n = config.arm.links.size();
    if (n == 0 ||
        config.q_start.size() != n || config.q_end.size() != n ||
        (!config.joint_pids.empty() && config.joint_pids.size() != n))
    {
        status.ok = false;
        status.error = "manipulator trajectory: config sizes must match the link count";
        return result;
    }
    const double dt = config.dt;
    const double total = config.move_time + config.settle_time;
    const int steps = static_cast<int>(std::ceil(total / dt));
    if (steps <= 0 || !(dt > 0.0))
    {
        status.ok = false;
        status.error = "manipulator trajectory: invalid time config";
        return result;
    }

    exd::engine::physics::robotics::ManipulatorState s;
    s.q = config.q_start;
    s.dq.assign(n, 0.0);

    // smooth reference: polynomial ramp with zero endpoint slope
    auto ramp = [&](double t01) {
        if (t01 <= 0.0) return 0.0;
        if (t01 >= 1.0) return 1.0;
        const double u = t01;
        return u * u * u * (10.0 - 15.0 * u + 6.0 * u * u);
    };

    std::vector<double> integral(n, 0.0), taus(n), ref(n);
    exd::engine::physics::robotics::PidGains default_pid;
    result.q_history.push_back(s.q);
    for (int it = 0; it < steps; ++it)
    {
        const double t = (it + 1) * dt;
        const double f = ramp(t / config.move_time);
        for (size_t j = 0; j < n; ++j)
        {
            const double r = config.q_start[j] + (config.q_end[j] - config.q_start[j]) * f;
            ref[j] = r;
            const auto& g = config.joint_pids.empty() ? default_pid : config.joint_pids[j];
            taus[j] = exd::engine::physics::robotics::pd_torque(
                s.q[j], s.dq[j], r, 0.0, g, dt, integral[j]);
        }
        if (!exd::engine::physics::robotics::step_manipulator(config.arm, s, taus, dt, status))
        {
            result.ok = false;
            status.error = "manipulator trajectory: step failed";
            return result;
        }
        for (size_t j = 0; j < n; ++j)
            result.max_tracking_error = std::max(result.max_tracking_error,
                                                 std::fabs(s.q[j] - ref[j]));
        result.q_history.push_back(s.q);
    }
    result.final_error = 0.0;
    result.q_end_reached = s.q;
    for (size_t j = 0; j < n; ++j)
        result.final_error = std::max(result.final_error, std::fabs(s.q[j] - config.q_end[j]));
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::robotics
