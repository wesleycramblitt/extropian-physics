#include <exd/engine/presets/robotics/arm_trajectory.hpp>

#include <cmath>

namespace exd::engine::presets::robotics {

ArmTrajectoryResult run_arm_trajectory(const ArmTrajectoryConfig& config)
{
    ArmTrajectoryResult result;
    core::ModelStatus& status = result.status;
    exd::engine::physics::robotics::ArmState s;
    s.q1 = config.q1_start;
    s.q2 = config.q2_start;
    const double dt = config.dt;
    const double total = config.move_time + config.settle_time;
    const int steps = static_cast<int>(std::ceil(total / dt));
    if (steps <= 0 || !(dt > 0.0))
    {
        status.ok = false;
        status.error = "arm trajectory: invalid time config";
        return result;
    }

    // smooth reference: quintic-ish ramp (velocity vanishes at both ends)
    auto ramp = [&](double t01) {
        // 0 → 1 over [0, move_time]; 1 afterwards; polynomial with zero
        // endpoint slope (minimum-jerk fragment)
        if (t01 <= 0.0) return 0.0;
        if (t01 >= 1.0) return 1.0;
        const double u = t01;
        return u * u * u * (10.0 - 15.0 * u + 6.0 * u * u);
    };
    auto ref = [&](double t, double q0, double q1) {
        return q0 + (q1 - q0) * ramp(t / config.move_time);
    };

    double i1 = 0.0, i2 = 0.0;
    double t = 0.0;
    result.q1_history.push_back(s.q1);
    result.q2_history.push_back(s.q2);
    result.t_history.push_back(t);
    for (int it = 0; it < steps; ++it)
    {
        t = (it + 1) * dt;
        const double r1 = ref(t, config.q1_start, config.q1_end);
        const double r2 = ref(t, config.q2_start, config.q2_end);
        const double t1 = exd::engine::physics::robotics::pd_torque(
            s.q1, s.dq1, r1, 0.0, config.joint1_pid, dt, i1);
        const double t2 = exd::engine::physics::robotics::pd_torque(
            s.q2, s.dq2, r2, 0.0, config.joint2_pid, dt, i2);
        if (!exd::engine::physics::robotics::step_arm(config.arm, s, t1, t2, dt, status))
        {
            result.ok = false;
            status.error = "arm trajectory: arm step failed";
            return result;
        }
        result.max_tracking_error = std::max(result.max_tracking_error,
                                             std::max(std::fabs(s.q1 - r1),
                                                      std::fabs(s.q2 - r2)));
        result.q1_history.push_back(s.q1);
        result.q2_history.push_back(s.q2);
        result.t_history.push_back(t);
    }
    result.final_error = std::max(std::fabs(s.q1 - config.q1_end),
                                  std::fabs(s.q2 - config.q2_end));
    result.q1_end_reached = s.q1;
    result.q2_end_reached = s.q2;
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::robotics
