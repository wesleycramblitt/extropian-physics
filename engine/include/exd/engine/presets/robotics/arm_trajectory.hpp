#pragma once

// ─────────────────────────────────────────────────────
// Robotic-arm trajectory preset (robotic-arm use case).
//
// Assembly: the 2R articulated arm (physics/robotics)
// with per-joint PD control along a smooth reference —
// a minimum-jerk style ramp between two poses.  All the
// machinery (joints, mass coupling, Coriolis, torque
// limits, stops, PD) lives in the physics module; this
// preset only wires a reference generator + controller
// loop and reports the tracking error.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/robotics/arm.hpp>

#include <vector>

namespace exd::engine::presets::robotics {

struct ArmTrajectoryConfig
{
    exd::engine::physics::robotics::TwoLinkArmConfig arm;
    exd::engine::physics::robotics::PidGains joint1_pid;
    exd::engine::physics::robotics::PidGains joint2_pid;
    double q1_start = 0.0, q2_start = 0.0;
    double q1_end = 1.0, q2_end = -0.6;
    double move_time = 2.0;            // s (smooth ramp duration)
    double settle_time = 1.0;          // s of dwell at the target
    double dt = 1e-3;
};

struct ArmTrajectoryResult
{
    bool ok = false;
    core::ModelStatus status;
    std::vector<double> q1_history, q2_history;
    std::vector<double> t_history;
    double max_tracking_error = 0.0;   // rad, over the whole run
    double final_error = 0.0;          // rad at the end of the dwell
    double q1_end_reached = 0.0, q2_end_reached = 0.0;
};

/// Run the arm along the reference trajectory under PD control.
ArmTrajectoryResult run_arm_trajectory(const ArmTrajectoryConfig& config);

} // namespace exd::engine::presets::robotics
