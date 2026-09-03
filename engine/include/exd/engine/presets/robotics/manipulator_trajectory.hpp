#pragma once

// ─────────────────────────────────────────────────────
// Manipulator trajectory preset (robotic-arm CLASS of
// problems).
//
// Assembly: the N-link serial manipulator (physics/
// robotics) with per-joint PD control along a smooth
// (minimum-jerk style) reference between two poses.
// The CONFIG is the user's data — any serial arm they
// import (link counts, geometry, inertia, joint limits,
// actuator ratings, controller gains) runs the same
// solver.  The physics all lives in the module; this
// preset only wires the reference generator + control
// loop and reports the tracking error.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/robotics/manipulator.hpp>

#include <vector>

namespace exd::engine::presets::robotics {

struct ManipulatorTrajectoryConfig
{
    exd::engine::physics::robotics::SerialManipulatorConfig arm;
    std::vector<exd::engine::physics::robotics::PidGains> joint_pids;
    std::vector<double> q_start;         // per joint (rad)
    std::vector<double> q_end;
    double move_time = 2.0;              // s (smooth ramp)
    double settle_time = 1.0;            // s of dwell at the target
    double dt = 1e-3;
};

struct ManipulatorTrajectoryResult
{
    bool ok = false;
    core::ModelStatus status;
    std::vector<std::vector<double>> q_history;   // per-step joint angles
    double max_tracking_error = 0.0;              // rad
    double final_error = 0.0;                     // rad
    std::vector<double> q_end_reached;
};

/// Run the arm along the reference trajectory under PD control.
ManipulatorTrajectoryResult run_manipulator_trajectory(
    const ManipulatorTrajectoryConfig& config);

} // namespace exd::engine::presets::robotics
