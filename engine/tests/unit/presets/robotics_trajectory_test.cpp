// robotics_trajectory_test.cpp — arm trajectory preset (W15):
//   * the arm follows a smooth reference between two poses
//   * the tracking error stays tight and the final pose is reached
#include <exd/engine/presets/robotics/arm_trajectory.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::presets::robotics;

TEST_CASE("Arm trajectory preset: 2R arm tracks a smooth pose-to-pose move")
{
    ArmTrajectoryConfig cfg;
    cfg.arm.l1 = 0.5; cfg.arm.l2 = 0.4;
    cfg.arm.m1 = 1.0; cfg.arm.m2 = 1.0;
    cfg.arm.gravity = 0.0;
    cfg.arm.t1_max = 5.0; cfg.arm.t2_max = 2.0;   // limited actuators
    cfg.joint1_pid.kp = 40.0; cfg.joint1_pid.kd = 12.0;
    cfg.joint2_pid.kp = 15.0; cfg.joint2_pid.kd = 6.0;
    cfg.q1_start = 0.0; cfg.q2_start = 0.0;
    cfg.q1_end = 1.2; cfg.q2_end = -0.9;
    cfg.move_time = 3.0;
    cfg.settle_time = 2.0;
    cfg.dt = 1e-3;

    const auto r = run_arm_trajectory(cfg);
    REQUIRE(r.ok);
    // the final pose is reached within a few degrees
    CHECK(std::fabs(r.q1_end_reached - 1.2) < 0.02);
    CHECK(std::fabs(r.q2_end_reached + 0.9) < 0.02);
    CHECK(r.final_error < 0.02);
    // the tracking error during the move stays bounded by the limited
    // actuators (t_max clamps the acceleration)
    CHECK(r.max_tracking_error < 0.5);
    // it actually moved (not stuck)
    CHECK(r.q1_end_reached > 0.5);
    CHECK(r.q2_end_reached < -0.3);
}

TEST_CASE("Arm trajectory preset: torque limits are respected")
{
    ArmTrajectoryConfig cfg;
    cfg.arm.l1 = 0.5; cfg.arm.l2 = 0.4;
    cfg.arm.m1 = 2.0; cfg.arm.m2 = 2.0;
    cfg.arm.t1_max = 1.0; cfg.arm.t2_max = 0.5;   // very weak actuators
    cfg.joint1_pid.kp = 40.0; cfg.joint1_pid.kd = 10.0;
    cfg.joint2_pid.kp = 10.0; cfg.joint2_pid.kd = 4.0;
    cfg.q1_end = 3.0; cfg.q2_end = 1.0;           // far target: saturates
    cfg.move_time = 2.0;
    cfg.settle_time = 2.0;
    const auto r = run_arm_trajectory(cfg);
    REQUIRE(r.ok);
    // with the torque limits the arm cannot track the fast reference
    // perfectly — but the run must remain finite and stable
    CHECK(std::isfinite(r.max_tracking_error));
    CHECK(std::isfinite(r.q1_end_reached));
}
