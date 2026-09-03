// robotics_trajectory_test.cpp — N-link manipulator trajectory preset (W15):
//   * any link count runs the same preset (2 and 3 links)
//   * the arm follows a smooth reference between two poses
//   * torque limits are respected (finite, stable runs)
#include <exd/engine/presets/robotics/manipulator_trajectory.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::presets::robotics;

TEST_CASE("Manipulator trajectory preset: 2-link pose-to-pose move")
{
    ManipulatorTrajectoryConfig cfg;
    cfg.arm.links = {{0.5, 1.0}, {0.4, 1.0}};
    cfg.arm.torque_max = {5.0, 2.0};           // limited actuators
    cfg.joint_pids = {{40.0, 12.0, 0.0, 10.0}, {15.0, 6.0, 0.0, 10.0}};
    cfg.q_start = {0.0, 0.0};
    cfg.q_end = {1.2, -0.9};
    cfg.move_time = 3.0;
    cfg.settle_time = 2.0;
    cfg.dt = 1e-3;

    const auto r = run_manipulator_trajectory(cfg);
    REQUIRE(r.ok);
    REQUIRE(r.q_end_reached.size() == 2);
    CHECK(std::fabs(r.q_end_reached[0] - 1.2) < 0.02);
    CHECK(std::fabs(r.q_end_reached[1] + 0.9) < 0.02);
    CHECK(r.final_error < 0.02);
    CHECK(r.max_tracking_error < 0.5);
    CHECK(r.q_end_reached[0] > 0.5);
    CHECK(r.q_end_reached[1] < -0.3);
}

TEST_CASE("Manipulator trajectory preset: 3-link move runs the same preset")
{
    ManipulatorTrajectoryConfig cfg;
    cfg.arm.links = {{0.4, 1.0}, {0.3, 0.8}, {0.2, 0.5}};
    cfg.arm.torque_max = {4.0, 3.0, 1.0};
    cfg.joint_pids = {{60.0, 15.0, 0.0, 10.0},
                      {40.0, 10.0, 0.0, 10.0},
                      {20.0, 6.0, 0.0, 10.0}};
    cfg.q_start = {0.2, -0.2, 0.3};
    cfg.q_end = {1.0, -0.8, 0.5};
    cfg.move_time = 3.0;
    cfg.settle_time = 2.0;
    cfg.dt = 1e-3;

    const auto r = run_manipulator_trajectory(cfg);
    REQUIRE(r.ok);
    REQUIRE(r.q_end_reached.size() == 3);
    CHECK(std::fabs(r.q_end_reached[0] - 1.0) < 0.03);
    CHECK(std::fabs(r.q_end_reached[1] + 0.8) < 0.03);
    CHECK(r.final_error < 0.03);
}
