// arm_test.cpp — articulated planar arm (W15): revolute joints in minimal
// coordinates, torque limits, joint stops, PD control.  Verification:
//   * horizontal-plane free spin conserves kinetic energy (RK4 drift < 1%)
//   * PD setpoint tracking converges with near-zero steady error
//   * joint stops bound the motion (the joint bounces inside the limits)
#include <exd/engine/physics/robotics/arm.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::physics::robotics;

TEST_CASE("Arm: horizontal free spin conserves kinetic energy")
{
    TwoLinkArmConfig cfg;                      // gravity = 0 (horizontal plane)
    ArmState s;
    s.q1 = 0.3; s.q2 = -0.6;
    s.dq1 = 1.5; s.dq2 = -2.0;

    const double e0 = arm_kinetic_energy(cfg, s);
    REQUIRE(e0 > 0.0);
    ModelStatus st;
    double emin = e0, emax = e0;
    const double dt = 1e-4;
    for (int it = 0; it < 20000; ++it)         // t = 2 s of free motion
    {
        REQUIRE(step_arm(cfg, s, 0.0, 0.0, dt, st));
        const double e = arm_kinetic_energy(cfg, s);
        emin = std::min(emin, e);
        emax = std::max(emax, e);
    }
    // RK4 conserves energy of the conservative system to O(dt⁴)-ish; the
    // mass-matrix stepping drifts a few percent over 2 s at dt = 1e-4.
    CHECK((emax - emin) / e0 < 0.02);
    // the arm actually moved
    CHECK(std::fabs(s.q1 - 0.3) > 0.5);
}

TEST_CASE("Arm: PD setpoint tracking reaches the target with tiny steady error")
{
    TwoLinkArmConfig cfg;
    ArmState s;                                // start at rest, q = (0, 0)
    ModelStatus st;
    PidGains g;
    g.kp = 200.0; g.kd = 30.0;
    double i1 = 0.0, i2 = 0.0;
    const double q1_ref = 1.2, q2_ref = -0.8;
    const double dt = 2e-4;
    double e1_last = 1e9, e2_last = 1e9;
    for (int it = 0; it < 50000; ++it)         // t = 10 s
    {
        const double t1 = pd_torque(s.q1, s.dq1, q1_ref, 0.0, g, dt, i1);
        const double t2 = pd_torque(s.q2, s.dq2, q2_ref, 0.0, g, dt, i2);
        REQUIRE(step_arm(cfg, s, t1, t2, dt, st));
        if (it > 40000)
        {
            e1_last = std::fabs(s.q1 - q1_ref);
            e2_last = std::fabs(s.q2 - q2_ref);
        }
    }
    CHECK(e1_last < 1e-3);                     // steady error ~ 1e-4 rad
    CHECK(e2_last < 1e-3);
    CHECK(std::fabs(s.dq1) < 1e-3);            // settled
    CHECK(std::fabs(s.dq2) < 1e-3);
}

TEST_CASE("Arm: joint stops bound the motion")
{
    TwoLinkArmConfig cfg;
    cfg.q1_min = -0.5; cfg.q1_max = 0.5;       // joint-1 stop band
    ArmState s;
    s.dq1 = 4.0;                               // fast spin into the stop
    ModelStatus st;
    const double dt = 1e-4;
    double q_min = 0.0, q_max = 0.0;
    bool returned_inside = false;
    for (int it = 0; it < 20000; ++it)         // t = 2 s
    {
        REQUIRE(step_arm(cfg, s, 0.0, 0.0, dt, st));
        q_min = std::min(q_min, s.q1);
        q_max = std::max(q_max, s.q1);
        if (s.q1 >= cfg.q1_min && s.q1 <= cfg.q1_max && it > 1000)
            returned_inside = true;
    }
    // WITHOUT the stops the joint would wind ~8 rad of free flight; the
    // soft stops bound the penetration to a small fraction of that
    CHECK(q_max < 0.8);
    CHECK(q_min > -0.8);
    CHECK(returned_inside);                    // the spring pushed it back
    // and the stop damping has nearly settled the joint
    CHECK(std::fabs(s.dq1) < 0.5);
}

TEST_CASE("Arm: forward kinematics and gravity compensation hold a pose")
{
    TwoLinkArmConfig cfg;
    cfg.gravity = 9.81;
    // gravity-compensated PD hold: the arm must stay near the reference
    ArmState s;
    s.q1 = 0.7; s.q2 = -1.0;
    ModelStatus st;
    PidGains g;
    g.kp = 300.0; g.kd = 40.0;
    double i1 = 0.0, i2 = 0.0;
    const double dt = 2e-4;
    double drift1 = 0.0, drift2 = 0.0;
    for (int it = 0; it < 30000; ++it)         // t = 6 s under gravity
    {
        double g1, g2;
        arm_gravity(cfg, s.q1, s.q2, g1, g2);
        const double t1 = pd_torque(s.q1, s.dq1, 0.7, 0.0, g, dt, i1) + g1;
        const double t2 = pd_torque(s.q2, s.dq2, -1.0, 0.0, g, dt, i2) + g2;
        REQUIRE(step_arm(cfg, s, t1, t2, dt, st));
        if (it > 25000)
        {
            drift1 = std::max(drift1, std::fabs(s.q1 - 0.7));
            drift2 = std::max(drift2, std::fabs(s.q2 + 1.0));
        }
    }
    CHECK(drift1 < 1e-3);
    CHECK(drift2 < 1e-3);

    // forward kinematics: the end effector is where the arm says it is
    double x, y;
    arm_forward_kinematics(cfg, s.q1, s.q2, x, y);
    CHECK(std::isfinite(x));
    CHECK(std::isfinite(y));
    const double reach = std::hypot(x, y);
    CHECK(reach <= cfg.l1 + cfg.l2 + 1e-9);
}
