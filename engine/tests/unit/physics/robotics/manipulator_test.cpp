// manipulator_test.cpp — N-link serial manipulator (W15): the GENERAL class
// (2 and 3 links here; any link count runs the same solver).  Verification:
//   * horizontal-plane free spin conserves kinetic energy
//   * PD setpoint tracking converges with near-zero steady error
//   * joint stops bound the motion
//   * gravity compensation holds a pose
#include <exd/engine/physics/robotics/manipulator.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::physics::robotics;

namespace {
SerialManipulatorConfig two_link()
{
    SerialManipulatorConfig cfg;
    cfg.links = {{0.5, 1.0}, {0.4, 1.0}};
    return cfg;
}
} // namespace

TEST_CASE("Manipulator: horizontal free spin conserves kinetic energy")
{
    auto cfg = two_link();                     // gravity = 0
    ManipulatorState s;
    s.q = {0.3, -0.6};
    s.dq = {1.5, -2.0};

    const double e0 = kinetic_energy(cfg, s);
    REQUIRE(e0 > 0.0);
    ModelStatus st;
    double emin = e0, emax = e0;
    const double dt = 1e-4;
    for (int it = 0; it < 20000; ++it)         // t = 2 s
    {
        REQUIRE(step_manipulator(cfg, s, {0.0, 0.0}, dt, st));
        const double e = kinetic_energy(cfg, s);
        emin = std::min(emin, e);
        emax = std::max(emax, e);
    }
    CHECK((emax - emin) / e0 < 0.02);
    CHECK(std::fabs(s.q[0] - 0.3) > 0.5);      // it actually moved
}

TEST_CASE("Manipulator: PD setpoint tracking reaches the target")
{
    auto cfg = two_link();
    ManipulatorState s;                        // at rest, q = (0, 0)
    s.q = {0.0, 0.0};
    s.dq = {0.0, 0.0};
    ModelStatus st;
    PidGains g;
    g.kp = 200.0; g.kd = 30.0;
    double i1 = 0.0, i2 = 0.0;
    const double q1_ref = 1.2, q2_ref = -0.8;
    const double dt = 2e-4;
    double e1_last = 1e9, e2_last = 1e9;
    for (int it = 0; it < 50000; ++it)
    {
        const double t1 = pd_torque(s.q[0], s.dq[0], q1_ref, 0.0, g, dt, i1);
        const double t2 = pd_torque(s.q[1], s.dq[1], q2_ref, 0.0, g, dt, i2);
        REQUIRE(step_manipulator(cfg, s, {t1, t2}, dt, st));
        if (it > 40000)
        {
            e1_last = std::fabs(s.q[0] - q1_ref);
            e2_last = std::fabs(s.q[1] - q2_ref);
        }
    }
    CHECK(e1_last < 1e-3);
    CHECK(e2_last < 1e-3);
    CHECK(std::fabs(s.dq[0]) < 1e-3);
    CHECK(std::fabs(s.dq[1]) < 1e-3);
}

TEST_CASE("Manipulator: joint stops bound the motion")
{
    auto cfg = two_link();
    cfg.q_min = {-0.5, -1.0};
    cfg.q_max = {0.5, 1.0};
    ManipulatorState s;
    s.q = {0.0, 0.0};
    s.dq = {4.0, 0.0};                         // fast spin into the stop
    ModelStatus st;
    const double dt = 1e-4;
    double q_min = 0.0, q_max = 0.0;
    bool returned_inside = false;
    for (int it = 0; it < 20000; ++it)         // t = 2 s
    {
        REQUIRE(step_manipulator(cfg, s, {0.0, 0.0}, dt, st));
        q_min = std::min(q_min, s.q[0]);
        q_max = std::max(q_max, s.q[0]);
        if (s.q[0] >= cfg.q_min[0] && s.q[0] <= cfg.q_max[0] && it > 1000)
            returned_inside = true;
    }
    CHECK(q_max < 0.8);                        // bounded vs ~8 rad free flight
    CHECK(q_min > -0.8);
    CHECK(returned_inside);
    CHECK(std::fabs(s.dq[0]) < 0.5);           // the stop damping settles it
}

TEST_CASE("Manipulator: gravity compensation holds a pose; 3-link general")
{
    // 3-link configuration (the general class, n = 3)
    SerialManipulatorConfig cfg;
    cfg.links = {{0.4, 1.0}, {0.3, 0.8}, {0.2, 0.5}};
    cfg.gravity = 9.81;
    ManipulatorState s;
    s.q = {0.7, -0.5, 0.9};
    s.dq = {0.0, 0.0, 0.0};
    ModelStatus st;
    PidGains g;
    g.kp = 300.0; g.kd = 40.0;
    std::vector<double> integral(3, 0.0), taus(3);
    double drift = 0.0;
    const double dt = 2e-4;
    for (int it = 0; it < 30000; ++it)         // t = 6 s
    {
        std::vector<double> gg;
        gravity(cfg, s.q, gg);
        for (size_t j = 0; j < 3; ++j)
            taus[j] = pd_torque(s.q[j], s.dq[j], s.q[j], 0.0, g, dt, integral[j]) + gg[j];
        REQUIRE(step_manipulator(cfg, s, taus, dt, st));
        if (it > 25000)
            for (size_t j = 0; j < 3; ++j)
                drift = std::max(drift, std::fabs(s.dq[j]));
    }
    CHECK(drift < 1e-3);                       // the pose holds under gravity

    double x, y;
    forward_kinematics(cfg, s.q, x, y);
    CHECK(std::isfinite(x));
    CHECK(std::isfinite(y));
    const double reach = std::hypot(x, y);
    CHECK(reach <= 0.9 + 1e-9);                // within the summed link length
}
