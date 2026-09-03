// pi_controller_test.cpp
// Unit tests for the PI controller with output clamping and anti-windup.

#include <exd/engine/physics/control/controller.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>

using namespace exd::engine::physics::control;

TEST_CASE("PI controller: proportional-only response equals kp*err")
{
    PiControllerConfig config;
    config.kp = 2.0;
    config.ki = 0.0; // pure proportional

    auto controller = make_pi_controller(config);
    REQUIRE(controller);
    CHECK(controller->name() == "pi");

    exd::engine::core::ModelStatus status;

    CHECK(controller->update(10.0, 7.0, 0.01, status) == doctest::Approx(6.0).epsilon(1e-12));
    CHECK(controller->update(10.0, 10.0, 0.01, status) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(controller->update(0.0, 5.0, 0.01, status) == doctest::Approx(-10.0).epsilon(1e-12));
    CHECK(status.ok);
}

TEST_CASE("PI controller: tracks a setpoint on a first-order plant")
{
    PiControllerConfig config;
    config.kp = 10.0;
    config.ki = 50.0;
    config.clamp_min = -100.0;
    config.clamp_max = 100.0;
    config.anti_windup = true;

    auto controller = make_pi_controller(config);
    REQUIRE(controller);

    exd::engine::core::ModelStatus status;
    constexpr double tau = 0.5;
    constexpr double dt = 0.001;
    constexpr double setpoint = 1.0;

    // Plant: dx/dt = (u - x)/tau, forward Euler over 3000 steps of dt.
    double x = 0.0;
    bool last_phase_stays = true;
    for (int i = 0; i < 3000; ++i)
    {
        const double u = controller->update(setpoint, x, dt, status);
        REQUIRE(status.ok);
        x += dt * (u - x) / tau;
        if (i >= 1500 && std::fabs(x - setpoint) > 0.02 * setpoint)
            last_phase_stays = false;
    }

    CHECK(x == doctest::Approx(setpoint).epsilon(0.02));
    CHECK(last_phase_stays); // once tracking, the state stays within 2%
}

TEST_CASE("PI controller: anti-windup bounds the integral while the output is saturated")
{
    PiControllerConfig config;
    config.kp = 0.0; // integral-only
    config.ki = 1.0;
    config.clamp_min = -1.0;
    config.clamp_max = 1.0;
    config.anti_windup = true;

    auto controller = make_pi_controller(config);
    REQUIRE(controller);

    exd::engine::core::ModelStatus status;
    constexpr double dt = 0.001;

    // Saturation phase: a constant error 5 (setpoint 5, measurement 0) keeps
    // the integral-only controller pegged at clamp_max. Anti-windup stops the
    // integral from racing upward while the output is clamped (without it the
    // integral would reach ki·err·time = 50 over these 10000 steps); the
    // returned output never grows past the clamp.
    double u = 0.0;
    double u_max = 0.0;
    for (int i = 0; i < 10000; ++i)
    {
        u = controller->update(5.0, 0.0, dt, status);
        REQUIRE(status.ok);
        u_max = std::max(u_max, u);
        if (i >= 2000) // comfortably after saturation (~step 200)
            CHECK(u == doctest::Approx(config.clamp_max).epsilon(1e-12));
    }
    CHECK(u == doctest::Approx(config.clamp_max).epsilon(1e-12));
    CHECK(u_max <= config.clamp_max + 1e-12); // never drives past the clamp

    // Recovery phase: reverse the error. With the integral bounded (~1.0) by
    // anti-windup, the output comes off the clamp immediately and the integral
    // unwinds to ~0 within a few hundred steps (without anti-windup the ~50
    // integral would take tens of thousands of steps to unwind).
    bool came_off_clamp = false;
    int unwind_steps = -1;
    for (int i = 0; i < 5000; ++i)
    {
        const double w = controller->update(5.0, 10.0, dt, status); // err = -5
        REQUIRE(status.ok);
        if (w < config.clamp_max)
            came_off_clamp = true;
        if (std::fabs(w) < 0.05)
        {
            unwind_steps = i + 1;
            break;
        }
    }
    CHECK(came_off_clamp);
    CHECK(unwind_steps > 0);
    CHECK(unwind_steps < 5000);

    // With the integral unwound near zero, removing the disturbance
    // (measurement back at the setpoint) leaves the output near zero.
    const double u_hold = controller->update(5.0, 5.0, dt, status);
    CHECK(std::fabs(u_hold) < 0.1);
    CHECK(status.ok);
}

TEST_CASE("PI controller: reset() clears the accumulated integral state")
{
    PiControllerConfig config;
    config.kp = 2.0;
    config.ki = 10.0;

    auto controller = make_pi_controller(config);
    REQUIRE(controller);

    exd::engine::core::ModelStatus status;
    constexpr double dt = 0.01;

    // Build up a large integral with a constant error of 10.
    for (int i = 0; i < 100; ++i)
        controller->update(10.0, 0.0, dt, status);
    REQUIRE(status.ok);

    const double u_built = controller->update(10.0, 0.0, dt, status);
    CHECK(u_built > 50.0); // dominated by the accumulated integral

    controller->reset();

    // After reset the controller responds like a freshly constructed one:
    // pure proportional plus a single integral step.
    const double u_after = controller->update(10.0, 0.0, dt, status);
    const double expected = 2.0 * 10.0 + 10.0 * 10.0 * dt;
    CHECK(u_after == doctest::Approx(expected).epsilon(1e-9));
    CHECK(status.ok);
}

TEST_CASE("PI controller: rejects non-positive dt")
{
    PiControllerConfig config;
    config.kp = 2.0;
    config.ki = 1.0;

    auto controller = make_pi_controller(config);
    REQUIRE(controller);

    exd::engine::core::ModelStatus status;

    SUBCASE("dt equal to zero")
    {
        CHECK(controller->update(1.0, 0.0, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK(status.error == "PI controller: dt must be positive");
    }

    SUBCASE("negative dt")
    {
        status = {};
        CHECK(controller->update(1.0, 0.0, -0.5, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("a rejected step does not advance the integral")
    {
        status = {};
        const double first = controller->update(10.0, 8.0, 0.01, status);
        CHECK(first == doctest::Approx(4.02).epsilon(1e-9));
        CHECK(status.ok);

        status = {};
        CHECK(controller->update(10.0, 8.0, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);

        status = {}; // caller-owned status: reset before the next call
        const double second = controller->update(10.0, 8.0, 0.01, status);
        CHECK(second == doctest::Approx(4.04).epsilon(1e-9));
        CHECK(status.ok);
    }
}