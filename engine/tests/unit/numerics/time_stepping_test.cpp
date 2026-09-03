// time_stepping_test.cpp
// Unit tests for the promoted public time-stepping infrastructure
// (TimeStepper, ConvergenceMonitor).

#include <exd/engine/numerics/time_stepping.hpp>

#include <doctest/doctest.h>

using namespace exd::engine::numerics;

TEST_CASE("TimeStepper: sensible defaults")
{
    TimeStepper ts;
    CHECK(ts.t == 0.0);
    CHECK(ts.dt == doctest::Approx(0.001));
    CHECK(ts.dt_min == doctest::Approx(1e-8));
    CHECK(ts.dt_max == doctest::Approx(1.0));
    CHECK(ts.cfl_target == doctest::Approx(0.8));
    CHECK(ts.cfl_current == 0.0);
    CHECK(ts.step_count == 0);
}

TEST_CASE("TimeStepper: CFL helpers")
{
    TimeStepper ts;
    ts.dt = 0.1;
    // Advection: u*dt/dx = 1 * 0.1 / 0.5 = 0.2.
    CHECK(ts.compute_cfl_u(1.0, 0.5) == doctest::Approx(0.2));
    // Diffusion: nu*dt/dx^2 = 0.01 * 0.1 / 0.25 = 0.004.
    CHECK(ts.compute_cfl_nu(0.01, 0.5) == doctest::Approx(0.004));
}

TEST_CASE("TimeStepper: adapt moves dt toward the CFL target")
{
    TimeStepper ts;
    ts.dt = 1.0;
    ts.cfl_target = 0.5;
    ts.adapt(/*u_max*/ 2.0, /*nu*/ 0.0, /*dx*/ 1.0);
    // cfl = 2*1/1 = 2, ratio = 0.5/2 -> dt = 1 * 0.25 = 0.25.
    CHECK(ts.cfl_current == doctest::Approx(2.0));
    CHECK(ts.dt == doctest::Approx(0.25));
}

TEST_CASE("TimeStepper: adapt clamps dt to [dt_min, dt_max]")
{
    TimeStepper ts;
    ts.dt = 1.0;
    ts.dt_min = 0.1;
    ts.adapt(/*u_max*/ 1e12, /*nu*/ 0.0, /*dx*/ 1.0);
    // cfl = 1e12, dt would collapse to ~1e-12 but clamps to dt_min.
    CHECK(ts.cfl_current == doctest::Approx(1e12));
    CHECK(ts.dt == doctest::Approx(0.1));
}

TEST_CASE("TimeStepper: advance increments t and step count")
{
    TimeStepper ts;
    ts.dt = 0.25;
    ts.advance();
    CHECK(ts.t == doctest::Approx(0.25));
    CHECK(ts.step_count == 1);
    ts.advance();
    CHECK(ts.t == doctest::Approx(0.5));
    CHECK(ts.step_count == 2);
}

TEST_CASE("ConvergenceMonitor: residual below tolerance flips converged")
{
    ConvergenceMonitor cm;
    CHECK_FALSE(cm.converged);

    CHECK_FALSE(cm.check(1e-3)); // >= tolerance 1e-6
    CHECK_FALSE(cm.converged);
    CHECK(cm.residual == doctest::Approx(1e-3));

    CHECK(cm.check(1e-7)); // < tolerance 1e-6
    CHECK(cm.converged);
    CHECK(cm.residual == doctest::Approx(1e-7));
}