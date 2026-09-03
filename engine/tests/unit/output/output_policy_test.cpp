// output_policy_test.cpp
// Deterministic cadence logic with injected clock.

#include <exd/engine/output/output_policy.hpp>

#include <doctest/doctest.h>

using namespace exd::engine::output;

TEST_CASE("output policy: pure step trigger")
{
    OutputPolicy p;
    p.every_n_steps = 10;
    for (uint64_t s = 0; s < 100; ++s)
        CHECK(p.step_triggered(s) == (s % 10 == 0));

    CHECK_FALSE(OutputPolicy{}.enabled());
    CHECK(OutputPolicy{}.step_triggered(0) == false); // disabled: never
    CHECK(OutputPolicy{0, 0.1}.enabled());
}

TEST_CASE("output scheduler: step mode is deterministic and clock-independent")
{
    OutputScheduler sched(OutputPolicy{5, 0.0}, 123.456);
    CHECK(sched.should_emit(0));   // step 0 lands on the grid
    CHECK_FALSE(sched.should_emit(1));
    CHECK_FALSE(sched.should_emit(2));
    CHECK(sched.should_emit(5));
    CHECK_FALSE(sched.should_emit(6));
}

TEST_CASE("output scheduler: wall-clock mode throttles to a cadence")
{
    OutputScheduler sched(OutputPolicy{0, 0.5}, 0.0);
    CHECK(sched.should_emit(0));   // first call emits (gap >= interval)
    CHECK_FALSE(sched.should_emit(1)); // 0.0s elapsed
    sched.set_now(0.25);
    CHECK_FALSE(sched.should_emit(2));
    sched.set_now(0.5);
    CHECK(sched.should_emit(3));
    sched.set_now(0.75);
    CHECK_FALSE(sched.should_emit(4));
    sched.set_now(1.0);
    CHECK(sched.should_emit(5));
}

TEST_CASE("output scheduler: OR semantics with both triggers")
{
    OutputScheduler sched(OutputPolicy{10, 0.25}, 0.0);
    CHECK(sched.should_emit(0));   // step grid
    CHECK(sched.should_emit(10));  // step grid
    CHECK_FALSE(sched.should_emit(11));
    sched.set_now(0.26);           // wall clock fires between step gridlines
    CHECK(sched.should_emit(11));
    sched.set_now(0.27);
    CHECK_FALSE(sched.should_emit(12));
    sched.set_now(0.51);           // wall clock again
    CHECK(sched.should_emit(12));
}

TEST_CASE("output scheduler: disabled policy never emits")
{
    OutputScheduler sched(OutputPolicy{}, 0.0);
    CHECK_FALSE(sched.should_emit(0));
    sched.set_now(1000.0);
    CHECK_FALSE(sched.should_emit(99999));
}
