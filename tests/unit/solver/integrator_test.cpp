// integrator_test.cpp
// Unit tests for the shared ODE integrators (integrate_step).

#include <exd/physics/model_status.hpp>
#include <exd/physics/solver/integrators.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

using namespace exd::physics::solver;

namespace {
// ─────────────────────────────────────────────────────────────
// Shared fixtures
// ─────────────────────────────────────────────────────────────

// y' = -y
const auto decay_deriv = [](std::span<const double> s, std::span<double> d, double) {
    d[0] = -s[0];
};

// y' = -100 y (stiff decay)
const auto stiff_deriv = [](std::span<const double> s, std::span<double> d, double) {
    d[0] = -100.0 * s[0];
};

// y' = y (exponential growth)
const auto growth_deriv = [](std::span<const double> s, std::span<double> d, double) {
    d[0] = s[0];
};

// Harmonic oscillator y'' = -y as state [y, v] -> [v, -y]
const auto oscillator_deriv = [](std::span<const double> s, std::span<double> d, double) {
    d[0] = s[1];
    d[1] = -s[0];
};

/// Integrate y' = -y from y(0) = 1 for `steps` steps of `dt`.
double run_decay(IntegrationMethod method, double dt, std::size_t steps)
{
    IntegratorConfig config;
    config.method = method;
    std::vector<double> state{1.0};
    exd::physics::ModelStatus status;
    for (std::size_t i = 0; i < steps; ++i)
    {
        REQUIRE(integrate_step(config, 0.0, dt, state, decay_deriv, status));
    }
    return state[0];
}

} // namespace

TEST_CASE("ForwardEuler is first-order accurate on y' = -y")
{
    const double y_coarse = run_decay(IntegrationMethod::ForwardEuler, 0.1, 10);
    const double y_fine = run_decay(IntegrationMethod::ForwardEuler, 0.05, 20);
    const double exact = std::exp(-1.0);
    const double err_coarse = std::fabs(y_coarse - exact);
    const double err_fine = std::fabs(y_fine - exact);
    CHECK(err_coarse > 1e-6);
    CHECK(err_fine > 1e-6);
    CHECK(err_coarse / err_fine == doctest::Approx(2.0).epsilon(0.2));
}

TEST_CASE("Heun is second-order accurate on y' = -y")
{
    const double y_coarse = run_decay(IntegrationMethod::Heun, 0.1, 10);
    const double y_fine = run_decay(IntegrationMethod::Heun, 0.05, 20);
    const double exact = std::exp(-1.0);
    const double err_coarse = std::fabs(y_coarse - exact);
    const double err_fine = std::fabs(y_fine - exact);
    CHECK(err_coarse > 1e-4);
    CHECK(err_fine > 1e-5);
    CHECK(err_coarse / err_fine == doctest::Approx(4.0).epsilon(0.2));
}

TEST_CASE("RK4 is fourth-order accurate on y' = -y")
{
    const double y_coarse = run_decay(IntegrationMethod::RK4, 0.1, 10);
    const double y_fine = run_decay(IntegrationMethod::RK4, 0.05, 20);
    const double exact = std::exp(-1.0);
    const double err_coarse = std::fabs(y_coarse - exact);
    const double err_fine = std::fabs(y_fine - exact);
    CHECK(err_coarse > 1e-7);
    CHECK(err_fine > 1e-9);
    CHECK(err_coarse / err_fine == doctest::Approx(16.0).epsilon(0.25));
}

TEST_CASE("RK4 conserves energy on the harmonic oscillator")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::RK4;

    std::vector<double> state{1.0, 0.0};
    exd::physics::ModelStatus status;
    const double dt = 0.01;
    for (int i = 0; i < 200; ++i)
    {
        REQUIRE(integrate_step(config, 0.0, dt, state, oscillator_deriv, status));
    }
    const double energy = state[0] * state[0] + state[1] * state[1];
    CHECK(energy == doctest::Approx(1.0).epsilon(1e-5));
}

TEST_CASE("SymplecticEuler: bounded energy drift on harmonic oscillator")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::SymplecticEuler;
    config.position_count = 1;

    std::vector<double> state{1.0, 0.0};
    exd::physics::ModelStatus status;
    const double dt = 0.05;
    for (int i = 0; i < 1000; ++i)
    {
        REQUIRE(integrate_step(config, 0.0, dt, state, oscillator_deriv, status));
    }
    const double energy = state[0] * state[0] + state[1] * state[1];
    CHECK(energy == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("Verlet: near-zero energy drift on harmonic oscillator")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::Verlet;
    config.position_count = 1;

    std::vector<double> state{1.0, 0.0};
    exd::physics::ModelStatus status;
    const double dt = 0.05;
    for (int i = 0; i < 1000; ++i)
    {
        REQUIRE(integrate_step(config, 0.0, dt, state, oscillator_deriv, status));
    }
    const double energy = state[0] * state[0] + state[1] * state[1];
    CHECK(energy == doctest::Approx(1.0).epsilon(0.005));
}

TEST_CASE("BackwardEuler: A-stable on stiff decay")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::BackwardEuler;
    // |lambda*dt| = 10 > 1, so under-relaxation is required for the
    // fixed-point iteration to contract.
    config.relaxation = 0.1;

    std::vector<double> state{1.0};
    exd::physics::ModelStatus status;
    REQUIRE(integrate_step(config, 0.0, 0.1, state, stiff_deriv, status));
    REQUIRE(status.ok);
    // BE: y1 = y0 / (1 + lambda*dt) = 1 / 11.
    CHECK(state[0] == doctest::Approx(1.0 / 11.0).epsilon(0.001));
}

TEST_CASE("CrankNicolson: stable on stiff decay")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::CrankNicolson;
    config.relaxation = 0.15; // contract fixed-point iteration for stiff case

    std::vector<double> state{1.0};
    exd::physics::ModelStatus status;
    REQUIRE(integrate_step(config, 0.0, 0.1, state, stiff_deriv, status));
    REQUIRE(status.ok);
    // CN: y1 = y0 * (1 - lambda*dt/2) / (1 + lambda*dt/2) = -2/3. Oscillatory
    // but bounded (|y1| <= 1).
    CHECK(std::fabs(state[0]) <= 1.0);
    CHECK(std::fabs(state[0]) == doctest::Approx(2.0 / 3.0).epsilon(0.01));
}

// y' = 1 (linear ramp): DP45 integrates it exactly, local error ~ 0, so the
// controller happily accepts larger steps.
const auto ramp_deriv = [](std::span<const double>, std::span<double> d, double) {
    d[0] = 1.0;
};

TEST_CASE("AdaptiveRK45: tracks exponential growth accurately")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::AdaptiveRK45;
    config.rel_tol = 1e-9;
    config.abs_tol = 1e-12;

    std::vector<double> state{1.0};
    exd::physics::ModelStatus status;

    double t = 0.0;
    double dt = 0.2;
    while (t < 2.0)
    {
        const double dt_try = std::min(dt, 2.0 - t);
        double dt_used = 0.0;
        REQUIRE(integrate_step(config, t, dt_try, state, growth_deriv, status, &dt_used));
        REQUIRE(status.ok);
        t += dt_used;
        CHECK(dt_used >= config.dt_min);
        CHECK(dt_used <= config.dt_max);
        CHECK(dt_used <= dt_try); // integrate_step only shrinks below the request
        dt = dt_used * 1.5;       // caller-side growth policy
    }

    // Accuracy despite aggressive requests: relative error control keeps the
    // global error bounded (per-step local error ~ tol * |y|).
    CHECK(t == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(state[0] == doctest::Approx(std::exp(2.0)).epsilon(1e-7));
}

TEST_CASE("AdaptiveRK45: grows the step on a low-error solution")
{
    IntegratorConfig config;
    config.method = IntegrationMethod::AdaptiveRK45;
    config.rel_tol = 1e-9;
    config.abs_tol = 1e-12;

    // y' = 1 -> DP45 error estimate ~ 0 -> every requested step is accepted
    // as-is, so the caller-side growth policy takes effect and dt keeps
    // growing (contrast: relative error control keeps dt constant for pure
    // exponentials because the local error scales with |y|).
    std::vector<double> state{0.0};
    exd::physics::ModelStatus status;

    double t = 0.0;
    double dt = 0.01;
    std::vector<double> dt_used_seq;
    while (t < 2.0)
    {
        const double dt_try = std::min(dt, 2.0 - t);
        double dt_used = 0.0;
        REQUIRE(integrate_step(config, t, dt_try, state, ramp_deriv, status, &dt_used));
        REQUIRE(status.ok);
        t += dt_used;
        dt_used_seq.push_back(dt_used);
        dt = dt_used * 1.5;
    }

    CHECK(t == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(state[0] == doctest::Approx(2.0).epsilon(1e-9)); // y = t integrated exactly

    REQUIRE(dt_used_seq.size() > 4);
    const double first_used = dt_used_seq.front();
    const double max_used = *std::max_element(dt_used_seq.begin(), dt_used_seq.end());
    CHECK(max_used > first_used);
    CHECK(max_used > 0.5); // visibly larger than the 0.01 start
}

TEST_CASE("Integrate step rejects invalid inputs without touching state")
{
    IntegratorConfig config; // default RK4

    SUBCASE("empty state")
    {
        std::vector<double> state;
        exd::physics::ModelStatus status;
        CHECK_FALSE(integrate_step(config, 0.0, 0.1, state, decay_deriv, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("dt equal to zero")
    {
        std::vector<double> state{1.0};
        exd::physics::ModelStatus status;
        CHECK_FALSE(integrate_step(config, 0.0, 0.0, state, decay_deriv, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(state[0] == doctest::Approx(1.0));
    }

    SUBCASE("non-callable derivative")
    {
        std::vector<double> state{1.0};
        exd::physics::ModelStatus status;
        DerivativeFn empty_deriv;
        CHECK_FALSE(integrate_step(config, 0.0, 0.1, state, empty_deriv, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(state[0] == doctest::Approx(1.0));
    }

    SUBCASE("position_count exceeds state size")
    {
        std::vector<double> state{1.0};
        exd::physics::ModelStatus status;
        config.position_count = 2;
        CHECK_FALSE(integrate_step(config, 0.0, 0.1, state, decay_deriv, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(state[0] == doctest::Approx(1.0));
    }

    SUBCASE("adaptive config rejects dt_min not positive")
    {
        std::vector<double> state{1.0};
        exd::physics::ModelStatus status;
        config.method = IntegrationMethod::AdaptiveRK45;
        config.dt_min = -1.0;
        CHECK_FALSE(integrate_step(config, 0.0, 0.1, state, growth_deriv, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(state[0] == doctest::Approx(1.0));
    }
}