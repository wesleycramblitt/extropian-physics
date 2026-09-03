#include <exd/engine/physics/rigid_body/quaternion.hpp>
#include <exd/engine/physics/rigid_body/rigid_body.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <memory>

using namespace exd::engine::physics::rigid_body;

namespace
{

std::unique_ptr<IRigidBodyDynamics> make_advance(RigidBodyIntegration method,
                                                 double mass,
                                                 double ixx,
                                                 double iyy,
                                                 double izz)
{
    RigidBodyConfig config;
    config.mass = mass;
    config.inertia_principal = {ixx, iyy, izz};
    return make_rigid_body_dynamics(method, config);
}

RigidBodyState rest_state()
{
    RigidBodyState state;
    state.position = {0.0, 0.0, 0.0};
    state.orientation = {1.0, 0.0, 0.0, 0.0};
    state.linear_velocity = {0.0, 0.0, 0.0};
    state.angular_velocity = {0.0, 0.0, 0.0};
    return state;
}

double vector_magnitude(const std::array<double, 3>& v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

std::array<double, 3> add(const std::array<double, 3>& a,
                          const std::array<double, 3>& b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

// Kinetic + linear-spring potential of the anchor-at-origin oscillator.
double oscillator_energy(const RigidBodyState& state, double mass, double stiffness)
{
    const double vx = state.linear_velocity[0];
    const double vy = state.linear_velocity[1];
    const double vz = state.linear_velocity[2];
    const double px = state.position[0];
    const double py = state.position[1];
    const double pz = state.position[2];
    return 0.5 * mass * (vx * vx + vy * vy + vz * vz)
         + 0.5 * stiffness * (px * px + py * py + pz * pz);
}

} // anonymous namespace

TEST_CASE("Rigid body free fall under constant gravity")
{
    constexpr double kDt = 0.01;
    constexpr int kSteps = 100; // t = 1 s
    const std::array<double, 3> kGravity = {0.0, 0.0, -9.81};

    SUBCASE("RK4 matches the continuum 0.5*g*t^2")
    {
        auto dynamics = make_advance(RigidBodyIntegration::RK4, 1.0, 1.0, 1.0, 1.0);
        REQUIRE(dynamics);
        CHECK(dynamics->name() == "rigid_body");

        RigidBodyState state = rest_state();
        for (int i = 0; i < kSteps; ++i)
        {
            RigidBodyForces loads;
            loads.force = force_gravity(1.0, kGravity);
            ModelStatus status;
            state = dynamics->advance(kDt, loads, state, status);
            REQUIRE(status.ok);
        }

        CHECK(state.position[2] == doctest::Approx(-4.905).epsilon(1e-3));
        CHECK(state.linear_velocity[2] == doctest::Approx(-9.81).epsilon(1e-3));
        CHECK(state.position[0] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(state.position[1] == doctest::Approx(0.0).epsilon(1e-12));
    }

    SUBCASE("SymplecticEuler tracks the free fall with bounded drift")
    {
        auto dynamics = make_advance(RigidBodyIntegration::SymplecticEuler, 1.0, 1.0, 1.0, 1.0);
        REQUIRE(dynamics);

        RigidBodyState state = rest_state();
        for (int i = 0; i < kSteps; ++i)
        {
            RigidBodyForces loads;
            loads.force = force_gravity(1.0, kGravity);
            ModelStatus status;
            state = dynamics->advance(kDt, loads, state, status);
            REQUIRE(status.ok);
        }

        // Semi-implicit positional error is +0.5·a·dt·T, about +0.049 here
        // (≈1% of the continuum value), so ε = 1e-2 covers it. Velocity is
        // integrated exactly and needs no slack beyond ε = 1e-2.
        CHECK(state.position[2] == doctest::Approx(-4.905).epsilon(1e-2));
        CHECK(state.linear_velocity[2] == doctest::Approx(-9.81).epsilon(1e-2));
        CHECK(state.position[0] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(state.position[1] == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("Rigid body torque-free spin about a principal axis")
{
    // No net torque and ω initially aligned with a principal axis means
    // ω × I·ω = 0: |ω| is conserved and the orientation precesses uniformly.
    constexpr double kSpinRate = 5.0;  // rad/s
    constexpr double kTotalTime = 0.1; // small so no wrap-around ambiguity
    constexpr int kSteps = 100;
    constexpr double kDt = kTotalTime / kSteps;

    const std::array<double, 3> kBodyAxis = {1.0, 0.0, 0.0};

    SUBCASE("RK4 preserves |omega| and rotates by omega*t within 1e-6")
    {
        auto dynamics = make_advance(RigidBodyIntegration::RK4, 1.0, 1.0, 1.0, 2.0);
        REQUIRE(dynamics);

        RigidBodyState state = rest_state();
        state.angular_velocity = {0.0, 0.0, kSpinRate};

        for (int i = 0; i < kSteps; ++i)
        {
            RigidBodyForces loads; // zero force and torque
            ModelStatus status;
            state = dynamics->advance(kDt, loads, state, status);
            REQUIRE(status.ok);
        }

        CHECK(vector_magnitude(state.angular_velocity) == doctest::Approx(kSpinRate).epsilon(1e-9));

        const std::array<double, 3> world_axis = quat_rotate(state.orientation, kBodyAxis);
        const double angle = std::atan2(world_axis[1], world_axis[0]);
        CHECK(angle == doctest::Approx(kSpinRate * kTotalTime).epsilon(1e-6));
    }

    SUBCASE("SymplecticEuler preserves |omega| and rotates by omega*t")
    {
        auto dynamics = make_advance(RigidBodyIntegration::SymplecticEuler, 1.0, 1.0, 1.0, 2.0);
        REQUIRE(dynamics);

        RigidBodyState state = rest_state();
        state.angular_velocity = {0.0, 0.0, kSpinRate};

        for (int i = 0; i < kSteps; ++i)
        {
            RigidBodyForces loads; // zero force and torque
            ModelStatus status;
            state = dynamics->advance(kDt, loads, state, status);
            REQUIRE(status.ok);
        }

        CHECK(vector_magnitude(state.angular_velocity) == doctest::Approx(kSpinRate).epsilon(1e-9));

        // Normalized first-order orientation integration leaves an angle error
        // of ~1e-6 at this step size, so a wider ε than RK4 is required.
        const std::array<double, 3> world_axis = quat_rotate(state.orientation, kBodyAxis);
        const double angle = std::atan2(world_axis[1], world_axis[0]);
        CHECK(angle == doctest::Approx(kSpinRate * kTotalTime).epsilon(1e-3));
    }
}

TEST_CASE("Rigid body spring oscillator conserves mechanical energy")
{
    // Anchor at the origin, m = 1 kg, k = 4 N/m, starting at p0 = (1,0,0)
    // at rest. No applied gravity: the preserved quantity asserted here is
    // kinetic + spring potential only (with constant gravity the equilibrium
    // shifts and E would include an extra g·p term).
    constexpr double kMass = 1.0;
    constexpr double kStiffness = 4.0;
    constexpr double kDt = 0.01;
    constexpr int kSteps = 500; // t = 5 s ≈ 1.6 periods at ω = 2 rad/s

    RigidBodyState state = rest_state();
    state.position = {1.0, 0.0, 0.0};
    const double energy0 = oscillator_energy(state, kMass, kStiffness);

    // Loads recomputed at the start of every step (advance() treats them as
    // constant for the duration of the step).
    const auto run = [&](RigidBodyIntegration method)
    {
        auto dynamics = make_advance(method, kMass, 1.0, 1.0, 1.0);
        double max_deviation = 0.0;
        for (int i = 0; i < kSteps; ++i)
        {
            RigidBodyForces loads;
            loads.force = add(force_gravity(kMass, {0.0, 0.0, 0.0}),
                              force_linear_spring(state.position, {0.0, 0.0, 0.0}, kStiffness));
            ModelStatus status;
            state = dynamics->advance(kDt, loads, state, status);
            REQUIRE(status.ok);

            const double deviation =
                std::fabs(oscillator_energy(state, kMass, kStiffness) - energy0);
            max_deviation = std::max(max_deviation, deviation);
        }
        return max_deviation;
    };

    SUBCASE("RK4 drift stays below 20% of E0")
    {
        // With loads frozen per step, RK4 on a spring reduces to the exact
        // constant-acceleration map whose phase-space expansion is
        // (1 + 0.5·(ω·dt)²)^N ≈ 1.10 at these settings; 20% covers it with
        // margin. A state-dependent-force RK4 (not expressible through the
        // advance() contract) would achieve O((ω·dt)⁵) here instead.
        const double max_deviation = run(RigidBodyIntegration::RK4);
        CHECK(max_deviation < 0.2 * energy0);
    }

    SUBCASE("SymplecticEuler drift stays below 5% of E0")
    {
        const double max_deviation = run(RigidBodyIntegration::SymplecticEuler);
        CHECK(max_deviation < 0.05 * energy0);
    }
}

TEST_CASE("Rigid body: non-positive mass is rejected and state unchanged")
{
    auto dynamics = make_advance(RigidBodyIntegration::RK4, 0.0, 1.0, 1.0, 1.0);
    REQUIRE(dynamics);

    RigidBodyState state = rest_state();
    state.position = {1.0, 2.0, 3.0};
    state.orientation = {0.5, 0.5, 0.5, 0.5};
    state.linear_velocity = {0.5, -0.5, 0.25};
    state.angular_velocity = {1.0, 2.0, 3.0};

    RigidBodyForces loads;
    loads.force = {10.0, 0.0, 0.0};

    ModelStatus status;
    const RigidBodyState next = dynamics->advance(0.01, loads, state, status);

    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(next.position == state.position);
    CHECK(next.orientation == state.orientation);
    CHECK(next.linear_velocity == state.linear_velocity);
    CHECK(next.angular_velocity == state.angular_velocity);
}

TEST_CASE("Rigid body: non-positive principal inertia is rejected and state unchanged")
{
    auto dynamics = make_advance(RigidBodyIntegration::SymplecticEuler, 1.0, 0.0, 1.0, 1.0);
    REQUIRE(dynamics);

    RigidBodyState state = rest_state();
    state.linear_velocity = {1.0, 0.0, 0.0};

    RigidBodyForces loads;

    ModelStatus status;
    const RigidBodyState next = dynamics->advance(0.01, loads, state, status);

    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(next.position == state.position);
    CHECK(next.orientation == state.orientation);
    CHECK(next.linear_velocity == state.linear_velocity);
    CHECK(next.angular_velocity == state.angular_velocity);
}

TEST_CASE("Rigid body: non-positive dt is rejected and state unchanged")
{
    auto dynamics = make_advance(RigidBodyIntegration::RK4, 1.0, 1.0, 1.0, 1.0);
    REQUIRE(dynamics);

    SUBCASE("zero dt")
    {
        RigidBodyState state = rest_state();
        state.linear_velocity = {1.0, 0.0, 0.0};

        RigidBodyForces loads;
        ModelStatus status;
        const RigidBodyState next = dynamics->advance(0.0, loads, state, status);

        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(next.linear_velocity == state.linear_velocity);
    }

    SUBCASE("negative dt")
    {
        RigidBodyState state = rest_state();
        state.linear_velocity = {1.0, 0.0, 0.0};

        RigidBodyForces loads;
        ModelStatus status;
        const RigidBodyState next = dynamics->advance(-0.1, loads, state, status);

        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(next.linear_velocity == state.linear_velocity);
    }
}

TEST_CASE("Rigid body load helpers: gravity and linear spring")
{
    const std::array<double, 3> gravity_force = force_gravity(2.0, {0.0, 0.0, -9.81});
    CHECK(gravity_force[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(gravity_force[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(gravity_force[2] == doctest::Approx(-19.62).epsilon(1e-12));

    const std::array<double, 3> spring_force =
        force_linear_spring({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 3.0);
    CHECK(spring_force[0] == doctest::Approx(-3.0).epsilon(1e-12));
    CHECK(spring_force[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(spring_force[2] == doctest::Approx(0.0).epsilon(1e-12));
}