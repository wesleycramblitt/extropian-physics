#include <exd/engine/physics/rigid_body/dynamics.hpp>
#include <exd/engine/physics/rigid_body/rotational_state.hpp>
#include <exd/engine/physics/rigid_body/status.hpp>

#include <doctest/doctest.h>

#include <memory>

using namespace exd::engine::physics::rigid_body;

TEST_CASE("Rigid rotor Euler: omega grows by M/J*dt and angle advances by omega*dt")
{
    RigidRotorConfig config;
    config.inertia = 10.0;
    config.integration = RotationalIntegration::Euler;

    auto dynamics = make_rigid_rotor_dynamics(config);
    REQUIRE(dynamics);
    CHECK(dynamics->name() == "rigid_rotor");

    RotationalState state;
    state.omega = 5.0;
    state.angle_rad = 0.0;

    ModelStatus status;
    const RotationalState next = dynamics->advance(0.1, 20.0, state, status);

    REQUIRE(status.ok);
    // alpha = 20/10 = 2; omega_new = 5 + 2*0.1 = 5.2
    CHECK(next.omega == doctest::Approx(5.2).epsilon(1e-12));
    // angle += omega_old*dt = 5*0.1 = 0.5
    CHECK(next.angle_rad == doctest::Approx(0.5).epsilon(1e-12));
}

TEST_CASE("Rigid rotor Heun: angle uses the averaged omega and leads Euler")
{
    RigidRotorConfig euler_config;
    euler_config.inertia = 10.0;
    euler_config.integration = RotationalIntegration::Euler;

    RigidRotorConfig heun_config;
    heun_config.inertia = 10.0;
    heun_config.integration = RotationalIntegration::Heun;

    auto euler = make_rigid_rotor_dynamics(euler_config);
    auto heun = make_rigid_rotor_dynamics(heun_config);
    REQUIRE(euler);
    REQUIRE(heun);

    RotationalState state;
    state.omega = 5.0;

    ModelStatus euler_status;
    ModelStatus heun_status;
    const RotationalState euler_next = euler->advance(0.1, 20.0, state, euler_status);
    const RotationalState heun_next = heun->advance(0.1, 20.0, state, heun_status);

    REQUIRE(euler_status.ok);
    REQUIRE(heun_status.ok);

    // Both integrations agree on the updated omega.
    CHECK(heun_next.omega == doctest::Approx(euler_next.omega).epsilon(1e-12));
    CHECK(heun_next.omega == doctest::Approx(5.2).epsilon(1e-12));
    // Heun: theta += (5 + 5.2)/2 * 0.1 = 0.51
    CHECK(heun_next.angle_rad == doctest::Approx(0.51).epsilon(1e-12));
    CHECK(heun_next.angle_rad > euler_next.angle_rad);
}

TEST_CASE("Rigid rotor: default configuration integrates with Heun")
{
    RigidRotorConfig config;
    config.inertia = 10.0; // RotationalIntegration::Heun is the default

    auto dynamics = make_rigid_rotor_dynamics(config);
    REQUIRE(dynamics);

    RotationalState state;
    state.omega = 0.0;

    ModelStatus status;
    const RotationalState next = dynamics->advance(0.1, 20.0, state, status);

    REQUIRE(status.ok);
    CHECK(next.omega == doctest::Approx(0.2).epsilon(1e-12));
    // Heun: theta += (0 + 0.2)/2 * 0.1 = 0.01
    CHECK(next.angle_rad == doctest::Approx(0.01).epsilon(1e-12));
}

TEST_CASE("Rigid rotor: non-positive inertia is rejected and state is unchanged")
{
    RigidRotorConfig config;
    config.inertia = 0.0;

    auto dynamics = make_rigid_rotor_dynamics(config);
    REQUIRE(dynamics);

    RotationalState state;
    state.omega = 3.0;
    state.angle_rad = 1.5;

    ModelStatus status;
    const RotationalState next = dynamics->advance(0.1, 20.0, state, status);

    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(next.omega == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(next.angle_rad == doctest::Approx(1.5).epsilon(1e-12));
}

TEST_CASE("Rigid rotor: non-positive dt is rejected and state is unchanged")
{
    RigidRotorConfig config;
    config.inertia = 10.0;

    auto dynamics = make_rigid_rotor_dynamics(config);
    REQUIRE(dynamics);

    RotationalState state;
    state.omega = 3.0;
    state.angle_rad = 1.5;

    ModelStatus status;
    const RotationalState next = dynamics->advance(0.0, 20.0, state, status);

    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(next.omega == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(next.angle_rad == doctest::Approx(1.5).epsilon(1e-12));
}

TEST_CASE("Rigid rotor: negative net moment decelerates")
{
    RigidRotorConfig config;
    config.inertia = 10.0;
    config.integration = RotationalIntegration::Euler;

    auto dynamics = make_rigid_rotor_dynamics(config);
    REQUIRE(dynamics);

    RotationalState state;
    state.omega = 5.0;

    ModelStatus status;
    const RotationalState next = dynamics->advance(0.1, -20.0, state, status);

    REQUIRE(status.ok);
    // alpha = -2; omega_new = 5 - 0.2 = 4.8
    CHECK(next.omega == doctest::Approx(4.8).epsilon(1e-12));
    CHECK(next.omega < state.omega);
    // Euler still advances the angle with the pre-step omega.
    CHECK(next.angle_rad == doctest::Approx(0.5).epsilon(1e-12));
}