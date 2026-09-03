// circuit_test.cpp
// Unit tests for the lumped DC machine models (Phase D part 3): the dynamic
// DcMotorModel (exact exponential RL armature integration) and the
// quasi-steady DcMotorMomentModel exposed as a exd::engine::physics::rigid_body::IMomentModel load,
// including an end-to-end coupled test through RotatingAssembly.

#include <exd/engine/physics/electromagnetics/circuit.hpp>

#include <exd/engine/physics/rigid_body/dynamics.hpp>
#include <exd/engine/physics/rigid_body/rotating_assembly.hpp>
#include <exd/engine/physics/rigid_body/rotational_state.hpp>
#include <exd/engine/physics/rigid_body/status.hpp>
#include <exd/engine/numerics/integrators.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <span>
#include <vector>

using namespace exd::engine::physics::electromagnetics;
using namespace exd::engine::physics::rigid_body;
using namespace exd::engine::numerics;

namespace
{
// ─────────────────────────────────────────────────────────────
// Shared fixtures
// ─────────────────────────────────────────────────────────────

// Constant "aero" drive used by the coupled RotatingAssembly test (N·m).
constexpr double kCoupledAeroTorque = 5.0;

} // anonymous namespace

TEST_CASE("DcMotorModel: RL armature current rises exponentially toward V/R")
{
    DcMotorConfig config;
    config.kt = 1.0;
    config.ke = 0.0; // no back-emf: pure RL circuit V = R·i + L·di/dt
    config.R = 10.0;
    config.L = 2.0; // τ = L/R = 0.2 s
    config.v_supply = 100.0;

    DcMotorModel motor(config);
    REQUIRE(motor.valid());
    CHECK(motor.current() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(motor.torque() == doctest::Approx(0.0).epsilon(1e-12));

    exd::engine::core::ModelStatus status;

    // i_inf = 100/10 = 10 A. One step of 0.1 s:
    //   i = 10·(1 − exp(−10·0.1/2)) = 10·(1 − e^{−0.5}) ≈ 3.93469 A.
    const double i1 = motor.step(0.1, 0.0, status);
    REQUIRE(status.ok);
    CHECK(i1 == doctest::Approx(10.0 * (1.0 - std::exp(-0.5))).epsilon(1e-9));
    CHECK(motor.current() == doctest::Approx(3.9346934).epsilon(1e-6));
    CHECK(motor.torque() == doctest::Approx(3.9346934).epsilon(1e-6));

    // Repeated steps follow the same exponential and approach steady state
    // (the decay is implicit in the rise: after 50 steps of 0.1 s = 5 s).
    for (int step = 0; step < 50; ++step)
    {
        motor.step(0.1, 0.0, status);
    }
    REQUIRE(status.ok);
    CHECK(motor.current() == doctest::Approx(10.0).epsilon(1e-3));
    CHECK(motor.torque() == doctest::Approx(10.0).epsilon(1e-3));
}

TEST_CASE("DcMotorModel: back-emf pulls the current to (V − ke·omega)/R at constant speed")
{
    DcMotorConfig config;
    config.kt = 1.0;
    config.ke = 1.0;
    config.R = 1.0;
    config.L = 0.1; // τ = L/R = 0.1 s
    config.v_supply = 10.0;

    DcMotorModel motor(config);
    REQUIRE(motor.valid());

    exd::engine::core::ModelStatus status;
    const double omega = 5.0; // target: (10 − 1·5)/1 = 5 A

    // 200 steps of 1 ms = 0.2 s = 2τ: the exponential is only ~86% there
    // (i = 5·(1 − e^{−2}) ≈ 4.32332 A), so check the analytic intermediate.
    for (int step = 0; step < 200; ++step)
    {
        motor.step(0.001, omega, status);
    }
    REQUIRE(status.ok);
    CHECK(motor.current() == doctest::Approx(5.0 * (1.0 - std::exp(-2.0))).epsilon(1e-6));

    // 2000 total steps = 2 s = 20τ: settled on the back-emf equilibrium.
    for (int step = 0; step < 1800; ++step)
    {
        motor.step(0.001, omega, status);
    }
    REQUIRE(status.ok);
    CHECK(motor.current() == doctest::Approx(5.0).epsilon(1e-3));
    CHECK(motor.torque() == doctest::Approx(5.0).epsilon(1e-3));
    CHECK(motor.torque() > 0.0); // motoring: current and torque are positive
}

TEST_CASE("DcMotorModel: L=0 evaluates the quasi-steady current in a single step")
{
    DcMotorConfig config;
    config.kt = 2.0;
    config.ke = 0.5;
    config.R = 4.0;
    config.L = 0.0; // quasi-steady path
    config.v_supply = 12.0;

    exd::engine::core::ModelStatus status;

    SUBCASE("motoring at omega = 10 rad/s")
    {
        DcMotorModel motor(config);
        const double torque = motor.step(0.1, 10.0, status);
        REQUIRE(status.ok);
        // i = (12 − 0.5·10)/4 = 1.75 A; T = kt·i = 2·1.75 = 3.5 N·m.
        CHECK(motor.current() == doctest::Approx(1.75).epsilon(1e-12));
        CHECK(torque == doctest::Approx(3.5).epsilon(1e-12));
        CHECK(torque > 0.0);
    }

    SUBCASE("motoring at stall: positive torque assists the drive")
    {
        DcMotorModel motor(config);
        const double torque = motor.step(0.1, 0.0, status);
        REQUIRE(status.ok);
        // i = 12/4 = 3.0 A; T = 6.0 N·m. At omega = 0 with v_supply > 0 the
        // torque is positive (motoring/assisting).
        CHECK(motor.current() == doctest::Approx(3.0).epsilon(1e-12));
        CHECK(torque == doctest::Approx(6.0).epsilon(1e-12));
        CHECK(torque > 0.0);
    }
}

TEST_CASE("DcMotorModel: generating at high speed brakes with a negative torque")
{
    DcMotorConfig config;
    config.kt = 2.0;
    config.ke = 0.5;
    config.R = 4.0;
    config.L = 0.0; // quasi-steady path
    config.v_supply = 0.0; // no supply: driven as a generator

    DcMotorModel motor(config);
    exd::engine::core::ModelStatus status;

    const double torque = motor.step(0.1, 20.0, status);
    REQUIRE(status.ok);
    // i = (0 − 0.5·20)/4 = −2.5 A; T = 2·(−2.5) = −5 N·m (opposing/braking).
    CHECK(motor.current() == doctest::Approx(-2.5).epsilon(1e-12));
    CHECK(torque == doctest::Approx(-5.0).epsilon(1e-12));
    CHECK(torque < 0.0);
}

TEST_CASE("DcMotorModel: reset zeroes the internal current and torque")
{
    DcMotorConfig config;
    config.kt = 1.0;
    config.ke = 0.0;
    config.R = 1.0;
    config.L = 0.0;
    config.v_supply = 10.0;

    DcMotorModel motor(config);
    exd::engine::core::ModelStatus status;

    motor.step(0.1, 0.0, status);
    REQUIRE(status.ok);
    CHECK(motor.current() == doctest::Approx(10.0).epsilon(1e-12));
    CHECK(motor.torque() == doctest::Approx(10.0).epsilon(1e-12));

    motor.reset();
    CHECK(motor.current() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(motor.torque() == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("DcMotorModel: invalid configuration and step size are rejected through status")
{
    exd::engine::core::ModelStatus status;

    SUBCASE("non-positive R is invalid")
    {
        DcMotorConfig config;
        config.kt = 1.0;
        config.ke = 0.0;
        config.R = 0.0;
        config.L = 0.0;
        DcMotorModel motor(config);
        CHECK_FALSE(motor.valid());
        const double torque = motor.step(0.1, 0.0, status);
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(torque == doctest::Approx(0.0).epsilon(1e-12));
    }

    SUBCASE("non-positive kt is invalid")
    {
        DcMotorConfig config;
        config.kt = 0.0;
        config.ke = 1.0;
        config.R = 1.0;
        config.L = 0.1;
        config.v_supply = 10.0;
        DcMotorModel motor(config);
        CHECK_FALSE(motor.valid());
        const double torque = motor.step(0.1, 0.0, status);
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(torque == doctest::Approx(0.0).epsilon(1e-12));
    }

    SUBCASE("negative inductance or back-emf constant is invalid")
    {
        DcMotorConfig neg_l;
        neg_l.kt = 1.0;
        neg_l.ke = 0.0;
        neg_l.R = 1.0;
        neg_l.L = -1.0;
        CHECK_FALSE(DcMotorModel(neg_l).valid());

        DcMotorConfig neg_ke;
        neg_ke.kt = 1.0;
        neg_ke.ke = -1.0;
        neg_ke.R = 1.0;
        neg_ke.L = 0.0;
        CHECK_FALSE(DcMotorModel(neg_ke).valid());
    }

    SUBCASE("non-positive dt is rejected and the stored torque is returned")
    {
        DcMotorConfig config;
        config.kt = 1.0;
        config.ke = 0.0;
        config.R = 1.0;
        config.L = 0.0;
        config.v_supply = 10.0;
        DcMotorModel motor(config);
        const double torque = motor.step(0.0, 0.0, status);
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
        CHECK(torque == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("make_dc_motor_moment: quasi-steady load law T(omega) = kt*(ke*omega − V)/R")
{
    DcMotorConfig config;
    config.kt = 1.0;
    config.ke = 1.0;
    config.R = 2.0;
    config.v_supply = 20.0;

    auto model = make_dc_motor_moment(config);
    REQUIRE(model);
    CHECK(model->name() == "dc_motor");

    RotationalState state;
    ModelStatus status;

    // IMomentModel convention: positive opposes rotation.
    state.omega = 0.0; // motoring stall: (0 − 20)/2 = −10 (assists)
    CHECK(model->moment(state, status) == doctest::Approx(-10.0).epsilon(1e-12));
    state.omega = 10.0; // (10 − 20)/2 = −5 (still assisting)
    CHECK(model->moment(state, status) == doctest::Approx(-5.0).epsilon(1e-12));
    state.omega = 20.0; // synchronous speed: no net armature current
    CHECK(model->moment(state, status) == doctest::Approx(0.0).epsilon(1e-12));
    state.omega = 30.0; // generating: (30 − 20)/2 = +5 opposing (regenerative braking)
    CHECK(model->moment(state, status) == doctest::Approx(5.0).epsilon(1e-12));

    REQUIRE(status.ok);
    CHECK(status.error.empty());
}

TEST_CASE("make_dc_motor_moment: non-positive R produces a null model")
{
    DcMotorConfig config;
    config.kt = 1.0;
    config.ke = 0.0;
    config.R = 0.0;
    CHECK_FALSE(make_dc_motor_moment(config));

    config.R = -2.0;
    CHECK_FALSE(make_dc_motor_moment(config));
}

TEST_CASE("Shared integrator reproduces the RL armature ODE")
{
    const double ke = 0.0;
    const double R = 10.0;
    const double L = 2.0;
    const double v_supply = 100.0;
    const double omega = 0.0;

    IntegratorConfig config;
    config.method = IntegrationMethod::RK4;

    // L·di/dt + R·i = v_supply − ke·ω  ⇒  di/dt = (v_supply − ke·ω − R·i)/L.
    auto armature_deriv = [=](std::span<const double> state,
                              std::span<double> dstate, double) {
        dstate[0] = (v_supply - ke * omega - R * state[0]) / L;
    };

    std::vector<double> current{0.0};
    exd::engine::core::ModelStatus status;
    REQUIRE(integrate_step(config, 0.0, 0.1, current, armature_deriv, status));

    // One RK4 step of the same ODE lands within 1e-3 of the exact
    // exponential solution used by DcMotorModel (3.93469 A).
    CHECK(current[0] == doctest::Approx(3.9347).epsilon(1e-3));
}

TEST_CASE("DcMotorMomentModel inside a RotatingAssembly: motor/generator coupling")
{
    // Constant "aero" drive; the provider captures nothing and the assembly
    // threads the rotational state from step to step.
    RotatingAssembly::ForceProvider provider =
        [](const RotationalState&, ModelStatus&) {
            AeroResult aero;
            aero.moments.valid = true;
            aero.moments.torque = kCoupledAeroTorque;
            aero.moments.axial_force = 0.0;
            return aero;
        };

    RigidRotorConfig dyn_config;
    dyn_config.inertia = 10.0;
    dyn_config.integration = RotationalIntegration::Heun;
    auto dynamics = make_rigid_rotor_dynamics(dyn_config);
    REQUIRE(dynamics);

    const double dt = 0.01;

    SUBCASE("motoring load assists the drive through a negative external moment")
    {
        // Motor (V = 20): T_load = (ω − 20)/2 — negative below synchronous
        // speed, so it ASSISTS the aero drive (net = aero − external grows).
        DcMotorConfig motor_config;
        motor_config.kt = 1.0;
        motor_config.ke = 1.0;
        motor_config.R = 2.0;
        motor_config.v_supply = 20.0;
        auto external = make_dc_motor_moment(motor_config);
        REQUIRE(external);

        RotatingAssembly assembly(provider, std::move(external), std::move(dynamics));

        RotationalState state; // stall: omega = 0
        const AssemblyStepResult result = assembly.step(dt, state);
        REQUIRE(result.ok);
        REQUIRE(result.status.ok);
        // T_load(0) = (0 − 20)/2 = −10 (assisting), aero 5: net = 5 − (−10) = 15
        CHECK(result.external_moment == doctest::Approx(-10.0).epsilon(1e-12));
        CHECK(result.net_moment == doctest::Approx(15.0).epsilon(1e-12));
        CHECK(result.state.omega == doctest::Approx(0.015).epsilon(1e-12));
    }

    SUBCASE("generator load (V = 0) spin-up: turbine drives the generator to equilibrium")
    {
        // Generator (V = 0): T_load = kt·ke·ω/R = ω/2 grows with omega, so the
        // equilibrium aero = T_load(ω*) ⇒ 5 = ω*/2 ⇒ ω* = 10 is ATTRACTING
        // (d(net)/dω = −0.5 < 0). The rotor spins up from rest and settles —
        // the physical wind-turbine-generator story, end to end.
        DcMotorConfig gen_config;
        gen_config.kt = 1.0;
        gen_config.ke = 1.0;
        gen_config.R = 2.0;
        gen_config.v_supply = 0.0;
        auto external = make_dc_motor_moment(gen_config);
        REQUIRE(external);

        RotatingAssembly assembly(provider, std::move(external), std::move(dynamics));

        RotationalState state; // from rest
        for (int step = 0; step < 15000; ++step)
        {
            const AssemblyStepResult result = assembly.step(dt, state);
            REQUIRE(result.ok);
            REQUIRE(result.status.ok);
            state = result.state;
        }

        CHECK(state.omega == doctest::Approx(10.0).epsilon(1e-3));
    }
}