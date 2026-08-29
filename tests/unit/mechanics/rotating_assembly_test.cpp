#include <exd/physics/mechanics/rotating_assembly.hpp>

#include <doctest/doctest.h>

#include <memory>

using namespace exd::physics::mechanics;

namespace
{

AeroResult make_load_result()
{
    AeroResult aero;
    aero.moments.valid = true;
    aero.moments.torque = 100.0;
    aero.moments.axial_force = 50.0;

    ElementForce3D blade;
    blade.r = 2.0;
    blade.ref = {2.0, 0.0, 0.0};
    blade.force = {0.0, 50.0, 50.0};
    aero.per_element.push_back(blade);
    return aero;
}

std::unique_ptr<IRotationalDynamics> make_dynamics()
{
    RigidRotorConfig config;
    config.inertia = 10.0;
    config.integration = RotationalIntegration::Heun;
    return make_rigid_rotor_dynamics(config);
}

} // anonymous namespace

TEST_CASE("RotatingAssembly: combines provider, external moment and dynamics in one step")
{
    AeroResult fixed = make_load_result();

    RotatingAssembly assembly(
        [&fixed](const RotationalState&, ModelStatus&) { return fixed; },
        make_constant_moment(ConstantMomentConfig{40.0}),
        make_dynamics());

    RotationalState state;
    state.omega = 10.0;

    const AssemblyStepResult result = assembly.step(0.1, state);

    REQUIRE(result.ok);
    REQUIRE(result.status.ok);
    CHECK(result.status.warnings.empty());

    // net = aero torque 100 - external 40 = 60; alpha = 6.
    CHECK(result.net_moment == doctest::Approx(60.0).epsilon(1e-12));
    CHECK(result.external_moment == doctest::Approx(40.0).epsilon(1e-12));
    CHECK(result.aero.torque == doctest::Approx(100.0).epsilon(1e-12));
    CHECK(result.aero.axial_force == doctest::Approx(50.0).epsilon(1e-12));

    // Powers use the pre-step omega = 10.
    CHECK(result.aero_power == doctest::Approx(1000.0).epsilon(1e-12));
    CHECK(result.mechanical_power == doctest::Approx(600.0).epsilon(1e-12));

    // Heun: omega_new = 10 + 6*0.1 = 10.6; theta = (10 + 10.6)/2 * 0.1 = 1.03.
    CHECK(result.state.omega == doctest::Approx(10.6).epsilon(1e-12));
    CHECK(result.state.angle_rad == doctest::Approx(1.03).epsilon(1e-12));

    REQUIRE(result.per_element.size() == 1);
    CHECK(result.per_element[0].r == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(result.per_element[0].force[1] == doctest::Approx(50.0).epsilon(1e-12));
}

TEST_CASE("RotatingAssembly: provider failure aborts the step")
{
    RotatingAssembly assembly(
        [](const RotationalState&, ModelStatus& status) {
            status.ok = false;
            status.error = "provider exploded";
            return AeroResult{};
        },
        make_constant_moment(ConstantMomentConfig{10.0}),
        make_dynamics());

    RotationalState state;
    state.omega = 2.0;

    const AssemblyStepResult result = assembly.step(0.1, state);

    CHECK_FALSE(result.ok);
    CHECK_FALSE(result.status.ok);
    CHECK(result.status.error == "provider exploded");
}

TEST_CASE("RotatingAssembly: invalid aero moments warn and run with zero torque")
{
    RotatingAssembly assembly(
        [](const RotationalState&, ModelStatus&) {
            AeroResult aero;
            aero.moments.valid = false; // degenerate axis: no usable torque
            return aero;
        },
        make_constant_moment(ConstantMomentConfig{0.0}),
        make_dynamics());

    RotationalState state;
    state.omega = 10.0;

    const AssemblyStepResult result = assembly.step(0.1, state);

    REQUIRE(result.ok);
    REQUIRE(result.status.ok);
    REQUIRE(result.status.warnings.size() == 1);
    CHECK(result.status.warnings[0] == "aero moments invalid (axis?)");

    // Zero effective torque and zero external: no acceleration; omega is
    // unchanged and the Heun angle still advances by omega*dt.
    CHECK(result.net_moment == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.aero_power == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.mechanical_power == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.state.omega == doctest::Approx(10.0).epsilon(1e-12));
    CHECK(result.state.angle_rad == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("RotatingAssembly: state progresses over multiple steps")
{
    RotatingAssembly assembly(
        [](const RotationalState&, ModelStatus&) {
            AeroResult aero;
            aero.moments.valid = true;
            aero.moments.torque = 100.0;
            return aero;
        },
        make_constant_moment(ConstantMomentConfig{40.0}),
        make_dynamics());

    RotationalState state;
    state.omega = 0.0;
    state.angle_rad = 0.0;

    const AssemblyStepResult step1 = assembly.step(0.1, state);
    REQUIRE(step1.ok);
    CHECK(step1.state.omega == doctest::Approx(0.6).epsilon(1e-12));

    const AssemblyStepResult step2 = assembly.step(0.1, step1.state);
    REQUIRE(step2.ok);
    // omega increases by alpha*dt = 6*0.1 = 0.6 each step.
    CHECK(step2.state.omega == doctest::Approx(1.2).epsilon(1e-12));
    // Heun angle accumulation:
    // step1: theta = (0 + 0.6)/2 * 0.1 = 0.03
    // step2: theta = 0.03 + (0.6 + 1.2)/2 * 0.1 = 0.12
    CHECK(step1.state.angle_rad == doctest::Approx(0.03).epsilon(1e-12));
    CHECK(step2.state.angle_rad == doctest::Approx(0.12).epsilon(1e-12));
}