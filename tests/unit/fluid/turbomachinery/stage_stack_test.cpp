// stage_stack_test.cpp
// Multi-stage propagation: single-stage equivalence, the two-stage reheat
// effect (inter-stage T0 rise degrades the second stage's pressure ratio at
// equal work), and bit-exact determinism.

#include <doctest/doctest.h>

#include <exd/physics/fluid/turbomachinery/stage.hpp>
#include <exd/physics/fluid/turbomachinery/stage_stack.hpp>
#include <exd/physics/thermo/eos.hpp>

#include <cmath>
#include <memory>
#include <vector>

using exd::physics::ModelStatus;
using exd::physics::fluid::turbomachinery::StageConfig;
using exd::physics::fluid::turbomachinery::StageInlet;
using exd::physics::fluid::turbomachinery::StageResult;
using exd::physics::fluid::turbomachinery::StageStackConfig;
using exd::physics::fluid::turbomachinery::StageStackResult;
using exd::physics::fluid::turbomachinery::solve_stage;
using exd::physics::fluid::turbomachinery::solve_stage_stack;
using exd::physics::fluid::turbomachinery::validate_stage_stack_config;

namespace
{

std::unique_ptr<exd::physics::thermo::IEos> make_air()
{
    return exd::physics::thermo::make_ideal_gas({287.05, 1.4});
}

StageStackConfig default_single_stage()
{
    StageStackConfig stack;
    stack.stages.resize(1);
    return stack;
}

} // namespace

TEST_CASE("stage_stack: single-stage stack equals single stage")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    const double omega = 500.0;
    const double mdot = 2.0;

    StageStackConfig stack = default_single_stage();
    StageInlet inlet;
    ModelStatus stack_status;
    const StageStackResult stack_result =
        solve_stage_stack(stack, inlet, omega, mdot, *eos, stack_status);
    REQUIRE(stack_result.ok);
    REQUIRE(stack_result.per_stage.size() == 1);

    ModelStatus stage_status;
    const StageResult stage_result =
        solve_stage(stack.stages.front(), inlet, omega, mdot, *eos, stage_status);
    REQUIRE(stage_result.ok);

    const StageResult& s1 = stack_result.per_stage.front();
    CHECK(s1.p0_out == doctest::Approx(stage_result.p0_out).epsilon(1e-9));
    CHECK(s1.T0_out == doctest::Approx(stage_result.T0_out).epsilon(1e-9));
    CHECK(s1.c_theta_out == doctest::Approx(stage_result.c_theta_out).epsilon(1e-9));
    CHECK(s1.static_p == doctest::Approx(stage_result.static_p).epsilon(1e-9));
    CHECK(s1.static_T == doctest::Approx(stage_result.static_T).epsilon(1e-9));
    CHECK(s1.static_rho == doctest::Approx(stage_result.static_rho).epsilon(1e-9));
    CHECK(s1.delta_h0 == doctest::Approx(stage_result.delta_h0).epsilon(1e-9));
    CHECK(s1.torque == doctest::Approx(stage_result.torque).epsilon(1e-9));
    CHECK(s1.power == doctest::Approx(stage_result.power).epsilon(1e-9));
    CHECK(s1.pi == doctest::Approx(stage_result.pi).epsilon(1e-9));
    CHECK(s1.tau == doctest::Approx(stage_result.tau).epsilon(1e-9));
    CHECK(s1.mach_rel_le == doctest::Approx(stage_result.mach_rel_le).epsilon(1e-9));
    CHECK(s1.work_coefficient == doctest::Approx(stage_result.work_coefficient).epsilon(1e-9));
    CHECK(s1.flow_coefficient == doctest::Approx(stage_result.flow_coefficient).epsilon(1e-9));
    CHECK(s1.choked == stage_result.choked);

    CHECK(stack_result.p0_out == doctest::Approx(stage_result.p0_out).epsilon(1e-9));
    CHECK(stack_result.T0_out == doctest::Approx(stage_result.T0_out).epsilon(1e-9));
    CHECK(stack_result.c_theta_out == doctest::Approx(stage_result.c_theta_out).epsilon(1e-9));
    CHECK(stack_result.total_pi == doctest::Approx(stage_result.pi).epsilon(1e-9));
    CHECK(stack_result.total_delta_h0 == doctest::Approx(stage_result.delta_h0).epsilon(1e-9));
    CHECK(stack_result.total_torque == doctest::Approx(stage_result.torque).epsilon(1e-9));
    CHECK(stack_result.total_power == doctest::Approx(stage_result.power).epsilon(1e-9));
    CHECK(stack_result.mach_rel_max == doctest::Approx(stage_result.mach_rel_le).epsilon(1e-9));
}

TEST_CASE("stage_stack: two-stage reheat effect")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Identical default stages. With alpha_1 == beta_2 the Euler work per
    // stage is exactly u^2 (the c_a terms cancel), so the stack sum is exact.
    StageStackConfig stack;
    stack.stages.resize(2);
    const double omega = 500.0;
    const double mdot = 2.0;
    const double u = omega * stack.stages.front().geometry.r_mean();

    StageInlet inlet;
    ModelStatus status;
    const StageStackResult result = solve_stage_stack(stack, inlet, omega, mdot, *eos, status);
    REQUIRE(result.ok);
    REQUIRE(result.per_stage.size() == 2);

    const StageResult& stage1 = result.per_stage[0];
    const StageResult& stage2 = result.per_stage[1];

    CHECK(stage1.delta_h0 == doctest::Approx(u * u).epsilon(1e-9));
    CHECK(stage2.delta_h0 == doctest::Approx(u * u).epsilon(1e-9));
    CHECK(result.total_delta_h0 == doctest::Approx(2.0 * u * u).epsilon(1e-9));

    // Reheat effect: equal work on a hotter inlet yields tau_2 < tau_1 and
    // therefore pi_2 < pi_1; the two-stage total pressure ratio is below the
    // square of the single-stage ratio.
    CHECK(stage2.tau < stage1.tau);
    CHECK(stage2.pi < stage1.pi);
    CHECK(result.total_pi < stage1.pi * stage1.pi);

    // Torque/power sum on one shaft.
    CHECK(result.total_torque == doctest::Approx(stage1.torque + stage2.torque).epsilon(1e-9));
    CHECK(result.total_power == doctest::Approx(stage1.power + stage2.power).epsilon(1e-9));
    CHECK(result.total_power == doctest::Approx(result.total_torque * omega).epsilon(1e-9));
}

TEST_CASE("stage_stack: determinism")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageStackConfig stack;
    stack.stages.resize(3);
    stack.stages[1].geometry.r_hub = 0.10;
    stack.stages[1].geometry.r_tip = 0.15;
    stack.stages[1].geometry.alpha_1_rad = 0.45;
    stack.stages[1].geometry.beta_2_rad = 0.25;

    StageInlet inlet;
    const double omega = 700.0;
    const double mdot = 2.5;

    ModelStatus status_a;
    const StageStackResult a = solve_stage_stack(stack, inlet, omega, mdot, *eos, status_a);
    ModelStatus status_b;
    const StageStackResult b = solve_stage_stack(stack, inlet, omega, mdot, *eos, status_b);

    REQUIRE(a.ok);
    REQUIRE(b.ok);
    REQUIRE(a.per_stage.size() == b.per_stage.size());

    CHECK(a.p0_out == b.p0_out);
    CHECK(a.T0_out == b.T0_out);
    CHECK(a.c_theta_out == b.c_theta_out);
    CHECK(a.total_pi == b.total_pi);
    CHECK(a.total_delta_h0 == b.total_delta_h0);
    CHECK(a.total_torque == b.total_torque);
    CHECK(a.total_power == b.total_power);
    CHECK(a.mach_rel_max == b.mach_rel_max);
    for (size_t i = 0; i < a.per_stage.size(); ++i)
    {
        const StageResult& x = a.per_stage[i];
        const StageResult& y = b.per_stage[i];
        CHECK(x.p0_out == y.p0_out);
        CHECK(x.T0_out == y.T0_out);
        CHECK(x.delta_h0 == y.delta_h0);
        CHECK(x.torque == y.torque);
        CHECK(x.pi == y.pi);
        CHECK(x.mach_rel_le == y.mach_rel_le);
        CHECK(!std::isnan(x.delta_h0));
        CHECK(!std::isnan(x.T0_out));
    }
}

TEST_CASE("stage_stack: empty stack validation")
{
    std::string error;
    std::vector<std::string> warnings;
    StageStackConfig empty;
    CHECK_FALSE(validate_stage_stack_config(empty, error, warnings));
    CHECK(error.find("at least one") != std::string::npos);

    StageStackConfig bad;
    bad.stages.resize(1);
    bad.stages[0].geometry.r_hub = 0.3;
    bad.stages[0].geometry.r_tip = 0.2;
    CHECK_FALSE(validate_stage_stack_config(bad, error, warnings));
    CHECK(error.find("stage 0") != std::string::npos);

    StageStackConfig ok_stack = default_single_stage();
    CHECK(validate_stage_stack_config(ok_stack, error, warnings));
}