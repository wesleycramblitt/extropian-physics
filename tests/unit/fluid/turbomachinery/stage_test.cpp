// stage_test.cpp
// Single-stage mean-line Euler closure: exact Euler-work identities,
// closure consistency checks at a meaningful Mach number, envelope/choke
// warnings, and domain validation.

#include <doctest/doctest.h>

#include <exd/physics/fluid/turbomachinery/stage.hpp>
#include <exd/physics/fluid/turbomachinery/stage_stack.hpp>
#include <exd/physics/thermo/eos.hpp>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using exd::physics::ModelStatus;
using exd::physics::fluid::turbomachinery::StageConfig;
using exd::physics::fluid::turbomachinery::StageGeometryConfig;
using exd::physics::fluid::turbomachinery::StageInlet;
using exd::physics::fluid::turbomachinery::StageResult;
using exd::physics::fluid::turbomachinery::solve_stage;
using exd::physics::fluid::turbomachinery::validate_stage_config;

namespace
{

std::unique_ptr<exd::physics::thermo::IEos> make_air()
{
    return exd::physics::thermo::make_ideal_gas({287.05, 1.4});
}

// Independent reimplementation of the stage inlet density-velocity fixed
// point (the closure's step 2). Used as a cross-check, not by the library.
struct InletState
{
    bool converged = false;
    double c_a = 0.0;
    double c_w1 = 0.0;
    double T1 = 0.0;
    double p1 = 0.0;
    double rho1 = 0.0;
};

InletState inlet_closure(const exd::physics::thermo::IEos& eos,
                         double p0, double T0, double alpha_1,
                         double mdot, double area)
{
    const double cp = eos.specific_heat_cp();
    const double gamma = eos.gamma();
    ModelStatus status;
    const double rho0 = eos.density(p0, T0, status);
    double c_a = mdot / (rho0 * area);
    for (int iter = 0; iter < 200; ++iter)
    {
        const double c_w1 = c_a * std::tan(alpha_1);
        const double c1_sq = c_a * c_a + c_w1 * c_w1;
        const double T1 = T0 - c1_sq / (2.0 * cp);
        const double p1 = p0 * std::pow(T1 / T0, gamma / (gamma - 1.0));
        const double rho1 = eos.density(p1, T1, status);
        const double c_a_new = mdot / (rho1 * area);
        if (std::fabs(c_a_new - c_a) / c_a < 1e-12)
        {
            return InletState{true, c_a_new, c_w1, T1, p1, rho1};
        }
        c_a = c_a_new;
    }
    return InletState{};
}

bool has_substring(const std::vector<std::string>& messages, const std::string& needle)
{
    for (const std::string& message : messages)
    {
        if (message.find(needle) != std::string::npos)
            return true;
    }
    return false;
}

} // namespace

TEST_CASE("stage: Euler work exact on synthetic triangles")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageConfig config;
    config.geometry.r_hub = 0.18;
    config.geometry.r_tip = 0.22;
    config.geometry.alpha_1_rad = 0.30;
    config.geometry.beta_2_rad = 0.55;

    const double r_mean = config.geometry.r_mean();
    const double area = config.geometry.flow_area();
    const double omega = 500.0;
    const double mdot = 2.0;
    const double u = omega * r_mean;

    StageInlet inlet;
    ModelStatus status;
    const StageResult result = solve_stage(config, inlet, omega, mdot, *eos, status);
    REQUIRE(result.ok);
    REQUIRE(status.ok);

    // Independent triangle reconstruction.
    const InletState st = inlet_closure(*eos, inlet.p0, inlet.T0,
                                        config.geometry.alpha_1_rad, mdot, area);
    REQUIRE(st.converged);
    const double c_w1 = st.c_w1;
    const double c_w2 = u + st.c_a * std::tan(config.geometry.beta_2_rad);
    const double delta_c_w = c_w2 - c_w1;

    const double expected_delta_h0 = u * delta_c_w;
    const double expected_torque = mdot * r_mean * delta_c_w;
    const double expected_power = expected_torque * omega;

    CHECK(result.delta_h0 == doctest::Approx(expected_delta_h0).epsilon(1e-9));
    CHECK(result.torque == doctest::Approx(expected_torque).epsilon(1e-9));
    CHECK(result.power == doctest::Approx(expected_power).epsilon(1e-9));
    CHECK(result.c_theta_out == doctest::Approx(c_w2).epsilon(1e-9));
    CHECK(result.work_coefficient == doctest::Approx(expected_delta_h0 / (u * u)).epsilon(1e-9));
}

TEST_CASE("stage: zero swirl change gives zero work")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // With u = 0 (omega = 0 allowed by the closure) and alpha = beta = 0 the
    // absolute swirl change is exactly zero: delta_h0, torque and power
    // vanish and the total-pressure ratio is exactly 1.
    StageConfig config;
    config.geometry.alpha_1_rad = 0.0;
    config.geometry.beta_2_rad = 0.0;

    StageInlet inlet;
    ModelStatus status;
    const StageResult zero = solve_stage(config, inlet, 0.0, 2.0, *eos, status);
    CHECK(zero.ok);
    CHECK(std::fabs(zero.delta_h0) < 1e-12);
    CHECK(zero.pi == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(std::fabs(zero.torque) < 1e-12);
    CHECK(std::fabs(zero.power) < 1e-12);

    // NOTE (spec deviation): "alpha_1 == 0 and beta_2 == 0 -> delta_h0 == 0"
    // is NOT the closure's physics at a nonzero shaft speed: a beta_2 = 0
    // rotor with no inlet swirl adds absolute swirl equal to u, so EITHER
    // u = 0 or beta_2 tuned so c_w2 = c_w1 produces zero work. Record the
    // closure's actual value here so the deviation is explicit.
    const StageResult rotating = solve_stage(config, inlet, 500.0, 2.0, *eos, status);
    CHECK(std::fabs(rotating.delta_h0 - (500.0 * config.geometry.r_mean()) *
                                            (500.0 * config.geometry.r_mean())) < 1e-6);

    // True zero-absolute-swirl-change at a nonzero speed: beta_2 chosen (from
    // the converged inlet) so the rotor re-establishes the inlet swirl.
    StageConfig tuned = config;
    const InletState st = inlet_closure(*eos, inlet.p0, inlet.T0,
                                        config.geometry.alpha_1_rad, 2.0,
                                        config.geometry.flow_area());
    REQUIRE(st.converged);
    const double u = 500.0 * config.geometry.r_mean();
    tuned.geometry.beta_2_rad = std::atan((st.c_w1 - u) / st.c_a);
    const StageResult tuned_r = solve_stage(tuned, inlet, 500.0, 2.0, *eos, status);
    CHECK(tuned_r.ok);
    CHECK(tuned_r.pi == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("stage: geometry-emergent sense reversal (full closure)")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Stage A compresses (+work).
    const double omega = 600.0;
    const double mdot = 2.0;
    const double r_hub = 0.10;
    const double r_tip = 0.15;
    const double u = omega * 0.5 * (r_hub + r_tip);

    StageConfig stage_a;
    stage_a.geometry.r_hub = r_hub;
    stage_a.geometry.r_tip = r_tip;
    stage_a.geometry.alpha_1_rad = 0.35;
    stage_a.geometry.beta_2_rad = 0.55;

    StageInlet inlet;
    ModelStatus status;
    StageResult a = solve_stage(stage_a, inlet, omega, mdot, *eos, status);
    REQUIRE(a.ok);
    REQUIRE(a.delta_h0 > 0.0);

    // Recover stage A's converged inlet triangle (independent loop), then
    // build the exact undo stage B in the ideal shared-c_a frame:
    //   B LE absolute angle  = A exit absolute angle
    //   B TE relative angle  = sets B exit absolute swirl back to A's inlet
    const InletState st = inlet_closure(*eos, inlet.p0, inlet.T0,
                                        stage_a.geometry.alpha_1_rad, mdot,
                                        stage_a.geometry.flow_area());
    REQUIRE(st.converged);
    const double c_w1A = st.c_w1;
    const double c_w2A = u + st.c_a * std::tan(stage_a.geometry.beta_2_rad);

    StageConfig stage_b;
    stage_b.geometry = stage_a.geometry;
    stage_b.geometry.alpha_1_rad = std::atan(c_w2A / st.c_a);
    stage_b.geometry.beta_2_rad = std::atan((c_w1A - u) / st.c_a);

    exd::physics::fluid::turbomachinery::StageStackConfig stack;
    stack.stages = {stage_a, stage_b};
    ModelStatus stack_status;
    auto stack_result =
        exd::physics::fluid::turbomachinery::solve_stage_stack(stack, inlet, omega, mdot, *eos, stack_status);
    REQUIRE(stack_result.ok);
    REQUIRE(stack_result.per_stage.size() == 2);

    const double dhA = stack_result.per_stage[0].delta_h0;
    const double dhB = stack_result.per_stage[1].delta_h0;

    // Work reverses sense. NOTE (spec deviation): exact cancellation to 1e-9
    // is impossible under the full compressible closure because c_a depends on
    // the stage's own inlet state (angle AND total pressure), so a B stage
    // designed from A's c_a has a different converged c_a. The residual here
    // is ~2.2e-3 relative; assertions use a tolerance that admits it.
    CHECK(dhA > 0.0);
    CHECK(dhB < 0.0);
    CHECK(dhB == doctest::Approx(-dhA).epsilon(5e-3));
    CHECK(stack_result.total_delta_h0 == doctest::Approx(0.0).epsilon(1e-2 * dhA));

    // Total enthalpy returns to the inlet to within the same residual (~14
    // mK here; a strict 1e-9 'returns' is not achievable, see above).
    CHECK(std::fabs(stack_result.T0_out - inlet.T0) < 0.05);

    // Polytropes are not reciprocal: enthalpy is recovered but entropy rose,
    // so total pressure drops below the inlet value.
    CHECK(stack_result.total_pi < 1.0);
    ModelStatus sstatus;
    const double s_in = eos->specific_entropy(inlet.p0, inlet.T0, sstatus);
    const double s_out = eos->specific_entropy(stack_result.p0_out, stack_result.T0_out, sstatus);
    CHECK(s_out - s_in > 0.0);

    // Isentropic limit: eta = 1 makes the pair reversible (total_pi ~ 1).
    // The small residual again comes from the density-velocity c_a shift.
    StageConfig stage_a_eta1 = stage_a;
    StageConfig stage_b_eta1 = stage_b;
    stage_a_eta1.loss.polytropic_efficiency = 1.0;
    stage_b_eta1.loss.polytropic_efficiency = 1.0;
    exd::physics::fluid::turbomachinery::StageStackConfig stack_r;
    stack_r.stages = {stage_a_eta1, stage_b_eta1};
    ModelStatus reversible_status;
    auto reversible =
        exd::physics::fluid::turbomachinery::solve_stage_stack(stack_r, inlet, omega, mdot, *eos, reversible_status);
    REQUIRE(reversible.ok);
    CHECK(reversible.total_pi == doctest::Approx(1.0).epsilon(5e-3));
}

TEST_CASE("stage: closed-form closure at meaningful Mach")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Small flow area: A = pi*(0.15^2 - 0.10^2) ~ 0.0393 m^2, r_mean = 0.125.
    StageConfig config;
    config.geometry.r_hub = 0.10;
    config.geometry.r_tip = 0.15;
    config.geometry.alpha_1_rad = 0.45;
    config.geometry.beta_2_rad = 0.25;

    const double cp = eos->specific_heat_cp();
    const double mdot = 4.0;
    const double omega = 900.0;

    StageInlet inlet;
    ModelStatus status;
    const StageResult result = solve_stage(config, inlet, omega, mdot, *eos, status);
    REQUIRE(result.ok);

    const InletState st = inlet_closure(*eos, inlet.p0, inlet.T0,
                                        config.geometry.alpha_1_rad, mdot,
                                        config.geometry.flow_area());
    REQUIRE(st.converged);

    const double c1_sq = st.c_a * st.c_a + st.c_w1 * st.c_w1;
    const double c2_sq = st.c_a * st.c_a + result.c_theta_out * result.c_theta_out;

    // (a) Static/total consistency at both stations.
    CHECK(cp * st.T1 + 0.5 * c1_sq == doctest::Approx(cp * inlet.T0).epsilon(1e-6));
    CHECK(cp * result.static_T + 0.5 * c2_sq == doctest::Approx(cp * result.T0_out).epsilon(1e-6));

    // (b) Continuity. The closure iterates the INLET density-velocity pair,
    // so inlet continuity holds to ~1e-10 relative. The exit static state is
    // purely algebraic, hence the exit density x c_a x A differs by the
    // stage's own compression residual (~1e-3); both are asserted here.
    CHECK(st.rho1 * st.c_a * config.geometry.flow_area() ==
          doctest::Approx(mdot).epsilon(1e-9));
    const double exit_mdot = result.static_rho * st.c_a * config.geometry.flow_area();
    CHECK(exit_mdot == doctest::Approx(mdot).epsilon(2e-3));

    // (c) Energy: cp*(T02 - T01) = delta_h0.
    CHECK(cp * (result.T0_out - inlet.T0) == doctest::Approx(result.delta_h0).epsilon(1e-9));

    // (d) Inside the subsonic validity envelope, positive static temperature.
    CHECK(result.mach_rel_le < 0.7);
    CHECK(result.static_T > 0.0);
    CHECK(result.static_rho > 0.0);
}

TEST_CASE("stage: envelope warnings")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageConfig config;
    config.geometry.r_hub = 0.10;
    config.geometry.r_tip = 0.15;
    config.geometry.alpha_1_rad = 0.45;
    config.geometry.beta_2_rad = 0.25;

    StageInlet inlet;
    const double mdot = 4.0;

    // mach_rel_le ~ 0.84 here: exceeds the 0.7 envelope, still subsonic.
    ModelStatus status;
    StageResult envelope = solve_stage(config, inlet, 2500.0, mdot, *eos, status);
    REQUIRE(envelope.ok);
    CHECK(!envelope.choked);
    CHECK(!status.warnings.empty());
    CHECK(has_substring(status.warnings, "0.7"));

    // Choked: mach_rel_le >> 1 at the rotor LE.
    ModelStatus status_choke;
    StageResult choked = solve_stage(config, inlet, 6000.0, mdot, *eos, status_choke);
    REQUIRE(choked.ok);
    CHECK(choked.choked);
    CHECK(has_substring(status_choke.warnings, "choked"));
}

TEST_CASE("stage: domain validation")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    std::string error;
    std::vector<std::string> warnings;

    // r_hub >= r_tip -> fatal configuration error.
    StageConfig bad_geom;
    bad_geom.geometry.r_hub = 0.30;
    bad_geom.geometry.r_tip = 0.20;
    CHECK_FALSE(validate_stage_config(bad_geom, error, warnings));
    CHECK(error.find("r_hub") != std::string::npos);

    // polytropic_efficiency outside (0, 1] -> fatal configuration error.
    StageConfig bad_eta_lo;
    bad_eta_lo.loss.polytropic_efficiency = 0.0;
    CHECK_FALSE(validate_stage_config(bad_eta_lo, error, warnings));
    CHECK(error.find("polytropic_efficiency") != std::string::npos);

    StageConfig bad_eta_hi;
    bad_eta_hi.loss.polytropic_efficiency = 1.5;
    CHECK_FALSE(validate_stage_config(bad_eta_hi, error, warnings));

    // Non-fatal envelope warning for a strongly tapered annulus config.
    StageConfig tapered;
    tapered.geometry.r_hub = 0.15;
    tapered.geometry.r_tip = 0.40;
    CHECK(validate_stage_config(tapered, error, warnings));
    CHECK_FALSE(warnings.empty());

    // Operating-point domain: mdot <= 0 fails through the status channel.
    StageConfig ok_config;
    StageInlet inlet;
    ModelStatus status;
    StageResult bad = solve_stage(ok_config, inlet, 500.0, 0.0, *eos, status);
    CHECK(!bad.ok);
    CHECK(!status.ok);
    CHECK(status.error.find("mdot") != std::string::npos);

    // Same operating-point checks for omega < 0 and non-positive inlet state.
    ModelStatus status_neg;
    StageResult bad_omega = solve_stage(ok_config, inlet, -1.0, 2.0, *eos, status_neg);
    CHECK(!bad_omega.ok);

    StageInlet bad_inlet;
    bad_inlet.T0 = 0.0;
    ModelStatus status_inlet;
    StageResult bad_T = solve_stage(ok_config, bad_inlet, 500.0, 2.0, *eos, status_inlet);
    CHECK(!bad_T.ok);
}