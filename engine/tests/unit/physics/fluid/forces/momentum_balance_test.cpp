#include <exd/engine/physics/fluid/forces/force_evaluator.hpp>
#include <exd/engine/physics/rigid_body/rotational_state.hpp>
#include <exd/engine/physics/rigid_body/status.hpp>
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

using namespace exd::engine::physics::fluid::forces;
using namespace exd::engine::physics::rigid_body;

namespace {

BladeGeometry make_turbine_blade()
{
    BladeGeometry b;
    b.r_hub = 0.5;
    b.r_tip = 5.0;
    b.blade_count = 3;
    b.z_rotor = 0.0;

    const int n = 20;
    const double dr = (b.r_tip - b.r_hub) / static_cast<double>(n);
    for (int i = 0; i < n; ++i)
    {
        BladeStation st;
        st.r = b.r_hub + (static_cast<double>(i) + 0.5) * dr;
        st.dr = dr;
        st.chord = 0.4;
        st.twist_deg = -2.0;
        st.thickness_ratio = 0.12;
        st.airfoil = "naca0012";
        b.stations.push_back(st);
    }
    return b;
}

RotationAxis axis_z()
{
    RotationAxis axis;
    axis.origin = {0.0, 0.0, 0.0};
    axis.direction = {0.0, 0.0, 1.0};
    return axis;
}

SurfaceFlow make_freestream()
{
    SurfaceFlow f;
    f.points = {{0.0, 0.0, 0.0}};
    f.normals = {{0.0, 0.0, 1.0}};
    f.velocity = {{0.0, 0.0, -10.0}}; // uniform inflow, velocity[0] is read
    f.shear_traction = {{0.0, 0.0, 0.0}};
    f.pressure = {101325.0};
    f.area = {1.0};
    f.element_index = {0};
    f.density = 1.225;
    f.viscosity = 1.81e-5;
    f.p_ref = 101325.0;
    return f;
}

std::unique_ptr<IForceEvaluator> make_momentum_evaluator(const MomentumBalanceConfig& config,
                                                         ModelStatus& status)
{
    ForceEvaluatorParams params;
    params.type = ForceEvaluatorType::MomentumBalance;
    params.momentum = config;
    params.polars = nullptr; // use built-in polars
    return make_force_evaluator(params, status);
}

} // anonymous namespace

TEST_CASE("momentum balance: power-extracting turbine sense at TSR ~ 3")
{
    const auto blade = make_turbine_blade();
    const auto flow = make_freestream();

    ModelStatus status;
    auto eval = make_momentum_evaluator(MomentumBalanceConfig{}, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    const double omega = 6.0; // rad/s → TSR ≈ 6·5/10 = 3
    eval->compute(blade, flow, omega, axis_z(), per_element, status);
    REQUIRE(status.ok);
    REQUIRE(per_element.size() == 60u); // 20 stations × 3 blades

    const int B = blade.blade_count;
    for (std::size_t i = 0; i < per_element.size(); ++i)
    {
        // Station-major ordering: r repeats every blade_count entries.
        CHECK(per_element[i].r ==
              doctest::Approx(blade.stations[i / static_cast<std::size_t>(B)].r).epsilon(1e-12));
        CHECK(std::isfinite(per_element[i].force[0]));
        CHECK(std::isfinite(per_element[i].force[1]));
        CHECK(std::isfinite(per_element[i].force[2]));
        CHECK(std::isfinite(per_element[i].moment[0]));
    }

    const auto result = integrate_moment(per_element, axis_z());
    REQUIRE(result.valid);

    // Positive torque about +z (power-extracting sense with omega > 0).
    CHECK(result.torque > 0.0);

    // Thrust acts into the inflow: with inflow (0,0,-10), the axial force is
    // negative (Σ F·e_z = -T_total). Sanity bound |T| < 0.5·rho·v²·π·R².
    const double rho = 1.225;
    const double v = 10.0;
    const double R = 5.0;
    const double A = M_PI * R * R;
    const double thrust_bound = 0.5 * rho * v * v * A;
    CHECK(result.axial_force < 0.0);
    CHECK(std::fabs(result.axial_force) > 0.0);
    CHECK(std::fabs(result.axial_force) < thrust_bound);

    // Cp = Q·omega / (0.5·rho·v³·π·R²) in a loose physical band.
    const double cp = result.torque * omega / (0.5 * rho * v * v * v * A);
    CHECK(cp > 0.05);
    CHECK(cp < 0.593);

    // Pinned from the reference implementation (built-in naca0012 polars).
    CHECK(result.torque == doctest::Approx(487.94).epsilon(0.01));
}

TEST_CASE("momentum balance: extreme conditions stay finite and non-fatal")
{
    const auto blade = make_turbine_blade();
    const auto flow = make_freestream();

    MomentumBalanceConfig config;
    config.max_iterations = 2; // exhaust the iteration budget → non-convergence warning
    config.under_relaxation = 0.25;

    ModelStatus status;
    auto eval = make_momentum_evaluator(config, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, flow, 0.01, axis_z(), per_element, status);
    // Warnings are non-fatal: status stays ok and values stay finite.
    REQUIRE(status.ok);
    REQUIRE(per_element.size() == 60u);
    for (const auto& f : per_element)
    {
        CHECK(std::isfinite(f.force[0]));
        CHECK(std::isfinite(f.force[1]));
        CHECK(std::isfinite(f.force[2]));
        CHECK(std::isfinite(f.moment[0]));
        CHECK(std::isfinite(f.moment[1]));
        CHECK(std::isfinite(f.moment[2]));
    }
}