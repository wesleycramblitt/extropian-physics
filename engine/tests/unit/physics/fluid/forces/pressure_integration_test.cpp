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

BladeGeometry make_single_blade()
{
    BladeGeometry b;
    b.r_hub = 1.0;
    b.r_tip = 2.0;
    b.blade_count = 1;
    b.z_rotor = 0.0;

    BladeStation st;
    st.r = 1.0;
    st.dr = 1.0;
    st.chord = 1.0;
    st.twist_deg = 0.0;
    st.thickness_ratio = 0.12;
    b.stations = {st};
    return b;
}

BladeGeometry make_two_station_blade()
{
    BladeGeometry b;
    b.r_hub = 1.0;
    b.r_tip = 3.0;
    b.blade_count = 1;
    b.z_rotor = 0.0;

    BladeStation s0;
    s0.r = 1.0;
    s0.dr = 1.0;
    s0.chord = 1.0;
    s0.twist_deg = 0.0;
    s0.thickness_ratio = 0.12;

    BladeStation s1;
    s1.r = 2.0;
    s1.dr = 1.0;
    s1.chord = 1.0;
    s1.twist_deg = 0.0;
    s1.thickness_ratio = 0.12;

    b.stations = {s0, s1};
    return b;
}

RotationAxis axis_z()
{
    RotationAxis axis;
    axis.origin = {0.0, 0.0, 0.0};
    axis.direction = {0.0, 0.0, 1.0};
    return axis;
}

SurfaceFlow make_flow(const BladeSurface& surface,
                      const std::vector<double>& pressures,
                      double p_ref,
                      const std::array<double, 3>& shear = {0, 0, 0})
{
    SurfaceFlow f;
    f.points = surface.points;
    f.normals = surface.normals;
    f.area = surface.areas;
    f.element_index = surface.element_index;
    const std::size_t n = f.points.size();
    f.velocity.assign(n, std::array<double, 3>{0, 0, -10});
    f.shear_traction.assign(n, shear);
    f.pressure = pressures;
    f.density = 1.225;
    f.viscosity = 1.81e-5;
    f.p_ref = p_ref;
    return f;
}

std::unique_ptr<IForceEvaluator> make_pressure_evaluator(bool include_shear,
                                                         ModelStatus& status)
{
    ForceEvaluatorParams params;
    params.type = ForceEvaluatorType::PressureIntegration;
    params.pressure.include_shear = include_shear;
    return make_force_evaluator(params, status);
}

} // anonymous namespace

TEST_CASE("pressure integration: uniform pressure produces zero force")
{
    const auto blade = make_single_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());
    REQUIRE(surface.points.size() == 2u);

    const double p = 101325.0;
    auto flow = make_flow(surface, {p, p}, p); // p_i == p_ref everywhere

    ModelStatus status;
    auto eval = make_pressure_evaluator(true, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
    REQUIRE(status.ok);
    REQUIRE(per_element.size() == 1u);
    CHECK(per_element[0].force[0] == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(per_element[0].force[1] == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(per_element[0].force[2] == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("pressure integration: pressure jump yields delta-p force")
{
    const auto blade = make_single_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    const double p_ref = 101325.0;
    // Upper p = p_ref + 100 (normal +e_z), lower p = p_ref (normal -e_z):
    // net = 100·(area=0.5)·e_z + 0·0.5·(-e_z) = 50·e_z.
    auto flow = make_flow(surface, {p_ref + 100.0, p_ref}, p_ref);

    ModelStatus status;
    auto eval = make_pressure_evaluator(true, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
    REQUIRE(status.ok);
    REQUIRE(per_element.size() == 1u);
    CHECK(per_element[0].force[0] == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(per_element[0].force[1] == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(per_element[0].force[2] == doctest::Approx(50.0).epsilon(1e-9));
    // Diagnostic split: the whole force is pressure.
    CHECK(per_element[0].force_pressure[2] == doctest::Approx(50.0).epsilon(1e-9));
    CHECK(per_element[0].force_shear[2] == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("pressure integration: wall shear contributes only when enabled")
{
    const auto blade = make_single_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    const double p_ref = 101325.0;
    const std::array<double, 3> shear = {3.0, 0.0, 0.0};
    auto flow = make_flow(surface, {p_ref, p_ref}, p_ref, shear);

    // include_shear = true: 3.0 Pa on both sides × 0.5 m² each → 3.0 N on x.
    {
        ModelStatus status;
        auto eval = make_pressure_evaluator(true, status);
        REQUIRE(status.ok);
        std::vector<ElementForce3D> per_element;
        eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
        REQUIRE(status.ok);
        REQUIRE(per_element.size() == 1u);
        CHECK(per_element[0].force[0] == doctest::Approx(3.0).epsilon(1e-12));
        CHECK(per_element[0].force[1] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(per_element[0].force[2] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(per_element[0].force_shear[0] == doctest::Approx(3.0).epsilon(1e-12));
        CHECK(per_element[0].force_pressure[0] == doctest::Approx(0.0).epsilon(1e-12));
    }

    // include_shear = false: no shear contribution at all.
    {
        ModelStatus status;
        auto eval = make_pressure_evaluator(false, status);
        REQUIRE(status.ok);
        std::vector<ElementForce3D> per_element;
        eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
        REQUIRE(status.ok);
        REQUIRE(per_element.size() == 1u);
        CHECK(per_element[0].force[0] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(per_element[0].force_shear[0] == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("pressure integration: two stations produce independent forces")
{
    const auto blade = make_two_station_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());
    REQUIRE(surface.points.size() == 4u);

    const double p_ref = 101325.0;
    // Point order: station 0 upper, station 0 lower, station 1 upper, station 1 lower.
    std::vector<double> pressures;
    for (std::size_t i = 0; i < surface.points.size(); ++i)
    {
        double p = p_ref;
        if (surface.element_index[i] == 0 && i % 2 == 0) p = p_ref + 100.0;
        if (surface.element_index[i] == 1 && i % 2 == 0) p = p_ref - 50.0;
        pressures.push_back(p);
    }
    auto flow = make_flow(surface, pressures, p_ref);

    ModelStatus status;
    auto eval = make_pressure_evaluator(true, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
    REQUIRE(status.ok);
    REQUIRE(per_element.size() == 2u);
    // Station 0: Δp = 100 on the upper side → +50 N.
    CHECK(per_element[0].r == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(per_element[0].force[2] == doctest::Approx(50.0).epsilon(1e-9));
    // Station 1: Δp = -50 on the upper side → -25 N.
    CHECK(per_element[1].r == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(per_element[1].force[2] == doctest::Approx(-25.0).epsilon(1e-9));
}

TEST_CASE("pressure integration: invalid SurfaceFlow sets error status")
{
    const auto blade = make_single_blade();

    SurfaceFlow bad;
    bad.points = {{0, 0, 0}, {1, 0, 0}};
    bad.normals = {{0, 0, 1}}; // mismatched length → valid() == false
    bad.velocity = {{0, 0, -10}};
    bad.shear_traction = {{0, 0, 0}};
    bad.pressure = {101325.0};
    bad.area = {1.0};
    bad.element_index = {0};
    CHECK_FALSE(bad.valid());

    ModelStatus status;
    auto eval = make_pressure_evaluator(true, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, bad, 0.0, axis_z(), per_element, status);
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(per_element.empty());
}

TEST_CASE("pressure integration: element_index out of range sets error status")
{
    const auto blade = make_single_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    auto flow = make_flow(surface, {101325.0, 101325.0}, 101325.0);
    flow.element_index = {0, 7}; // station 7 does not exist

    ModelStatus status;
    auto eval = make_pressure_evaluator(true, status);
    REQUIRE(status.ok);

    std::vector<ElementForce3D> per_element;
    eval->compute(blade, flow, 0.0, axis_z(), per_element, status);
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(per_element.empty());
}