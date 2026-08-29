#include <exd/physics/mechanics/rotational_state.hpp>

#include <doctest/doctest.h>

#include <array>
#include <vector>

using namespace exd::physics::mechanics;

namespace
{

void require_array_approx(const std::array<double, 3>& actual,
                          const std::array<double, 3>& expected)
{
    CHECK(actual[0] == doctest::Approx(expected[0]).epsilon(1e-12));
    CHECK(actual[1] == doctest::Approx(expected[1]).epsilon(1e-12));
    CHECK(actual[2] == doctest::Approx(expected[2]).epsilon(1e-12));
}

ElementForce3D make_force(std::array<double, 3> ref,
                          std::array<double, 3> force,
                          std::array<double, 3> moment = {0.0, 0.0, 0.0})
{
    ElementForce3D f;
    f.ref = ref;
    f.force = force;
    f.moment = moment;
    return f;
}

} // anonymous namespace

TEST_CASE("normalize: divides a non-unit vector and reports success")
{
    std::array<double, 3> v = {0.0, 0.0, 2.0};
    REQUIRE(normalize(v));
    CHECK(v[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(v[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(v[2] == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("normalize: zero vector reports failure and leaves input unchanged")
{
    std::array<double, 3> v = {0.0, 0.0, 0.0};
    CHECK_FALSE(normalize(v));
    CHECK(v[0] == 0.0);
    CHECK(v[1] == 0.0);
    CHECK(v[2] == 0.0);
}

TEST_CASE("integrate_moment: pure axial force yields torque 0 and axial force 10")
{
    RotationAxis axis; // default: origin {0,0,0}, direction {0,0,1}
    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 0.0}, {0.0, 0.0, 10.0}),
    };

    const MomentResult result = integrate_moment(forces, axis);

    REQUIRE(result.valid);
    CHECK(result.torque == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.axial_force == doctest::Approx(10.0).epsilon(1e-12));
    require_array_approx(result.total_force, {0.0, 0.0, 10.0});
    // (ref x force) = {0, -20, 0}
    require_array_approx(result.total_moment, {0.0, -20.0, 0.0});
}

TEST_CASE("integrate_moment: tangential force produces torque 20 about +z")
{
    RotationAxis axis;
    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 0.0}, {0.0, 10.0, 0.0}),
    };

    const MomentResult result = integrate_moment(forces, axis);

    REQUIRE(result.valid);
    CHECK(result.torque == doctest::Approx(20.0).epsilon(1e-12));
    CHECK(result.axial_force == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("integrate_moment: multiple elements sum and element moments add to total moment only")
{
    RotationAxis axis;
    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 5.0}),
        make_force({0.0, 3.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}),
    };

    const MomentResult result = integrate_moment(forces, axis);

    REQUIRE(result.valid);
    // Torque: 20 from element 1, -3 from element 2; element moments (5, 2) excluded.
    CHECK(result.torque == doctest::Approx(17.0).epsilon(1e-12));
    CHECK(result.axial_force == doctest::Approx(0.0).epsilon(1e-12));
    require_array_approx(result.total_force, {1.0, 10.0, 0.0});
    // total_moment includes element moments: (0,0,20)+(0,0,5)+(0,0,-3)+(2,0,0).
    require_array_approx(result.total_moment, {2.0, 0.0, 22.0});
}

TEST_CASE("integrate_moment: non-unit axis direction gives the same scalar moments")
{
    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 0.0}, {0.0, 10.0, 0.0}),
    };

    RotationAxis unit_axis;
    unit_axis.direction = {0.0, 0.0, 1.0};

    RotationAxis scaled_axis;
    scaled_axis.direction = {0.0, 0.0, 2.0};

    const MomentResult unit_result = integrate_moment(forces, unit_axis);
    const MomentResult scaled_result = integrate_moment(forces, scaled_axis);

    REQUIRE(unit_result.valid);
    REQUIRE(scaled_result.valid);
    CHECK(scaled_result.torque == doctest::Approx(unit_result.torque).epsilon(1e-12));
    CHECK(scaled_result.torque == doctest::Approx(20.0).epsilon(1e-12));
    CHECK(scaled_result.axial_force == doctest::Approx(unit_result.axial_force).epsilon(1e-12));
    require_array_approx(scaled_result.total_moment, unit_result.total_moment);
}

TEST_CASE("integrate_moment: degenerate axis direction reports invalid with zeros")
{
    RotationAxis axis;
    axis.direction = {0.0, 0.0, 0.0};

    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 0.0}, {0.0, 10.0, 0.0}),
    };

    const MomentResult result = integrate_moment(forces, axis);

    CHECK_FALSE(result.valid);
    CHECK(result.torque == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.axial_force == doctest::Approx(0.0).epsilon(1e-12));
    require_array_approx(result.total_force, {0.0, 0.0, 0.0});
    require_array_approx(result.total_moment, {0.0, 0.0, 0.0});
}

TEST_CASE("integrate_moment: torque is invariant to the axis origin along the axis")
{
    // The application point, force and rotated origin all move by {0,0,5},
    // so the lever arm (ref - origin) is unchanged.
    std::vector<ElementForce3D> forces = {
        make_force({2.0, 0.0, 5.0}, {0.0, 10.0, 0.0}),
    };

    RotationAxis origin_axis; // origin {0,0,0}
    RotationAxis offset_axis;
    offset_axis.origin = {0.0, 0.0, 5.0};

    const MomentResult origin_result = integrate_moment(forces, origin_axis);
    const MomentResult offset_result = integrate_moment(forces, offset_axis);

    REQUIRE(origin_result.valid);
    REQUIRE(offset_result.valid);
    CHECK(offset_result.torque == doctest::Approx(origin_result.torque).epsilon(1e-12));
    CHECK(offset_result.torque == doctest::Approx(20.0).epsilon(1e-12));
    CHECK(offset_result.axial_force == doctest::Approx(0.0).epsilon(1e-12));
}