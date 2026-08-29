#include <exd/physics/mechanics/quaternion.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>

using namespace exd::physics::mechanics;

namespace
{

std::array<double, 4> unit_quat(double w, double x, double y, double z)
{
    const double length = std::sqrt(w * w + x * x + y * y + z * z);
    return {
        w / length,
        x / length,
        y / length,
        z / length,
    };
}

void require_vector4_approx(const std::array<double, 4>& actual,
                            const std::array<double, 4>& expected,
                            double epsilon)
{
    CHECK(actual[0] == doctest::Approx(expected[0]).epsilon(epsilon));
    CHECK(actual[1] == doctest::Approx(expected[1]).epsilon(epsilon));
    CHECK(actual[2] == doctest::Approx(expected[2]).epsilon(epsilon));
    CHECK(actual[3] == doctest::Approx(expected[3]).epsilon(epsilon));
}

void require_vector3_approx(const std::array<double, 3>& actual,
                            const std::array<double, 3>& expected,
                            double epsilon)
{
    CHECK(actual[0] == doctest::Approx(expected[0]).epsilon(epsilon));
    CHECK(actual[1] == doctest::Approx(expected[1]).epsilon(epsilon));
    CHECK(actual[2] == doctest::Approx(expected[2]).epsilon(epsilon));
}

} // anonymous namespace

TEST_CASE("quat_multiply: identity is the multiplicative identity")
{
    const std::array<double, 4> q = unit_quat(0.8, 0.3, -0.4, 0.5);
    const std::array<double, 4> identity = {1.0, 0.0, 0.0, 0.0};

    require_vector4_approx(quat_multiply(q, identity), q, 1e-12);
    require_vector4_approx(quat_multiply(identity, q), q, 1e-12);
}

TEST_CASE("quat_rotate: identity orientation leaves a vector unchanged")
{
    const std::array<double, 3> v = {1.0, 2.0, -3.0};
    require_vector3_approx(quat_rotate({1.0, 0.0, 0.0, 0.0}, v), v, 1e-12);
}

TEST_CASE("quat_rotate: 90 degrees about z maps x to y and y to -x")
{
    // Half-angle representation: q = (cos 45°, 0, 0, sin 45°).
    const double half = std::sqrt(0.5);
    const std::array<double, 4> q = {half, 0.0, 0.0, half};

    require_vector3_approx(quat_rotate(q, {1.0, 0.0, 0.0}), {0.0, 1.0, 0.0}, 1e-12);
    require_vector3_approx(quat_rotate(q, {0.0, 1.0, 0.0}), {-1.0, 0.0, 0.0}, 1e-12);
    // The rotation axis is invariant.
    require_vector3_approx(quat_rotate(q, {0.0, 0.0, 1.0}), {0.0, 0.0, 1.0}, 1e-12);
}

TEST_CASE("quat_rotate: 180 degrees about z maps x to -x")
{
    // q = cos(180°/2) + z·sin(180°/2) = (0, 0, 0, 1).
    require_vector3_approx(quat_rotate({0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}),
                          {-1.0, 0.0, 0.0}, 1e-12);
}

TEST_CASE("quat_normalize: divides by length and leaves the zero quaternion unchanged")
{
    std::array<double, 4> q = {2.0, 0.0, 0.0, 0.0};
    quat_normalize(q);
    require_vector4_approx(q, {1.0, 0.0, 0.0, 0.0}, 1e-12);

    std::array<double, 4> zero = {0.0, 0.0, 0.0, 0.0};
    quat_normalize(zero);
    require_vector4_approx(zero, {0.0, 0.0, 0.0, 0.0}, 1e-12);
}

TEST_CASE("quat_multiply: is associative on unit quaternions")
{
    const std::array<double, 4> a = unit_quat(1.0, 0.2, -0.3, 0.4);
    const std::array<double, 4> b = unit_quat(0.8, -0.1, 0.5, 0.2);
    const std::array<double, 4> c = unit_quat(0.5, 0.6, 0.3, -0.5);

    const std::array<double, 4> left = quat_multiply(quat_multiply(a, b), c);
    const std::array<double, 4> right = quat_multiply(a, quat_multiply(b, c));

    require_vector4_approx(left, right, 1e-12);
}

TEST_CASE("quat_integrate_body: constant body-frame spin rotates about the spin axis")
{
    // Identity orientation with ω = (0, 0, 2) rad/s for dt = 0.1 s should
    // rotate the body x-axis toward +y by about 2·0.1 = 0.2 rad.
    const std::array<double, 4> q =
        quat_integrate_body({1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 2.0}, 0.1);

    // The returned quaternion is unit.
    const double norm =
        std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    CHECK(norm == doctest::Approx(1.0).epsilon(1e-12));

    // First-order integration is accurate to O(dt²): ε = 1e-3 is ample
    // for the ~6.6e-4 rotation-angle error of this step size.
    const std::array<double, 3> world = quat_rotate(q, {1.0, 0.0, 0.0});
    CHECK(world[0] == doctest::Approx(std::cos(0.2)).epsilon(1e-3));
    CHECK(world[1] == doctest::Approx(std::sin(0.2)).epsilon(1e-3));
    CHECK(world[2] == doctest::Approx(0.0).epsilon(1e-12));
}