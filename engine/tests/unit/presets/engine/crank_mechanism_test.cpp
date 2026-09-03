// crank_mechanism_test.cpp
// Analytic slider-crank kinematics + equivalent inertia.

#include <doctest/doctest.h>

#include <cmath>

#include "../../../../src/presets/engine/engine_internal.hpp"

using namespace exd::engine::presets::engine;

namespace
{
constexpr double PI = 3.14159265358979323846;

EngineGeometryConfig fixture()
{
    EngineGeometryConfig g;
    g.crank_radius = 0.05;
    g.rod_length = 0.20;
    g.bore = 0.086;
    g.clearance_volume = 1.0e-4;
    g.piston_mass = 0.5;
    g.flywheel_inertia = 0.02;
    return g;
}
} // anonymous namespace

TEST_CASE("crank: TDC and BDC extremes")
{
    auto g = fixture();

    // TDC: theta = 0 → x = l + r, velocity factor zero.
    auto k0 = crank_kinematics(0.0, g);
    CHECK(k0.x == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(k0.dx_dtheta == doctest::Approx(0.0).epsilon(1e-9));

    // BDC: theta = pi → x = l − r, velocity factor zero.
    auto kp = crank_kinematics(PI, g);
    CHECK(kp.x == doctest::Approx(0.15).epsilon(1e-9));
    CHECK(kp.dx_dtheta == doctest::Approx(0.0).epsilon(1e-9));

    // Mid-stroke: x = sqrt(l² − r²) since cos(pi/2)=0, moving down (dx/dθ < 0).
    auto kq = crank_kinematics(0.5 * PI, g);
    CHECK(kq.x == doctest::Approx(std::sqrt(0.2 * 0.2 - 0.05 * 0.05)).epsilon(1e-9));
    CHECK(kq.dx_dtheta < 0.0);

    // Volume: TDC = clearance, BDC = clearance + A·stroke.
    const double A = piston_area(g);
    const double stroke = 2.0 * g.crank_radius;
    CHECK(cylinder_volume(0.0, g) == doctest::Approx(g.clearance_volume).epsilon(1e-9));
    CHECK(cylinder_volume(PI, g) == doctest::Approx(g.clearance_volume + A * stroke).epsilon(1e-6));
}

TEST_CASE("crank: equivalent inertia is flywheel at TDC/BDC and larger mid-stroke")
{
    auto g = fixture();
    CHECK(crank_kinematics(0.0, g).j_eq == doctest::Approx(g.flywheel_inertia).epsilon(1e-9));
    CHECK(crank_kinematics(PI, g).j_eq == doctest::Approx(g.flywheel_inertia).epsilon(1e-9));
    CHECK(crank_kinematics(0.5 * PI, g).j_eq > g.flywheel_inertia);
    // dJ/dθ is odd in theta (mirror through TDC).
    const double a = crank_kinematics(0.25 * PI, g).dj_dtheta;
    const double b = crank_kinematics(-0.25 * PI, g).dj_dtheta;
    CHECK(a == doctest::Approx(-b).epsilon(1e-9));
}

TEST_CASE("crank: dJ/dθ matches finite difference of J_eq")
{
    auto g = fixture();
    const double theta = 0.7;
    const double h = 1e-5;
    const double jp = crank_kinematics(theta + h, g).j_eq;
    const double jm = crank_kinematics(theta - h, g).j_eq;
    CHECK(crank_kinematics(theta, g).dj_dtheta == doctest::Approx((jp - jm) / (2 * h)).epsilon(1e-3));
}

TEST_CASE("crank: d²x/dθ² matches finite difference of dx/dθ")
{
    auto g = fixture();
    const double theta = 1.1;
    const double h = 1e-5;
    const double dp = crank_kinematics(theta + h, g).dx_dtheta;
    const double dm = crank_kinematics(theta - h, g).dx_dtheta;
    CHECK(crank_kinematics(theta, g).d2x_dtheta2 == doctest::Approx((dp - dm) / (2 * h)).epsilon(1e-3));
}