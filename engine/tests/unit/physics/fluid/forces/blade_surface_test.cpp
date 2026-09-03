#include <exd/engine/physics/fluid/forces/force_evaluator.hpp>
#include <exd/engine/physics/rigid_body/rotational_state.hpp>
#include <exd/engine/physics/rigid_body/status.hpp>
#include <doctest/doctest.h>

#include <array>
#include <cmath>

using namespace exd::engine::physics::fluid::forces;
using namespace exd::engine::physics::rigid_body;

namespace {

BladeGeometry make_three_blade()
{
    BladeGeometry b;
    b.r_hub = 1.0;
    b.r_tip = 3.0;
    b.blade_count = 3;
    b.z_rotor = 0.5;

    BladeStation s0;
    s0.r = 1.5;
    s0.dr = 1.0;
    s0.chord = 0.4;
    s0.twist_deg = 0.0;
    s0.thickness_ratio = 0.12;

    BladeStation s1;
    s1.r = 2.5;
    s1.dr = 1.0;
    s1.chord = 0.4;
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

double length(const std::array<double, 3>& v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

} // anonymous namespace

TEST_CASE("blade surface: counts, normals, areas and side separation at zero azimuth")
{
    const auto blade = make_three_blade();
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    // 2 sides × 2 stations × 3 blades = 12 sample points.
    REQUIRE(surface.points.size() == 12u);
    REQUIRE(surface.normals.size() == 12u);
    REQUIRE(surface.areas.size() == 12u);
    REQUIRE(surface.element_index.size() == 12u);

    for (std::size_t i = 0; i < 12; ++i)
    {
        CHECK(length(surface.normals[i]) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(surface.areas[i] == doctest::Approx(0.2).epsilon(1e-12)); // dr·chord/2
        // Zero twist: normals are exactly ±(0,0,1).
        CHECK(std::fabs(surface.normals[i][0]) < 1e-12);
        CHECK(std::fabs(surface.normals[i][1]) < 1e-12);
        CHECK(std::fabs(std::fabs(surface.normals[i][2]) - 1.0) < 1e-12);
        // Blade 0 at θ=0: station 0 mid at (1.5, 0, 0.5).
        if (surface.element_index[i] == 0 && i < 2)
        {
            CHECK(surface.points[i][0] == doctest::Approx(1.5).epsilon(1e-9));
            CHECK(surface.points[i][1] == doctest::Approx(0.0).epsilon(1e-9));
        }
    }

    // Upper/lower pair separation = thickness_ratio·chord = 0.048 along z.
    CHECK(surface.points[0][2] == doctest::Approx(0.5 + 0.024).epsilon(1e-9)); // upper
    CHECK(surface.points[1][2] == doctest::Approx(0.5 - 0.024).epsilon(1e-9)); // lower
    CHECK(surface.points[0][2] - surface.points[1][2] ==
          doctest::Approx(0.048).epsilon(1e-9));
}

TEST_CASE("blade surface: azimuth rotation places blades on their theta_b")
{
    const auto blade = make_three_blade();
    const auto surface = build_blade_surface(blade, M_PI / 2.0, axis_z());

    // Blade 0 at θ = π/2: e_r = (0,1,0) → x ≈ 0, y = r > 0.
    // Blade 1 at θ = π/2 + 2π/3 = 7π/6, blade 2 at 11π/6.
    bool blade0_found = false;
    bool blade1_found = false;
    bool blade2_found = false;
    for (const auto& p : surface.points)
    {
        const double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
        if (r == doctest::Approx(1.5).epsilon(1e-9))
        {
            if (std::fabs(p[0]) < 1e-9 && p[1] > 0.0) blade0_found = true;
            if (p[0] < 0.0 && p[1] < 0.0) blade1_found = true;
            if (p[0] > 0.0 && p[1] < 0.0) blade2_found = true;
        }
    }
    CHECK(blade0_found);
    CHECK(blade1_found);
    CHECK(blade2_found);

    // Explicit azimuth of blade 1 at station radius 1.5:
    // p = 1.5 · e_r(7π/6) = (-0.75·√3, -0.75, ·).
    bool blade1_pos = false;
    for (const auto& p : surface.points)
    {
        const double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
        if (r == doctest::Approx(1.5).epsilon(1e-9) &&
            p[0] == doctest::Approx(1.5 * std::cos(7.0 * M_PI / 6.0)).epsilon(1e-6) &&
            p[1] == doctest::Approx(1.5 * std::sin(7.0 * M_PI / 6.0)).epsilon(1e-6))
        {
            blade1_pos = true;
        }
    }
    CHECK(blade1_pos);
}

TEST_CASE("blade surface: twist rotates the face normal and chord toward +axis")
{
    auto blade = make_three_blade();
    for (auto& st : blade.stations) st.twist_deg = 10.0;
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    // Blade 0 at θ=0: e_t = (0,1,0), e_z = (0,0,1).
    // n_hat = -sin(10°)·e_t + cos(10°)·e_z.
    const double sb = std::sin(10.0 * M_PI / 180.0);
    const double cb = std::cos(10.0 * M_PI / 180.0);
    const auto& n = surface.normals[0]; // blade 0, station 0, upper side

    CHECK(length(n) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(n[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(n[2] == doctest::Approx(cb).epsilon(1e-12));  // normal z-component = cos(10°)
    CHECK(n[1] == doctest::Approx(-sb).epsilon(1e-12)); // in-plane = -sin(10°)

    // Chord direction c_hat = cos(10°)·e_t + sin(10°)·e_z is orthogonal to n_hat
    // in the e_t-e_z plane, so c_hat·e_z = sin(10°) = -n_hat·e_t.
    CHECK(-n[1] == doctest::Approx(sb).epsilon(1e-12));
}

TEST_CASE("blade surface: single-blade geometry still works")
{
    auto blade = make_three_blade();
    blade.blade_count = 1;
    const auto surface = build_blade_surface(blade, 0.0, axis_z());

    REQUIRE(surface.points.size() == 4u); // 2 sides × 2 stations
    CHECK(surface.element_index[0] == 0); // station 0 upper
    CHECK(surface.element_index[1] == 0); // station 0 lower
    CHECK(surface.element_index[2] == 1); // station 1 upper
    CHECK(surface.element_index[3] == 1); // station 1 lower
}

TEST_CASE("blade surface: degenerate axis yields error and empty surface")
{
    const auto blade = make_three_blade();
    RotationAxis degenerate;
    degenerate.origin = {0.0, 0.0, 0.0};
    degenerate.direction = {0.0, 0.0, 0.0};

    ModelStatus status;
    const RotorFrame frame = make_rotor_frame(0.0, degenerate, status);
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(frame.e_z == std::array<double, 3>{0, 0, 1}); // default frame

    const auto surface = build_blade_surface(blade, 0.0, degenerate);
    CHECK(surface.points.empty());
    CHECK(surface.normals.empty());
    CHECK(surface.areas.empty());
    CHECK(surface.element_index.empty());
}