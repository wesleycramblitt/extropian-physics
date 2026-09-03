#include <doctest/doctest.h>

#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_solver.hpp>

#include <cmath>

using namespace exd::geometry;
using namespace exd::engine::physics::fluid::reduced_order::bem;

namespace {

FlowPath make_flow_path(const std::vector<exd::math::Vec2f>& shroud_points)
{
    FlowPath f;
    f.hub_points    = {{0.0f, 0.2f}, {2.0f, 0.2f}};
    f.shroud_points = shroud_points;
    f.tip_clearance = {0.0f, 0.0f, 0.02f, "m", false};
    return f;
}

BladeRow rotor()
{
    BladeRow r;
    r.type = BladeRowType::Rotor;
    r.blade_count = {3, 1, 200, "", false};
    r.rotational_speed = {1200, 1, 100000, "rpm", false};
    r.leading_edge_hub    = {0.8f, 0.2f};
    r.leading_edge_shroud = {0.8f, 1.0f};
    r.trailing_edge_hub   = {0.864f, 0.2f};
    r.trailing_edge_shroud= {0.836f, 1.0f};
    r.sections = {BladeSection{0.0f}, BladeSection{1.0f}};
    r.sections[0].stagger = {9.0f, -90.0f, 90.0f, "deg", false};
    r.sections[1].stagger = {-1.0f, -90.0f, 90.0f, "deg", false};
    r.tip_feature = TipFeature::None;
    return r;
}

AirfoilPolar make_synthetic_polar()
{
    AirfoilPolar p;
    p.name = "ideal";
    p.re = 1e6;
    for (int i = -45; i <= 45; ++i)
    {
        const double a = static_cast<double>(i);
        const double arad = a * M_PI / 180.0;
        double cl = 2.0 * M_PI * arad;
        if (cl > 2.0 * M_PI * 0.5) cl = 2.0 * M_PI * 0.5;
        if (cl < -2.0 * M_PI * 0.5) cl = -2.0 * M_PI * 0.5;
        p.alpha_deg.push_back(a);
        p.cl.push_back(cl);
        p.cd.push_back(0.0);
    }
    return p;
}

TurbineDefinition make_turbine(const std::vector<exd::math::Vec2f>& shroud_points)
{
    TurbineDefinition t;
    t.flow_path = make_flow_path(shroud_points);
    t.blade_rows = {rotor()};
    return t;
}

BEMSolverConfig base_config()
{
    BEMSolverConfig c;
    c.element_count = 32;
    c.k_duct = 0.0;
    c.under_relaxation = 0.25;
    c.reference_area = ReferenceArea::RotorDisk;
    c.include_flow_field = false;
    c.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};
    return c;
}

OperatingConditions base_conditions()
{
    OperatingConditions o;
    o.v_inf = 10.0;
    o.rho = 1.225;
    o.mu = 1.81e-5;
    return o;
}

} // namespace

TEST_CASE("duct: K_duct=0 gives M_duct=1 exactly")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    auto config = base_config();
    auto res = solve_turbine(make_turbine({{0.0f, 1.0f}, {2.0f, 1.0f}}),
                             base_conditions(), polars, config);
    REQUIRE(res.valid);
    CHECK(res.duct.m_duct == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(res.duct.v_rotor == doctest::Approx(base_conditions().v_inf).epsilon(1e-12));
}

TEST_CASE("duct: K_duct=1 gives M_duct = A_u/A_r")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    auto config = base_config();
    config.k_duct = 1.0;

    // Converging shroud: upstream radius 1.2, rotor radius 1.0.
    auto res = solve_turbine(make_turbine({{0.0f, 1.2f}, {0.825f, 1.0f}, {2.0f, 1.0f}}),
                             base_conditions(), polars, config);
    REQUIRE(res.valid);
    const double area_ratio = (1.2 * 1.2) / (1.0 * 1.0);
    CHECK(res.duct.m_duct == doctest::Approx(area_ratio).epsilon(1e-6));
    CHECK(res.duct.area_ratio == doctest::Approx(area_ratio).epsilon(1e-6));
}

TEST_CASE("duct: converging shroud raises Cp vs cylindrical")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    auto config = base_config();
    config.k_duct = 0.5;

    auto res_cyl = solve_turbine(make_turbine({{0.0f, 1.0f}, {2.0f, 1.0f}}),
                                 base_conditions(), polars, config);
    auto res_conv = solve_turbine(make_turbine({{0.0f, 1.2f}, {0.825f, 1.0f}, {2.0f, 1.0f}}),
                                  base_conditions(), polars, config);
    REQUIRE(res_cyl.valid);
    REQUIRE(res_conv.valid);
    REQUIRE(res_cyl.converged);
    REQUIRE(res_conv.converged);

    CHECK(res_conv.rotor.cp > res_cyl.rotor.cp * 1.01);
}

TEST_CASE("duct: diverging shroud lowers Cp vs cylindrical")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    auto config = base_config();
    config.k_duct = 0.5;

    auto res_cyl = solve_turbine(make_turbine({{0.0f, 1.0f}, {2.0f, 1.0f}}),
                                 base_conditions(), polars, config);
    auto res_div = solve_turbine(make_turbine({{0.0f, 0.8f}, {0.825f, 1.0f}, {2.0f, 1.0f}}),
                                 base_conditions(), polars, config);
    REQUIRE(res_cyl.valid);
    REQUIRE(res_div.valid);
    REQUIRE(res_cyl.converged);
    REQUIRE(res_div.converged);

    CHECK(res_div.rotor.cp < res_cyl.rotor.cp * 0.99);
}

TEST_CASE("duct: pathological diverging geometry gives invalid result")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    auto config = base_config();
    config.k_duct = 1.0;

    // A_u/A_r = (0.1/1.0)^2 = 0.01 -> M_duct = 0.01 <= 0.02.
    auto res = solve_turbine(make_turbine({{0.0f, 0.1f}, {0.825f, 1.0f}, {2.0f, 1.0f}}),
                             base_conditions(), polars, config);
    CHECK_FALSE(res.valid);
    CHECK(res.error.find("M_duct") != std::string::npos);
}
