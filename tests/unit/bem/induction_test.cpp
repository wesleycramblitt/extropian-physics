#include <doctest/doctest.h>

#include <exd/geometry/turbine.hpp>
#include <exd/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>

#include <cmath>
#include <iostream>
#include <string>

using namespace exd::geometry;
using namespace exd::physics::fluid::reduced_order::bem;

namespace {

FlowPath cylindrical_shroud()
{
    FlowPath f;
    f.hub_points    = {{0.0f, 0.2f}, {2.0f, 0.2f}};
    f.shroud_points = {{0.0f, 1.0f}, {2.0f, 1.0f}};
    f.tip_clearance = {0.0f, 0.0f, 0.02f, "m", false};
    return f;
}

BladeRow rotor_betz()
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

AirfoilPolar make_ideal_polar()
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

} // namespace

TEST_CASE("induction: high-loading case engages Buhl branch")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    BladeRow r = rotor_betz();
    r.blade_count.value = 4;   // higher solidity -> a > 0.4
    r.rotational_speed.value = 800; // lower rpm -> higher loading
    turbine.blade_rows = {r};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.under_relaxation = 0.35;
    config.glauert_threshold = 0.4;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    bool buhl_engaged = false;
    for (const auto& rs : res.radial)
    {
        CHECK(std::isfinite(rs.induction_axial));
        CHECK(std::isfinite(rs.induction_tangential));
        CHECK(rs.induction_axial < 1.0);
        CHECK(rs.converged);
        if (rs.induction_axial > 0.4) buhl_engaged = true;
    }
    CHECK(buhl_engaged);
}

TEST_CASE("induction: light-loading case converges with small a")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    BladeRow r = rotor_betz();
    r.rotational_speed.value = 2000.0f; // high rpm -> light loading
    turbine.blade_rows = {r};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 16;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);
    for (const auto& rs : res.radial)
    {
        CHECK(rs.induction_axial >= 0.0);
        CHECK(rs.induction_axial < 0.75);
        CHECK(rs.converged);
    }
}

TEST_CASE("induction: tangential denom guard does not crash")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    BladeRow r = rotor_betz();
    r.sections[0].stagger.value = 85.0f;
    r.sections[1].stagger.value = 85.0f;
    turbine.blade_rows = {r};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 8;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    for (const auto& rs : res.radial)
    {
        CHECK(std::isfinite(rs.induction_tangential));
        CHECK(rs.induction_tangential >= 0.0);
        CHECK(rs.induction_tangential <= 1.0);
    }
}

TEST_CASE("induction: clamping warnings appear for deep-stall-like polar")
{
    AirfoilPolar p;
    p.name = "negative";
    p.re = 1e6;
    for (int i = -10; i <= 10; ++i)
    {
        p.alpha_deg.push_back(static_cast<double>(i));
        p.cl.push_back(-0.5); // negative Cn will clamp axial induction to 0
        p.cd.push_back(0.01);
    }

    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    BladeRow r = rotor_betz();
    turbine.blade_rows = {r};

    PolarDatabase polars;
    polars.add(p);

    BEMSolverConfig config;
    config.element_count = 8;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "negative"}, {1.0, "negative"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    bool clamp_warned = false;
    for (const auto& w : res.warnings)
        if (w.find("clamped") != std::string::npos) clamp_warned = true;
    CHECK(clamp_warned);
}
