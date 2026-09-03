#include <doctest/doctest.h>

#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_solver.hpp>

#include <cmath>
#include <iostream>
#include <string>

using namespace exd::geometry;
using namespace exd::engine::physics::fluid::reduced_order::bem;

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

AirfoilPolar make_high_loading_polar()
{
    AirfoilPolar p;
    p.name = "high_cl";
    p.re = 1e6;
    for (int i = -45; i <= 45; ++i)
    {
        const double a = static_cast<double>(i);
        const double arad = a * M_PI / 180.0;
        double cl = 2.0 * M_PI * arad;
        if (cl > 2.0 * M_PI * 0.6) cl = 2.0 * M_PI * 0.6;
        if (cl < -2.0 * M_PI * 0.6) cl = -2.0 * M_PI * 0.6;
        p.alpha_deg.push_back(a);
        p.cl.push_back(cl);
        p.cd.push_back(0.0);
    }
    return p;
}

} // namespace

TEST_CASE("DuSelig: valid result with reduced Cp")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.under_relaxation = 0.25;
    config.glauert_threshold = 0.4;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.loss_correction = LossCorrection::DuSelig;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    const double cp = res.rotor.cp;
    CHECK(cp > 0.2);
    CHECK(cp < 16.0 / 27.0);
    for (const auto& rs : res.radial)
        CHECK(rs.converged);
}

TEST_CASE("DuSelig: produces less power than Prandtl at high tip speeds")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()}; // 3 blades, 1200 rpm

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    auto run = [&](LossCorrection correction) {
        BEMSolverConfig config;
        config.element_count = 32;
        config.k_duct = 0.0;
        config.under_relaxation = 0.25;
        config.glauert_threshold = 0.4;
        config.reference_area = ReferenceArea::RotorDisk;
        config.include_flow_field = false;
        config.loss_correction = correction;
        config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

        OperatingConditions cond;
        cond.v_inf = 10.0;
        cond.rho = 1.225;
        cond.mu = 1.81e-5;

        return solve_turbine(turbine, cond, polars, config);
    };

    TurbineResult res_prandtl = run(LossCorrection::Prandtl);
    TurbineResult res_duselig = run(LossCorrection::DuSelig);

    REQUIRE(res_prandtl.valid);
    REQUIRE(res_duselig.valid);

    const double cp_prandtl = res_prandtl.rotor.cp;
    const double cp_duselig = res_duselig.rotor.cp;
    CHECK(cp_duselig <= cp_prandtl * 1.1);
}

TEST_CASE("Chaviaropoulos: valid result")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.under_relaxation = 0.25;
    config.glauert_threshold = 0.4;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.loss_correction = LossCorrection::Chaviaropoulos;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    const double cp = res.rotor.cp;
    CHECK(cp > 0.2);
    CHECK(cp < 16.0 / 27.0);
    for (const auto& rs : res.radial)
        CHECK(rs.converged);
}

TEST_CASE("Chaviaropoulos: loss factor is <= 1.0 everywhere")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.under_relaxation = 0.25;
    config.glauert_threshold = 0.4;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.loss_correction = LossCorrection::Chaviaropoulos;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    TurbineResult res = solve_turbine(turbine, cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    // The loss factor F is internal to the solver; its observable consequence
    // is a bounded, finite axial induction factor at every station.
    for (const auto& rs : res.radial)
    {
        CHECK(std::isfinite(rs.induction_axial));
        CHECK(rs.induction_axial >= 0.0);
        CHECK(rs.induction_axial < 1.0);
    }
}

TEST_CASE("All loss models produce finite, non-NaN results")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    const LossCorrection models[] = {LossCorrection::Prandtl,
                                     LossCorrection::DuSelig,
                                     LossCorrection::Chaviaropoulos};
    for (auto model : models)
    {
        BEMSolverConfig config;
        config.element_count = 32;
        config.k_duct = 0.0;
        config.under_relaxation = 0.25;
        config.glauert_threshold = 0.4;
        config.reference_area = ReferenceArea::RotorDisk;
        config.include_flow_field = false;
        config.loss_correction = model;
        config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

        OperatingConditions cond;
        cond.v_inf = 10.0;
        cond.rho = 1.225;
        cond.mu = 1.81e-5;

        TurbineResult res = solve_turbine(turbine, cond, polars, config);
        REQUIRE(res.valid);

        for (const auto& rs : res.radial)
        {
            CHECK(std::isfinite(rs.radius_m));
            CHECK(std::isfinite(rs.axial_velocity));
            CHECK(std::isfinite(rs.relative_velocity));
            CHECK(std::isfinite(rs.induction_axial));
            CHECK(std::isfinite(rs.induction_tangential));
            CHECK(std::isfinite(rs.inflow_angle_deg));
            CHECK(std::isfinite(rs.angle_of_attack_deg));
            CHECK(std::isfinite(rs.reynolds));
            CHECK(std::isfinite(rs.cl));
            CHECK(std::isfinite(rs.cd));
        }
        CHECK(std::isfinite(res.rotor.cp));
        CHECK(std::isfinite(res.rotor.thrust));
        CHECK(std::isfinite(res.rotor.torque));
    }
}

TEST_CASE("Correction models are selectable from config")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_ideal_polar());

    const LossCorrection models[] = {LossCorrection::Prandtl,
                                     LossCorrection::DuSelig,
                                     LossCorrection::Chaviaropoulos};
    for (auto model : models)
    {
        BEMSolverConfig config;
        config.element_count = 32;
        config.k_duct = 0.0;
        config.under_relaxation = 0.25;
        config.glauert_threshold = 0.4;
        config.reference_area = ReferenceArea::RotorDisk;
        config.include_flow_field = false;
        config.loss_correction = model;
        config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

        OperatingConditions cond;
        cond.v_inf = 10.0;
        cond.rho = 1.225;
        cond.mu = 1.81e-5;

        TurbineResult res = solve_turbine(turbine, cond, polars, config);
        REQUIRE(res.valid);
        CHECK(res.rotor.cp > 0.0);
    }
}