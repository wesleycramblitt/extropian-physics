#include <doctest/doctest.h>

#include <exd/geometry/turbine.hpp>
#include <exd/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>

#include <cmath>
#include <iostream>

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
    // Tapered, twisted rotor: keeps local angle of attack and induction moderate.
    r.leading_edge_hub    = {0.8f, 0.2f};
    r.leading_edge_shroud = {0.8f, 1.0f};
    r.trailing_edge_hub   = {0.864f, 0.2f};   // chord_hub ~0.064
    r.trailing_edge_shroud= {0.836f, 1.0f};   // chord_tip ~0.036
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

} // namespace

TEST_CASE("BEM Betz-limit test: synthetic ideal polar, cylindrical shroud")
{
    TurbineDefinition turbine;
    turbine.flow_path = cylindrical_shroud();
    turbine.blade_rows = {rotor_betz()};

    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;       // no duct acceleration
    config.under_relaxation = 0.25;
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

    for (const auto& rs : res.radial)
    {
        CHECK(rs.induction_axial > 0.0);
        CHECK(rs.induction_axial < 0.7); // tip loss drives the outermost stations higher
        CHECK(rs.converged);
    }

    const double cp = res.rotor.cp;
    CHECK(cp > 0.3);
    CHECK(cp < (16.0 / 27.0) * 0.999);

    CHECK(res.rotor.thrust > 0.0);
    CHECK(res.rotor.torque > 0.0);
    CHECK(res.rotor.power > 0.0);

    // Pinned golden values from first validated run (config above).
    std::cout << "Betz run: Cp=" << cp
              << " T=" << res.rotor.thrust
              << " Q=" << res.rotor.torque
              << " P=" << res.rotor.power << "\n";

    CHECK(cp == doctest::Approx(0.467292).epsilon(1e-3)); // pinned from first validated run
    CHECK(res.rotor.thrust == doctest::Approx(127.451).epsilon(1e-2)); // pinned
    CHECK(res.rotor.power == doctest::Approx(899.176).epsilon(1e-2)); // pinned

    for (const auto& w : res.warnings)
        std::cout << "WARNING: " << w << "\n";
    REQUIRE(res.warnings.empty());
}
