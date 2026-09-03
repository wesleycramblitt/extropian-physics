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

TurbineDefinition make_turbine()
{
    TurbineDefinition t;
    t.flow_path = cylindrical_shroud();
    t.blade_rows = {rotor_betz()};
    return t;
}

} // namespace

TEST_CASE("hull: hull_cd=0 gives zero drag and net_thrust == rotor thrust")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.hull_cd = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);
    CHECK(res.hull.drag == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.hull.pressure_drag == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.hull.viscous_drag == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.system.net_thrust == doctest::Approx(res.rotor.thrust).epsilon(1e-9));
}

TEST_CASE("hull: hull_cd=0.5 drag matches exact frontal-area formula")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.hull_cd = 0.5;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const double A_f = M_PI * 1.0 * 1.0;
    const double expected_drag = 0.5 * cond.rho * cond.v_inf * cond.v_inf * config.hull_cd * A_f;
    CHECK(res.hull.drag == doctest::Approx(expected_drag).epsilon(1e-9));
    CHECK(res.hull.reference_area == doctest::Approx(A_f).epsilon(1e-9));
}

TEST_CASE("hull: viscous drag is positive and pressure drag is non-negative")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.hull_cd = 0.2;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);
    CHECK(res.hull.viscous_drag > 0.0);
    CHECK(res.hull.pressure_drag >= 0.0);
    CHECK(res.hull.drag == doctest::Approx(res.hull.pressure_drag + res.hull.viscous_drag).epsilon(1e-9));
}

TEST_CASE("hull: net thrust equals rotor thrust minus hull drag")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.hull_cd = 0.2;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = false;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);
    CHECK(res.system.net_thrust == doctest::Approx(res.rotor.thrust - res.hull.drag).epsilon(1e-9));
}
