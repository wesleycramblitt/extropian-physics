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

std::string data_dir()
{
    return std::string(EXT_PHYSICS_DATA_DIR) + "/airfoils";
}

// Note: the canonical annulus fixture from blade_geometry_test has very high
// solidity (24 blades, chord 0.4) and does not converge cleanly with the
// built-in polars.  We use the canonical Betz-limit rotor geometry instead,
// which is also a canonical fixture from the other tests and gives a clean
// converged solve with the built-in NACA4412 polar.
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

TurbineDefinition make_turbine()
{
    TurbineDefinition t;
    t.flow_path = cylindrical_shroud();
    t.blade_rows = {rotor_betz()};
    return t;
}

OperatingConditions base_conditions()
{
    OperatingConditions o;
    o.v_inf = 10.0;
    o.rho   = 1.225;
    o.mu    = 1.81e-5;
    // Tuned for the built-in NACA4412 polar: moderate positive angles of attack
    // across the span, clean convergence, no warnings.
    o.rpm_override = 400.0;
    return o;
}

BEMSolverConfig base_config()
{
    BEMSolverConfig c;
    c.element_count = 32;
    c.k_duct = 0.0;
    c.under_relaxation = 0.25;
    c.reference_area = ReferenceArea::RotorDisk;
    c.include_flow_field = false;
    c.airfoils = {{0.0, "naca4412"}, {1.0, "naca4412"}};
    return c;
}

} // namespace

TEST_CASE("integration: full solve with builtin polars is valid and converged")
{
    PolarDatabase polars;
    polars.add_builtin_polars();

    auto res = solve_turbine(make_turbine(), base_conditions(), polars, base_config());
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    for (const auto& w : res.warnings)
        std::cout << "integration warning: " << w << "\n";
    CHECK(res.warnings.empty());
    CHECK(res.rotor.cp > 0.0);
    CHECK(std::isfinite(res.rotor.cp));
}

TEST_CASE("integration: CSV polar database matches builtin database")
{
    PolarDatabase builtin;
    builtin.add_builtin_polars();

    PolarDatabase csv;
    REQUIRE(csv.load_directory(data_dir()));

    auto res_builtin = solve_turbine(make_turbine(), base_conditions(), builtin, base_config());
    auto res_csv     = solve_turbine(make_turbine(), base_conditions(), csv, base_config());

    REQUIRE(res_builtin.valid);
    REQUIRE(res_csv.valid);
    REQUIRE(res_builtin.converged);
    REQUIRE(res_csv.converged);

    CHECK(res_builtin.rotor.cp == doctest::Approx(res_csv.rotor.cp).epsilon(1e-9));
    CHECK(res_builtin.rotor.ct == doctest::Approx(res_csv.rotor.ct).epsilon(1e-9));
    CHECK(res_builtin.rotor.thrust == doctest::Approx(res_csv.rotor.thrust).epsilon(1e-9));
    CHECK(res_builtin.rotor.power == doctest::Approx(res_csv.rotor.power).epsilon(1e-9));
}

TEST_CASE("integration: element_count sweep shows converging Cp")
{
    PolarDatabase polars;
    polars.add_builtin_polars();

    std::vector<double> cps;
    for (uint32_t n : {16, 32, 64, 128})
    {
        auto config = base_config();
        config.element_count = n;
        auto res = solve_turbine(make_turbine(), base_conditions(), polars, config);
        REQUIRE(res.valid);
        REQUIRE(res.converged);
        cps.push_back(res.rotor.cp);
        std::cout << "Cp(n=" << n << ") = " << res.rotor.cp << "\n";
    }

    CHECK(std::fabs(cps[2] - cps[3]) < 0.01 * std::max(std::fabs(cps[3]), 1e-6));
}
