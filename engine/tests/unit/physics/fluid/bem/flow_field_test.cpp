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

// Locate grid indices relative to the rotor plane.  z is ordered upstream
// (high z) -> downstream (low z), so the nearest upstream point has the
// largest index with z > z_r, and the nearest downstream point has the
// smallest index with z < z_r.
struct GridIndices
{
    std::size_t upstream = 0;
    std::size_t disk = 0;
    std::size_t downstream = 0;
};

GridIndices locate_indices(const FlowFieldGrid& grid, double z_r)
{
    GridIndices idx;
    const std::size_t nz = grid.z.size();

    // Locate the rotor plane on the grid.  The solver places z_r at a grid
    // station, so the upstream/downstream neighbors are immediately adjacent.
    double best = std::fabs(grid.z[0] - z_r);
    idx.disk = 0;
    for (std::size_t i = 1; i < nz; ++i)
    {
        const double d = std::fabs(grid.z[i] - z_r);
        if (d < best)
        {
            best = d;
            idx.disk = i;
        }
    }

    idx.upstream = (idx.disk > 0) ? idx.disk - 1 : 0;
    idx.downstream = (idx.disk + 1 < nz) ? idx.disk + 1 : nz - 1;
    return idx;
}

} // namespace

TEST_CASE("flow_field: grid sizes and R_max match specification")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);
    REQUIRE(res.converged);

    const auto& g = res.flow_field;
    CHECK(g.z.size() == config.field_axial_points);
    CHECK(g.r.size() == config.field_radial_points);
    CHECK(g.velocity.size() == g.z.size() * g.r.size());
    CHECK(g.pressure.size() == g.z.size() * g.r.size());
    CHECK(g.r.back() >= res.radial.back().radius_m);
    CHECK(g.r.back() >= 1.0);
}

TEST_CASE("flow_field: upstream station nearest z_r gives V ≈ V_rotor inside duct")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const auto idx = locate_indices(res.flow_field, res.duct.v_rotor > 0.0 ? 0.7 : 0.7);
    // z_r for this fixture is 0.825, but locate just needs a consistent z_r.
    // Recompute exact z_r from geometry: 0.5*(0.5*(0.8+0.8) + 0.5*(0.864+0.836)) = 0.825.
    const auto gidx = locate_indices(res.flow_field, 0.825);
    const std::size_t i_up = gidx.upstream;
    const std::size_t nr = res.flow_field.r.size();

    // Inside the cylindrical duct at r ≈ 0.5 (outside the hub).
    const std::size_t j = nr / 2;
    const double V = res.flow_field.velocity[i_up * nr + j];
    CHECK(V == doctest::Approx(res.duct.v_rotor).epsilon(1e-3));
}

TEST_CASE("flow_field: disk plane matches BEM velocity at element radii")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const auto gidx = locate_indices(res.flow_field, 0.825);
    const std::size_t i_disk = gidx.disk;
    const std::size_t nr = res.flow_field.r.size();

    const double r_hub = res.radial.front().radius_m - 0.5 * (res.radial[1].radius_m - res.radial[0].radius_m);
    const double r_tip = res.radial.back().radius_m + 0.5 * (res.radial.back().radius_m - res.radial[res.radial.size() - 2].radius_m);
    const double dr_grid = (nr > 1) ? (res.flow_field.r.back() - res.flow_field.r.front()) / (nr - 1) : 0.0;

    for (const auto& rs : res.radial)
    {
        // Find the nearest radial grid point to this element radius.
        std::size_t j_best = 0;
        double best = std::fabs(res.flow_field.r[0] - rs.radius_m);
        for (std::size_t j = 1; j < nr; ++j)
        {
            const double d = std::fabs(res.flow_field.r[j] - rs.radius_m);
            if (d < best)
            {
                best = d;
                j_best = j;
            }
        }
        const double r_grid = res.flow_field.r[j_best];
        // Skip grid points that sit in the boundary layer where a(r) is
        // blended to zero at the hub or tip.
        if (r_grid < r_hub + dr_grid || r_grid > r_tip - dr_grid)
            continue;

        const double V = res.flow_field.velocity[i_disk * nr + j_best];
        const double expected = res.duct.v_rotor * (1.0 - rs.induction_axial);
        // The disk plane is placed at z_r and the radial grid does not align
        // exactly with element midpoints, so allow a generous relative tolerance.
        CHECK(V == doctest::Approx(expected).epsilon(0.15));
    }
}

TEST_CASE("flow_field: pressure jump across disk is present and negative")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const auto gidx = locate_indices(res.flow_field, 0.825);
    const std::size_t i_up = gidx.upstream;
    const std::size_t i_down = gidx.downstream;
    const std::size_t nr = res.flow_field.r.size();

    // Pick an in-annulus radius near mid-span.
    const std::size_t j = res.flow_field.r.size() / 2;
    const double r = res.flow_field.r[j];
    REQUIRE(r > res.radial.front().radius_m);
    REQUIRE(r < res.radial.back().radius_m);

    const double p_up = res.flow_field.pressure[i_up * nr + j];
    const double p_down = res.flow_field.pressure[i_down * nr + j];
    const double delta_p = p_down - p_up;

    // Interpolate the per-element pressure jump to this radius.
    double dp_elem = 0.0;
    for (std::size_t k = 0; k + 1 < res.radial.size(); ++k)
    {
        const auto& a = res.radial[k];
        const auto& b = res.radial[k + 1];
        if (r >= a.radius_m && r <= b.radius_m)
        {
            const double t = (r - a.radius_m) / (b.radius_m - a.radius_m);
            dp_elem = a.pressure_jump + t * (b.pressure_jump - a.pressure_jump);
            break;
        }
    }
    REQUIRE(dp_elem > 0.0);

    // The pressure drop is negative; allow for the Bernoulli kinetic term
    // across the disk (both sides use the freestream reference).
    CHECK(delta_p < 0.0);
    CHECK(std::fabs(delta_p + dp_elem) / dp_elem < 0.85);
}

TEST_CASE("flow_field: wake deficit grows downstream and stays within [a, 2a]")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const auto& g = res.flow_field;
    const std::size_t nz = g.z.size();
    const std::size_t nr = g.r.size();

    // Use a fixed in-annulus radius for the "centerline" sample (hub blocks r=0).
    const std::size_t j_center = static_cast<std::size_t>(0.5 * nr); // r ≈ 0.5
    const double r_center = g.r[j_center];

    const auto gidx = locate_indices(g, 0.825);
    const std::size_t i_disk = gidx.disk;
    const std::size_t i_far = nz - 1;

    const double V_disk = g.velocity[i_disk * nr + j_center];
    const double V_far = g.velocity[i_far * nr + j_center];
    const double V_rotor = res.duct.v_rotor;

    // Interpolate induction a(r) to the centerline radius.
    double a_center = 0.0;
    for (std::size_t k = 0; k + 1 < res.radial.size(); ++k)
    {
        const auto& ra = res.radial[k];
        const auto& rb = res.radial[k + 1];
        if (r_center >= ra.radius_m && r_center <= rb.radius_m)
        {
            const double t = (r_center - ra.radius_m) / (rb.radius_m - ra.radius_m);
            a_center = ra.induction_axial + t * (rb.induction_axial - ra.induction_axial);
            break;
        }
    }
    REQUIRE(a_center > 0.0);

    // Disk deficit = a*V_rotor; far-wake centerline deficit ≤ 2a*V_rotor.
    CHECK(V_rotor - V_far >= (V_rotor - V_disk) * 0.99);
    CHECK(V_rotor - V_far <= 2.0 * a_center * V_rotor * 1.01);

    // Radial profile at a fixed downstream station: velocity recovers toward
    // the annulus edge / freestream (outside the loaded disk a(r) = 0).
    const std::size_t i_mid = (i_disk + i_far) / 2;
    const std::size_t j_edge = nr - 1; // r = R_max, a = 0
    const double V_center = g.velocity[i_mid * nr + j_center];
    const double V_edge   = g.velocity[i_mid * nr + j_edge];
    CHECK(V_edge > V_center);

    // Wake spreads: centerline deficit is larger farther downstream.
    CHECK(V_far < V_center);
}

TEST_CASE("flow_field: far radial velocity returns to freestream")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    const auto& g = res.flow_field;
    const std::size_t nz = g.z.size();
    const std::size_t nr = g.r.size();

    // Check a few axial stations at the outermost radial station.
    for (std::size_t i : {std::size_t{0}, nz / 2, nz - 1})
    {
        const double V = g.velocity[i * nr + (nr - 1)];
        CHECK(V == doctest::Approx(cond.v_inf).epsilon(1e-3));
    }
}

TEST_CASE("flow_field: no NaN anywhere in the grid")
{
    PolarDatabase polars;
    polars.add(make_synthetic_polar());

    BEMSolverConfig config;
    config.element_count = 32;
    config.k_duct = 0.0;
    config.reference_area = ReferenceArea::RotorDisk;
    config.include_flow_field = true;
    config.field_radial_points = 64;
    config.airfoils = {{0.0, "ideal"}, {1.0, "ideal"}};

    OperatingConditions cond;
    cond.v_inf = 10.0;
    cond.rho = 1.225;
    cond.mu = 1.81e-5;

    auto res = solve_turbine(make_turbine(), cond, polars, config);
    REQUIRE(res.valid);

    for (double v : res.flow_field.velocity)
        CHECK(std::isfinite(v));
    for (double p : res.flow_field.pressure)
        CHECK(std::isfinite(p));
}
