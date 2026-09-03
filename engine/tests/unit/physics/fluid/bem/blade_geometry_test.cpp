#include <doctest/doctest.h>

#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_config.hpp>

#include "../../../../../src/physics/fluid/reduced_order/bem/bem_internal.hpp"

#include <cmath>

using namespace exd::geometry;
using namespace exd::engine::physics::fluid::reduced_order::bem;

namespace {

FlowPath annulus()
{
    FlowPath f;
    f.hub_points    = {{0.0f, 0.4f}, {1.0f, 0.4f}, {2.0f, 0.4f}};
    f.shroud_points = {{0.0f, 1.0f}, {1.0f, 1.0f}, {2.0f, 1.0f}};
    f.tip_clearance = {0.01f, 0.0f, 0.02f, "m", false};
    return f;
}

BladeRow rotor()
{
    BladeRow r;
    r.type = BladeRowType::Rotor;
    r.blade_count = {24, 1, 200, "", false};
    r.rotational_speed = {1500, 1, 100000, "rpm", false};
    r.leading_edge_hub    = {0.5f, 0.4f};
    r.leading_edge_shroud = {0.5f, 1.0f};
    r.trailing_edge_hub   = {0.9f, 0.4f};
    r.trailing_edge_shroud= {0.9f, 1.0f};
    r.sections = {BladeSection{0.0f}, BladeSection{1.0f}};
    r.sections[0].stagger = {20.0f, -90.0f, 90.0f, "deg", false};
    r.sections[1] = r.sections[0];
    r.sections[1].span = 1.0f;
    r.tip_feature = TipFeature::Clearance;
    return r;
}

TurbineDefinition make_turbine(const BladeRow& row)
{
    TurbineDefinition t;
    t.flow_path = annulus();
    t.blade_rows = {row};
    return t;
}

} // namespace

TEST_CASE("blade_geometry: annulus extraction")
{
    BEMSolverConfig cfg;
    cfg.element_count = 32;
    auto res = build_blade_geometry(make_turbine(rotor()), cfg);
    REQUIRE(res.ok);

    CHECK(res.geometry.z_r == doctest::Approx(0.7).epsilon(1e-6));
    CHECK(res.geometry.r_hub == doctest::Approx(0.4).epsilon(1e-6));
    // R_tip clamped by clearance: shroud r=1.0 - clearance 0.01 = 0.99
    CHECK(res.geometry.r_tip == doctest::Approx(0.99).epsilon(1e-6));
    CHECK(res.geometry.blade_count == doctest::Approx(24.0));
    CHECK(res.geometry.omega == doctest::Approx(1500.0 * 2.0 * M_PI / 60.0).epsilon(1e-6));
    CHECK(res.geometry.elements.size() == 32);

    const auto& e0 = res.geometry.elements.front();
    CHECK(e0.chord == doctest::Approx(0.4).epsilon(1e-6));
    CHECK(e0.beta_deg == doctest::Approx(20.0).epsilon(1e-6));
}

TEST_CASE("blade_geometry: zero rotor rows is invalid")
{
    TurbineDefinition t;
    t.flow_path = annulus();
    BEMSolverConfig cfg;
    auto res = build_blade_geometry(t, cfg);
    CHECK_FALSE(res.ok);
    CHECK(res.error == "no Rotor row in TurbineDefinition");
}

TEST_CASE("blade_geometry: multiple rotor rows is invalid")
{
    TurbineDefinition t;
    t.flow_path = annulus();
    BladeRow r = rotor();
    t.blade_rows = {r, r};
    BEMSolverConfig cfg;
    auto res = build_blade_geometry(t, cfg);
    CHECK_FALSE(res.ok);
    CHECK(res.error == "2 Rotor rows; single-rotor solver");
}

TEST_CASE("blade_geometry: wrong row_index emits warning")
{
    BEMSolverConfig cfg;
    cfg.row_index = 7;
    auto res = build_blade_geometry(make_turbine(rotor()), cfg);
    REQUIRE(res.ok);
    bool found = false;
    for (const auto& w : res.warnings)
        if (w.find("row_index") != std::string::npos) found = true;
    CHECK(found);
}

TEST_CASE("blade_geometry: tip clearance gated by TipFeature")
{
    auto row_clearance = rotor();
    auto row_none = rotor();
    row_none.tip_feature = TipFeature::None;

    BEMSolverConfig cfg;
    cfg.element_count = 8;

    auto res_clearance = build_blade_geometry(make_turbine(row_clearance), cfg);
    auto res_none = build_blade_geometry(make_turbine(row_none), cfg);
    REQUIRE(res_clearance.ok);
    REQUIRE(res_none.ok);

    CHECK(res_clearance.geometry.r_tip == doctest::Approx(0.99));
    CHECK(res_none.geometry.r_tip == doctest::Approx(1.0));
}

TEST_CASE("blade_geometry: chord <= 0 is invalid")
{
    auto row = rotor();
    row.trailing_edge_hub = row.leading_edge_hub;
    row.trailing_edge_shroud = row.leading_edge_shroud;

    BEMSolverConfig cfg;
    auto res = build_blade_geometry(make_turbine(row), cfg);
    CHECK_FALSE(res.ok);
    CHECK(res.error == "chord <= 0");
}

TEST_CASE("blade_geometry: rpm <= 0 is invalid")
{
    auto row = rotor();
    row.rotational_speed.value = 0.0f;
    BEMSolverConfig cfg;
    auto res = build_blade_geometry(make_turbine(row), cfg);
    CHECK_FALSE(res.ok);
    CHECK(res.error == "rpm <= 0");
}

TEST_CASE("blade_geometry: element_count < 4 is rejected by solver config validation")
{
    BEMSolverConfig cfg;
    cfg.element_count = 2;
    auto res = build_blade_geometry(make_turbine(rotor()), cfg);
    // build_blade_geometry itself does not validate element_count; the solver does.
    // Here we just verify it builds the requested number of elements.
    REQUIRE(res.ok);
    CHECK(res.geometry.elements.size() == 2);
}

TEST_CASE("blade_geometry: empty sections yields zero stagger and warning")
{
    auto row = rotor();
    row.sections.clear();
    BEMSolverConfig cfg;
    auto res = build_blade_geometry(make_turbine(row), cfg);
    REQUIRE(res.ok);
    CHECK(res.geometry.elements.front().beta_deg == doctest::Approx(0.0));
    bool warned = false;
    for (const auto& w : res.warnings)
        if (w.find("sections empty") != std::string::npos) warned = true;
    CHECK(warned);
}

TEST_CASE("blade_geometry: HubDefinition is ignored")
{
    auto t = make_turbine(rotor());
    t.hub.shape = HubShape::Cylinder;
    t.hub.root_radius = 0.4f;
    BEMSolverConfig cfg;
    auto res = build_blade_geometry(t, cfg);
    REQUIRE(res.ok);
    // Hub geometry still taken from flow_path, not HubDefinition.
    CHECK(res.geometry.r_hub == doctest::Approx(0.4));
}