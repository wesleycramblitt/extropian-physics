// turbine_builder_test.cpp
// Parametric TurbineDefinition builder: geometry mapping + validation.

#include <exd/physics/turbine/turbine_builder.hpp>

#include <doctest/doctest.h>

#include <cmath>

using namespace exd::physics::turbine;
using exd::physics::ModelStatus;

TEST_CASE("builder: maps parameters onto a single rotor row")
{
    TurbineBuilderConfig cfg;
    cfg.hub_radius = 0.2;
    cfg.tip_radius = 1.0;
    cfg.chord = 0.3;
    cfg.twist_hub_deg = 1.0;
    cfg.twist_tip_deg = -3.0;
    cfg.rpm = 45.0;
    cfg.blade_count = 4;
    cfg.section_count = 3;
    cfg.leading_edge_z = 0.4;
    cfg.duct_length = 1.5;

    ModelStatus status;
    auto t = make_turbine_definition(cfg, status);
    REQUIRE(status.ok);
    REQUIRE(t.blade_rows.size() == 1U);

    const auto& row = t.blade_rows[0];
    CHECK(row.type == exd::geometry::BladeRowType::Rotor);
    CHECK(row.blade_count.value == doctest::Approx(4.0));
    CHECK(row.rotational_speed.value == doctest::Approx(45.0));
    CHECK(row.leading_edge_hub[0] == doctest::Approx(0.4f).epsilon(1e-6));
    CHECK(row.leading_edge_hub[1] == doctest::Approx(0.2f).epsilon(1e-6));
    CHECK(row.leading_edge_shroud[1] == doctest::Approx(1.0f).epsilon(1e-6));
    CHECK(row.trailing_edge_shroud[0] == doctest::Approx(0.7f).epsilon(1e-6)); // LE + chord

    REQUIRE(row.sections.size() == 3U);
    CHECK(row.sections[0].span == doctest::Approx(0.0f).epsilon(1e-6));
    CHECK(row.sections[2].span == doctest::Approx(1.0f).epsilon(1e-6));
    CHECK(row.sections[0].chord.value == doctest::Approx(0.3f).epsilon(1e-6));
    // twist: linear hub → tip
    CHECK(row.sections[0].stagger.value == doctest::Approx(1.0f).epsilon(1e-6));
    CHECK(row.sections[1].stagger.value == doctest::Approx(-1.0f).epsilon(1e-6));
    CHECK(row.sections[2].stagger.value == doctest::Approx(-3.0f).epsilon(1e-6));

    CHECK(t.flow_path.hub_points.size() == 2U);
    CHECK(t.flow_path.hub_points[1][1] == doctest::Approx(0.2f).epsilon(1e-6));
    CHECK(t.flow_path.shroud_points[1][1] == doctest::Approx(1.0f).epsilon(1e-6));
}

TEST_CASE("builder: invalid inputs fail cleanly")
{
    ModelStatus status;

    auto bad = [&](TurbineBuilderConfig cfg, const char* needle)
    {
        ModelStatus st;
        auto t = make_turbine_definition(cfg, st);
        CHECK_FALSE(st.ok);
        CHECK(st.error.find(needle) != std::string::npos);
        CHECK(t.blade_rows.empty());
    };

    auto c1 = TurbineBuilderConfig{};
    c1.hub_radius = 0.0;
    bad(c1, "hub_radius");

    auto c2 = TurbineBuilderConfig{};
    c2.tip_radius = 0.1; // < hub 0.4
    bad(c2, "tip_radius");

    auto c3 = TurbineBuilderConfig{};
    c3.chord = 0.0;
    bad(c3, "chord");

    auto c4 = TurbineBuilderConfig{};
    c4.rpm = 0.0;
    bad(c4, "rpm");

    auto c5 = TurbineBuilderConfig{};
    c5.blade_count = 0;
    bad(c5, "blade_count");

    auto c6 = TurbineBuilderConfig{};
    c6.duct_length = 0.2; // < LE_z 0.3 + chord 0.5 + 0.1
    bad(c6, "duct_length");

    auto c7 = TurbineBuilderConfig{};
    c7.shroud_radius = 0.1; // < tip_radius
    bad(c7, "shroud_radius");
}

TEST_CASE("builder: open rotor defaults to tip radius shroud")
{
    ModelStatus status;
    auto t = make_turbine_definition(TurbineBuilderConfig{}, status);
    REQUIRE(status.ok);
    CHECK(t.flow_path.shroud_points[0][1] == doctest::Approx(2.0f).epsilon(1e-6));
}
