// field_test.cpp — spec §5/§6/§7: Field metadata, State, EntitySet.
#include <exd/engine/core/entity_set.hpp>
#include <exd/engine/core/field.hpp>
#include <exd/engine/core/material.hpp>
#include <exd/engine/core/state.hpp>
#include <doctest/doctest.h>

using namespace exd::engine::core;

TEST_CASE("Field: metadata + contiguous storage")
{
    ModelStatus st;
    FieldMetadata meta;
    meta.name = "temperature";
    meta.rank = FieldRank::Scalar;
    meta.components = 1;
    meta.units = units::kelvin;
    meta.location = FieldLocation::Node;
    meta.precision = FieldPrecision::F64;
    meta.domain = "thermal";
    meta.mesh = "rod";
    CHECK(meta.is_valid(st));

    Field f(meta, 21);
    CHECK(f.size() == 21);
    CHECK(f.entity_count() == 21);
    f.assign(300.0);
    CHECK(f.at(7) == 300.0);
    f.set(7, 350.0);
    CHECK(f.version() > 0);
    CHECK(f.validate(st));

    // vector field: components == 3
    FieldMetadata vmeta;
    vmeta.name = "velocity";
    vmeta.rank = FieldRank::Vector;
    vmeta.components = 3;
    Field vf(vmeta, 9);
    CHECK(vf.entity_count() == 3);
    CHECK(vf.entity(1).size() == 3);
}

TEST_CASE("Field: invalid metadata rejected")
{
    ModelStatus st;
    FieldMetadata bad;
    bad.name = "";                       // empty name
    CHECK(!bad.is_valid(st));
    CHECK(!st.ok);

    st = ModelStatus{};
    FieldMetadata rankbad;
    rankbad.name = "x";
    rankbad.rank = FieldRank::Vector;    // vector requires 3 components
    rankbad.components = 1;
    CHECK(!rankbad.is_valid(st));
}

TEST_CASE("EntitySet: contiguous indexed collection")
{
    ModelStatus st;
    EntitySet cells("cells", EntityKind::Cells, 1000);
    cells.add_attribute("volume", std::vector<double>(1000, 1.0e-6), st);
    REQUIRE(st.ok);
    CHECK(cells.count() == 1000);
    CHECK(cells.attribute("volume")->size() == 1000);
    CHECK(cells.validate(st));

    // wrong-size attribute rejected
    EntitySet bad("bad", EntityKind::Nodes, 10);
    bad.add_attribute("x", std::vector<double>(9), st);
    CHECK(!st.ok);
}

TEST_CASE("State: fields, entity sets, versioning, residency")
{
    ModelStatus st;
    State s("test_state");
    FieldMetadata meta;
    meta.name = "p";
    meta.rank = FieldRank::Scalar;
    meta.components = 1;
    meta.units = units::pascal;

    s.add_field(meta, 16, st);
    REQUIRE(st.ok);
    // duplicate rejected
    s.add_field(meta, 16, st);
    CHECK(!st.ok);
    st = ModelStatus{};

    s.add_entity_set(EntitySet("nodes", EntityKind::Nodes, 16), st);
    REQUIRE(st.ok);

    CHECK(s.field("p") != nullptr);
    CHECK(s.field("nope") == nullptr);
    CHECK(s.validate(st));
    s.set_space(MemorySpace::Cpu);
    CHECK(s.space() == MemorySpace::Cpu);
    s.sync_from_cpu();
    s.sync_to_cpu();   // CPU no-ops; GPU path is backend-owned (Phase 11)
    CHECK(s.version() > 0);
}

TEST_CASE("Materials: constitutive properties independent of discretization (§33)")
{
    using exd::engine::core::MaterialDatabase;
    using exd::engine::core::MaterialProperties;
    using exd::engine::core::MaterialType;
    using exd::engine::core::ModelStatus;
    ModelStatus st;
    MaterialDatabase db;
    MaterialProperties air;
    air.name = "air";
    air.type = MaterialType::Fluid;
    air.density = 1.225;
    air.dynamic_viscosity = 1.81e-5;
    air.thermal_conductivity = 0.0257;
    air.specific_heat = 1005.0;
    REQUIRE(db.add(air, st));
    CHECK(db.size() == 1);

    const auto* found = db.find("air");
    REQUIRE(found != nullptr);
    CHECK(found->kinematic_viscosity() == doctest::Approx(1.81e-5 / 1.225).epsilon(1e-12));
    CHECK(found->thermal_diffusivity() == doctest::Approx(0.0257 / (1.225 * 1005.0)).epsilon(1e-12));

    // invalid materials rejected; duplicates rejected
    MaterialProperties bad;
    bad.name = "bad";
    bad.density = 0.0;
    CHECK(!db.add(bad, st));
    CHECK(!st.ok);
    st = ModelStatus{};
    CHECK(!db.add(air, st));   // duplicate
    CHECK(!st.ok);
}
