// contracts_test.cpp — spec §17/§21/§22/§25: coupling contracts, units,
// conservation requirements, compatibility rules.
#include <exd/engine/coupling/contracts.hpp>
#include <exd/engine/coupling/rules.hpp>
#include <doctest/doctest.h>

using namespace exd::engine;
using namespace exd::engine::coupling;
using namespace exd::engine::core;

TEST_CASE("CouplingContract: fluid wall shear → rigid body torque (§17 example)")
{
    ModelStatus st;
    CouplingContract c;
    c.id = "shear_to_torque";
    c.source_domain = "fluid";
    c.source_quantity = "wall_shear";
    c.destination_domain = "rigid_body";
    c.destination_quantity = "torque";
    c.units = units::torque;            // Pa·m → N·m after surface integral
    c.rank = FieldRank::Vector;
    c.source_association = FieldLocation::Face;
    c.destination_association = FieldLocation::RigidBody;
    c.mapping = MappingKind::SurfaceToBody;
    c.interpolation = InterpolationKind::SurfaceIntegral;
    c.conservation = ConservationRequirement::Required;
    c.temporal = TemporalBehavior::Iterative;
    c.execution_interval = 0.0;         // every timestep
    CHECK(c.validate(st));
    CHECK(st.ok);
}

TEST_CASE("CouplingContract: temperature → force rejected without transformation")
{
    ModelStatus st;
    CouplingContract c;
    c.id = "bad";
    c.source_domain = "thermal";
    c.source_quantity = "temperature";
    c.destination_domain = "rigid_body";
    c.destination_quantity = "force";
    c.units = units::kelvin;            // T → F is dimensionally invalid (§25)
    c.rank = FieldRank::Scalar;
    CHECK(c.validate(st));
    // validation against the actual field metadata catches the mismatch
    FieldMetadata src;
    src.name = "temperature";
    src.rank = FieldRank::Scalar;
    src.units = units::kelvin;
    FieldMetadata dst;
    dst.name = "force";
    dst.rank = FieldRank::Vector;
    dst.units = units::force;
    CHECK(!contract_compatible(c, src, dst, st));   // rank + units both fail
    CHECK(!st.ok);
}

TEST_CASE("CouplingContract: conservation-required + non-conservative interpolation rejected")
{
    ModelStatus st;
    CouplingContract c;
    c.id = "conservation";
    c.source_domain = "fluid";
    c.source_quantity = "flux";
    c.destination_domain = "thermal";
    c.destination_quantity = "flux";
    c.units = units::watt;
    c.rank = FieldRank::Scalar;
    c.source_association = FieldLocation::Face;
    c.destination_association = FieldLocation::Face;
    c.mapping = MappingKind::MeshToMesh;
    c.interpolation = InterpolationKind::Linear;        // not conservative
    c.conservation = ConservationRequirement::Required;
    FieldMetadata src;
    src.name = "s";
    src.rank = FieldRank::Scalar;
    src.units = units::watt;
    FieldMetadata dst;
    dst.name = "d";
    dst.rank = FieldRank::Scalar;
    dst.units = units::watt;
    const bool compatible = contract_compatible(c, src, dst, st);
    CHECK(!compatible);
    CHECK(!st.ok);
    CHECK(!st.warnings.empty());       // diagnostic suggests conservative mapping
}

TEST_CASE("Compatibility rules: registry enforces module discretization support")
{
    ModelStatus st;
    RuleRegistry reg = default_rule_registry();
    ValidationContext ctx;
    ctx.module_name = "thermal";
    ctx.discretization = "FVM";        // not in supported set
    ctx.mesh_family = "structured";
    DiscretizationRule dr({"FDM"});
    CHECK(!dr.check(ctx, st));
    CHECK(!st.ok);
    st = ModelStatus{};

    ctx.discretization = "FDM";
    REQUIRE(reg.validate(ctx, st));
    CHECK(st.ok);
    // unstructured mesh + FDM rejected (§9)
    st = ModelStatus{};
    ctx.mesh_family = "unstructured";
    CHECK(!reg.validate(ctx, st));
}
