#pragma once

// ─────────────────────────────────────────────────────
// Coupling contracts (implementation_spec §17–§21, §25).
//
// Every coupling defines a contract: source, destination,
// quantity, units, field type/rank, spatial association,
// temporal behavior, mapping, interpolation, conservation
// requirement, sign convention, execution frequency,
// coupling strength.  Validation runs before execution —
// invalid combinations (wrong units, rank mismatch,
// unsupported mapping) are rejected with diagnostics.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/core/units.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>
#include <string>
#include <vector>

namespace exd::engine::coupling {

/// Temporal behavior of a coupling (spec §19).
enum class TemporalBehavior : uint8_t
{
    Explicit,    // evaluated once per exchange, current data
    Lagged,      // uses previous-timestep source data
    Staggered,   // half-step offset (operator splitting)
    Iterative,   // sub-iterated to convergence within the step
    Implicit,    // part of an implicit system solve
};

/// Spatial mapping family (spec §18).
enum class MappingKind : uint8_t
{
    SameMesh,        // identical association on the same mesh
    NodeToCell,
    CellToNode,
    FaceToCell,
    SurfaceToBody,   // surface integral onto a body
    ParticleToField,
    FieldToParticle,
    MeshToMesh,      // different meshes
};

/// Interpolation/projection method (spec §18).
enum class InterpolationKind : uint8_t
{
    Nearest,
    Linear,          // (bi/trilinear)
    HigherOrder,
    Conservative,    // conservative projection (flux/mass preserving)
    L2Projection,
    WeightedAverage,
    SurfaceIntegral,
    VolumeIntegral,
};

/// Conservation requirement (spec §17–§18, §21).
enum class ConservationRequirement : uint8_t
{
    None,
    Recommended,
    Required,        // arbitrary interpolation must not be silently allowed
};

constexpr const char* to_string(TemporalBehavior b)
{
    switch (b)
    {
    case TemporalBehavior::Explicit: return "explicit";
    case TemporalBehavior::Lagged: return "lagged";
    case TemporalBehavior::Staggered: return "staggered";
    case TemporalBehavior::Iterative: return "iterative";
    case TemporalBehavior::Implicit: return "implicit";
    }
    return "?";
}

constexpr const char* to_string(MappingKind m)
{
    switch (m)
    {
    case MappingKind::SameMesh: return "same_mesh";
    case MappingKind::NodeToCell: return "node_to_cell";
    case MappingKind::CellToNode: return "cell_to_node";
    case MappingKind::FaceToCell: return "face_to_cell";
    case MappingKind::SurfaceToBody: return "surface_to_body";
    case MappingKind::ParticleToField: return "particle_to_field";
    case MappingKind::FieldToParticle: return "field_to_particle";
    case MappingKind::MeshToMesh: return "mesh_to_mesh";
    }
    return "?";
}

/// A full coupling contract (§17).  `validate` checks internal consistency;
/// `compatible` checks the contract against the actual source/destination
/// field metadata.
struct CouplingContract
{
    std::string id;
    std::string source_domain;       // e.g. "fluid"
    std::string source_quantity;     // e.g. "wall_shear"
    std::string destination_domain;  // e.g. "rigid_body"
    std::string destination_quantity;

    // ── quantity contract (§21: quantity type, dimensional units, rank) ──
    core::Unit units;                            // of the transferred quantity
    core::FieldRank rank = core::FieldRank::Scalar;
    core::FieldLocation source_association = core::FieldLocation::Node;
    core::FieldLocation destination_association = core::FieldLocation::Node;

    // ── spatial/temporal contract (§18–§20) ──
    MappingKind mapping = MappingKind::SameMesh;
    InterpolationKind interpolation = InterpolationKind::Linear;
    std::vector<std::array<double, 3>> probe_points;   // where the exchange samples/writes
    ConservationRequirement conservation = ConservationRequirement::None;
    TemporalBehavior temporal = TemporalBehavior::Explicit;
    double execution_interval = 0.0;   // seconds between exchanges (0 = every step)
    double coupling_strength = 1.0;    // relaxation factor in (0, 1]
    bool sign_convention_reversed = false;

    /// Internal consistency of the contract itself.
    [[nodiscard]] bool validate(core::ModelStatus& status) const;
};

/// Spec §21/§25: reject invalid combinations before execution.
/// Verifies the contract fields against the actual source/destination
/// field metadata delivered by the domains.
[[nodiscard]] inline bool contract_compatible(const CouplingContract& c,
                                              const core::FieldMetadata& source_field,
                                              const core::FieldMetadata& destination_field,
                                              core::ModelStatus& status)
{
    if (!c.validate(status)) return false;

    if (c.rank != source_field.rank || c.rank != destination_field.rank)
    {
        status.ok = false;
        status.error = "coupling '" + c.id + "': rank mismatch — source " +
                       std::to_string(static_cast<int>(source_field.rank)) + ", destination " +
                       std::to_string(static_cast<int>(destination_field.rank)) +
                       ", contract " + std::to_string(static_cast<int>(c.rank));
        return false;
    }
    if (!core::units_compatible(c.units, source_field.units))
    {
        status.ok = false;
        status.error = "coupling '" + c.id + "': source field '" + source_field.name +
                       "' units mismatch (contract " + core::to_string(c.units) +
                       ", field " + core::to_string(source_field.units) + ")";
        return false;
    }
    if (!core::units_compatible(c.units, destination_field.units))
    {
        status.ok = false;
        status.error = "coupling '" + c.id + "': destination field '" + destination_field.name +
                       "' units mismatch (contract " + core::to_string(c.units) +
                       ", field " + core::to_string(destination_field.units) + ")";
        return false;
    }
    if (c.conservation == ConservationRequirement::Required &&
        c.interpolation != InterpolationKind::Conservative &&
        c.interpolation != InterpolationKind::SurfaceIntegral &&
        c.interpolation != InterpolationKind::VolumeIntegral &&
        c.interpolation != InterpolationKind::WeightedAverage)
    {
        status.ok = false;
        status.error = "coupling '" + c.id +
                       "': conservation-required transfer with non-conservative interpolation '"
                       + std::to_string(static_cast<int>(c.interpolation)) +
                       "' — specify a conservative projection (spec §18)";
        status.warnings.push_back("suggested mappings: conservative projection, "
                                  "surface/volume integration");
        return false;
    }
    if (c.coupling_strength <= 0.0 || c.coupling_strength > 1.0)
    {
        status.ok = false;
        status.error = "coupling '" + c.id + "': coupling_strength must be in (0, 1]";
        return false;
    }
    return true;
}

inline bool CouplingContract::validate(core::ModelStatus& status) const
{
    if (id.empty() || source_domain.empty() || source_quantity.empty() ||
        destination_domain.empty() || destination_quantity.empty())
    {
        status.ok = false;
        status.error = "coupling contract: id/source/destination must not be empty";
        return false;
    }
    if (source_domain == destination_domain && source_quantity == destination_quantity)
    {
        status.ok = false;
        status.error = "coupling '" + id + "': source == destination (no-op coupling)";
        return false;
    }
    if (mapping == MappingKind::SameMesh &&
        source_association != destination_association)
    {
        status.ok = false;
        status.error = "coupling '" + id +
                       "': same-mesh mapping requires identical associations";
        return false;
    }
    return true;
}

} // namespace exd::engine::coupling
