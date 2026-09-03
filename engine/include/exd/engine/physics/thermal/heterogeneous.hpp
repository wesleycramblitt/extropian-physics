#pragma once

// ─────────────────────────────────────────────────────
// Heterogeneous steady conduction (implementation_spec
// §33 materials, §52 "Heat Conduction" class).
//
// THE GENERAL CLASS, not a product preset: steady
// conduction in a heterogeneous solid with per-node
// conductivity, per-node volumetric sources, and
// per-face boundary kinds (fixed temperature or
// insulated):
//
//   −∇·(k(x)·∇T) = q(x)
//
// The user brings the DATA (from CAD bodies, measured
// fields, or region builders): per-node k and q overrides
// are the data-driven contract; box material/source
// regions are the convenience builder.  Whatever they
// import — a chiplet board, a heat sink, a reactor wall,
// a building envelope — becomes the same config and the
// same solver runs.
//
// Face conductivity uses the harmonic mean (exact for
// series conduction).  The operator is affine in the
// fixed-temperature faces; the linear part + CG solve it
// directly (spec §35 matrix-free).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>
#include <vector>

namespace exd::engine::physics::thermal {

/// Per-face boundary kind of the conduction domain.
enum class ConductionFaceKind : uint8_t
{
    FixedValue,   // T = face_value on the face
    Insulated,    // zero normal gradient (mirror ghost)
};

/// A box material region (e.g., one CAD body's bounding box).
struct ConductionRegion
{
    std::array<double, 3> center = {0, 0, 0};
    std::array<double, 3> half_extents = {0.1, 0.1, 0.1};
    double conductivity = 400.0;       // W/(m·K)
};

/// A box volumetric heat-source region (W/m³; from power: q = P/V).
struct SourceRegion
{
    std::array<double, 3> center = {0, 0, 0};
    std::array<double, 3> half_extents = {0.05, 0.05, 0.05};
    double volumetric_source = 0.0;    // W/m³
};

struct HeterogeneousConductionConfig
{
    mesh::StructuredGrid grid;                     // the domain
    double base_conductivity = 20.0;               // W/(m·K), outside regions

    std::vector<ConductionRegion> materials;       // region builder (CAD bodies)
    std::vector<SourceRegion> sources;             // region builder (heat sources)

    // per-face BC in the order {x−, x+, y−, y+, z−, z+}
    std::array<ConductionFaceKind, 6> face_kind = {
        ConductionFaceKind::FixedValue, ConductionFaceKind::FixedValue,
        ConductionFaceKind::FixedValue, ConductionFaceKind::FixedValue,
        ConductionFaceKind::FixedValue, ConductionFaceKind::FixedValue,
    };
    std::array<double, 6> face_value = {300, 300, 300, 300, 300, 300};  // K

    // ── data-driven overrides (the CAD/import contract): when non-empty,
    //    these per-node fields REPLACE the region-built fields. ──
    std::vector<double> conductivity_field;        // per node, W/(m·K)
    std::vector<double> source_field;              // per node, W/m³

    double tolerance = 1e-10;
    int max_iterations = 20000;
};

struct HeterogeneousConductionResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;
    mesh::StructuredScalarGrid temperature;        // K
    double total_power = 0.0;                      // ∫q dV (W)
    double sink_flux = 0.0;                        // net outward flux (W)
    double peak_temperature = 0.0;
    double peak_x = 0.0, peak_y = 0.0;
    double max_residual = 0.0;
};

/// Solve the heterogeneous steady conduction problem (deterministic, no
/// exceptions).  This is THE general "heat conduction" assembly — presets
/// and applications configure it with their own geometry/BCs/parameters.
HeterogeneousConductionResult solve_heterogeneous_conduction(
    const HeterogeneousConductionConfig& config);

} // namespace exd::engine::physics::thermal
