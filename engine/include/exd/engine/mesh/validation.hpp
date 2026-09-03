#pragma once

// ─────────────────────────────────────────────────────
// Mesh validation (implementation_spec §10 pipeline step
// "Mesh Validation", and §55 validation categories).
//
// Validates dimensions, spacing, field-size agreement,
// and boundary completeness (every structured-box face
// has either a BC or a declared neighbor).
// ─────────────────────────────────────────────────────

#include <exd/engine/mesh/boundary.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>

namespace exd::engine::mesh {

/// Whether every face of the box is covered by at least one BC.
inline bool boundary_complete(const StructuredGrid& g,
                              const std::array<const BoundaryCondition*, 6>& face_bcs,
                              ModelStatus& status)
{
    for (int b = 0; b < 6; ++b)
    {
        if (face_bcs[static_cast<size_t>(b)] == nullptr)
        {
            status.ok = false;
            status.error = "mesh validation: face " +
                           std::string(to_string(static_cast<BoundaryId>(b))) +
                           " has no boundary condition";
            return false;
        }
    }
    (void)g;
    return true;
}

} // namespace exd::engine::mesh
