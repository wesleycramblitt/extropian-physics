#pragma once

// ─────────────────────────────────────────────────────
// Structured mesh (implementation_spec §9–§10).
//
// The initial mesh family is Cartesian/structured: a
// regular node lattice defined by origin, spacing, dims.
// Geometry (origin/spacing), topology (indexing, neighbor
// relations), and metrics (cell volume, face area, face
// normal, centroid) are separate concerns over the same
// lattice.  The API does not assume a fixed topology for
// the rest of the engine — unstructured/AMR families land
// behind the same mesh notions later (spec §9 "future").
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::mesh {

/// Face orientation of a structured hexahedral lattice.
enum class Axis : uint8_t { X = 0, Y = 1, Z = 2 };

/// Boundary patch identifiers of a structured box domain.
enum class BoundaryId : uint8_t
{
    XNeg = 0, XPos = 1, YNeg = 2, YPos = 3, ZNeg = 4, ZPos = 5,
};

inline constexpr const char* to_string(BoundaryId b)
{
    switch (b)
    {
    case BoundaryId::XNeg: return "x-";
    case BoundaryId::XPos: return "x+";
    case BoundaryId::YNeg: return "y-";
    case BoundaryId::YPos: return "y+";
    case BoundaryId::ZNeg: return "z-";
    case BoundaryId::ZPos: return "z+";
    }
    return "?";
}

/// Regular structured grid: node lattice origin/spacing/dims (node-centered
/// values live in the grid's `values`; cell-centered data uses dims-1 cells
/// per axis).
struct StructuredGrid
{
    std::array<double, 3> origin = {0, 0, 0};        // node (0,0,0) (m)
    std::array<double, 3> spacing = {0.1, 0.1, 0.1}; // node spacing per axis (m), > 0
    std::array<int32_t, 3> dims = {0, 0, 0};         // node counts (>= 2 per axis)

    // ── metrics (spec §9) ──
    [[nodiscard]] std::array<double, 3> node_coords(int32_t i, int32_t j, int32_t k) const
    {
        return {origin[0] + i * spacing[0], origin[1] + j * spacing[1], origin[2] + k * spacing[2]};
    }
    [[nodiscard]] int32_t cell_count_x() const { return dims[0] - 1; }
    [[nodiscard]] int32_t cell_count_y() const { return dims[1] - 1; }
    [[nodiscard]] int32_t cell_count_z() const { return dims[2] - 1; }
    [[nodiscard]] size_t node_count() const
    {
        return static_cast<size_t>(dims[0]) * dims[1] * dims[2];
    }
    [[nodiscard]] size_t cell_count() const
    {
        return static_cast<size_t>(cell_count_x()) * cell_count_y() * cell_count_z();
    }
    [[nodiscard]] double cell_volume() const
    {
        return spacing[0] * spacing[1] * spacing[2];
    }
    [[nodiscard]] double face_area(Axis axis) const
    {
        switch (axis)
        {
        case Axis::X: return spacing[1] * spacing[2];
        case Axis::Y: return spacing[0] * spacing[2];
        case Axis::Z: return spacing[0] * spacing[1];
        }
        return 0.0;
    }
    [[nodiscard]] std::array<double, 3> face_normal(Axis axis, bool positive) const
    {
        const double s = positive ? 1.0 : -1.0;
        switch (axis)
        {
        case Axis::X: return {s, 0, 0};
        case Axis::Y: return {0, s, 0};
        case Axis::Z: return {0, 0, s};
        }
        return {0, 0, 0};
    }
    [[nodiscard]] std::array<double, 3> cell_centroid(int32_t i, int32_t j, int32_t k) const
    {
        return node_coords(i + 1, j + 1, k + 1);   // cell (i,j,k) spans nodes i..i+1 etc.
    }
    [[nodiscard]] std::array<double, 3> node_centroid(int32_t i, int32_t j, int32_t k) const
    {
        return node_coords(i, j, k);
    }

    /// Neighbor node index along axis (+1/-1); returns -1 when out of range
    /// (boundary node — the boundary layer resolves BCs).
    [[nodiscard]] int32_t neighbor(int32_t i, int32_t j, int32_t k, Axis axis, bool positive) const
    {
        const int32_t d = positive ? 1 : -1;
        if (axis == Axis::X) return (i + d >= 0 && i + d < dims[0]) ? i + d : -1;
        if (axis == Axis::Y) return (j + d >= 0 && j + d < dims[1]) ? j + d : -1;
        return (k + d >= 0 && k + d < dims[2]) ? k + d : -1;
    }

    [[nodiscard]] bool validate(ModelStatus& status) const;
};

/// Scalar field on a structured grid (flat nx·ny·nz, index i + nx·(j + ny·k)).
/// Inherits the lattice (origin/spacing/dims) so existing solver code keeps
/// the direct `.origin/.spacing/.dims` access pattern.
struct StructuredScalarGrid : StructuredGrid
{
    std::vector<double> values;
    [[nodiscard]] bool validate(ModelStatus& status) const;
};

/// Vector field on a structured grid (3·nx·ny·nz flat).
struct StructuredVectorGrid : StructuredGrid
{
    std::vector<double> values;  // 3 · nx·ny·nz
    [[nodiscard]] bool validate(ModelStatus& status) const;
};

inline bool StructuredGrid::validate(ModelStatus& status) const
{
    for (int a = 0; a < 3; ++a)
    {
        if (dims[a] < 2)
        {
            status.ok = false;
            status.error = "mesh: dims[" + std::to_string(a) + "] must be >= 2";
            return false;
        }
        if (!(spacing[a] > 0.0))
        {
            status.ok = false;
            status.error = "mesh: spacing[" + std::to_string(a) + "] must be > 0";
            return false;
        }
    }
    return true;
}

inline bool StructuredScalarGrid::validate(ModelStatus& status) const
{
    if (!StructuredGrid::validate(status)) return false;
    if (values.size() != node_count())
    {
        status.ok = false;
        status.error = "mesh scalar field: size " + std::to_string(values.size()) +
                       " != node count " + std::to_string(node_count());
        return false;
    }
    return true;
}

inline bool StructuredVectorGrid::validate(ModelStatus& status) const
{
    if (!StructuredGrid::validate(status)) return false;
    if (values.size() != 3 * node_count())
    {
        status.ok = false;
        status.error = "mesh vector field: size " + std::to_string(values.size()) +
                       " != 3 × node count " + std::to_string(3 * node_count());
        return false;
    }
    return true;
}

} // namespace exd::engine::mesh
