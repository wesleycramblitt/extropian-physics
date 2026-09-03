#pragma once

// ─────────────────────────────────────────────────────
// Mesh generation (implementation_spec §10).
//
// Pipeline: Geometry → Domain Analysis → Mesh Generation →
// Boundary Classification → Mesh Validation → Simulation
// Mesh.  Initial capabilities: uniform refinement, local
// refinement (box region), boundary refinement (halo), and
// resolution controls.  AMR / error-driven refinement are
// roadmap (spec §10 "future").
// ─────────────────────────────────────────────────────

#include <exd/engine/mesh/structured.hpp>

#include <algorithm>
#include <array>
#include <string_view>

namespace exd::engine::mesh {

/// Build a uniform Cartesian grid from origin/extent/resolution.
/// `n_per_axis` is the node count per axis (>= 2).
inline StructuredGrid make_structured_grid(const std::array<double, 3>& origin,
                                           const std::array<double, 3>& extent,
                                           const std::array<int32_t, 3>& n_per_axis)
{
    StructuredGrid g;
    g.origin = origin;
    g.dims = n_per_axis;
    for (int a = 0; a < 3; ++a)
        g.spacing[a] = (n_per_axis[a] > 1) ? extent[a] / (n_per_axis[a] - 1) : extent[a];
    return g;
}

/// A box region selector (half-open node index ranges; negatives = full axis).
struct BoxRegion
{
    std::array<int32_t, 3> lo = {0, 0, 0};
    std::array<int32_t, 3> hi = {-1, -1, -1};  // -1 → dims[a] - 1
};

inline bool in_region(const StructuredGrid& g, const BoxRegion& r,
                      int32_t i, int32_t j, int32_t k)
{
    const std::array<int32_t, 3> idx = {i, j, k};
    for (int a = 0; a < 3; ++a)
    {
        const int32_t lo = r.lo[a];
        const int32_t hi = (r.hi[a] < 0) ? g.dims[a] - 1 : r.hi[a];
        if (idx[a] < lo || idx[a] > hi) return false;
    }
    return true;
}

/// Uniform refinement: subdivide every cell by `factor` per axis and
/// resample node values with trilinear interpolation (values unchanged if
/// no source values provided — pure mesh refinement).
inline StructuredGrid refine(const StructuredGrid& g, int factor)
{
    StructuredGrid out;
    out.origin = g.origin;
    for (int a = 0; a < 3; ++a)
    {
        out.dims[a] = (g.dims[a] - 1) * factor + 1;
        out.spacing[a] = g.spacing[a] / factor;
    }
    return out;
}

/// Resample a scalar node field onto a refined grid (trilinear).
inline std::vector<double> refine_values(const StructuredGrid& src,
                                         const std::vector<double>& src_values,
                                         int factor)
{
    StructuredGrid dst = refine(src, factor);
    std::vector<double> out(dst.node_count(), 0.0);
    for (int32_t k2 = 0; k2 < dst.dims[2]; ++k2)
    for (int32_t j2 = 0; j2 < dst.dims[1]; ++j2)
    for (int32_t i2 = 0; i2 < dst.dims[0]; ++i2)
    {
        const double x = dst.origin[0] + i2 * dst.spacing[0];
        const double y = dst.origin[1] + j2 * dst.spacing[1];
        const double z = dst.origin[2] + k2 * dst.spacing[2];
        const double fx = (x - src.origin[0]) / src.spacing[0];
        const double fy = (y - src.origin[1]) / src.spacing[1];
        const double fz = (z - src.origin[2]) / src.spacing[2];
        const int32_t i0 = std::clamp(static_cast<int32_t>(fx), 0, src.dims[0] - 2);
        const int32_t j0 = std::clamp(static_cast<int32_t>(fy), 0, src.dims[1] - 2);
        const int32_t k0 = std::clamp(static_cast<int32_t>(fz), 0, src.dims[2] - 2);
        const double tx = std::clamp(fx - i0, 0.0, 1.0);
        const double ty = std::clamp(fy - j0, 0.0, 1.0);
        const double tz = std::clamp(fz - k0, 0.0, 1.0);
        const size_t nx = static_cast<size_t>(src.dims[0]);
        const size_t ny = static_cast<size_t>(src.dims[1]);
        auto V = [&](int32_t i, int32_t j, int32_t k) {
            return src_values[static_cast<size_t>(i) + nx * (static_cast<size_t>(j) + ny * k)];
        };
        const double c000 = V(i0, j0, k0), c100 = V(i0 + 1, j0, k0);
        const double c010 = V(i0, j0 + 1, k0), c110 = V(i0 + 1, j0 + 1, k0);
        const double c001 = V(i0, j0, k0 + 1), c101 = V(i0 + 1, j0, k0 + 1);
        const double c011 = V(i0, j0 + 1, k0 + 1), c111 = V(i0 + 1, j0 + 1, k0 + 1);
        const double c00 = c000 + tx * (c100 - c000);
        const double c01 = c001 + tx * (c101 - c001);
        const double c10 = c010 + tx * (c110 - c010);
        const double c11 = c011 + tx * (c111 - c011);
        const double c0 = c00 + ty * (c10 - c00);
        const double c1 = c01 + ty * (c11 - c01);
        const size_t idx = static_cast<size_t>(i2) +
                           static_cast<size_t>(dst.dims[0]) *
                               (static_cast<size_t>(j2) +
                                static_cast<size_t>(dst.dims[1]) * k2);
        out[idx] = c0 + tz * (c1 - c0);
    }
    return out;
}

inline bool validate_mesh(const StructuredGrid& g, ModelStatus& status)
{
    return g.validate(status);
}

} // namespace exd::engine::mesh
