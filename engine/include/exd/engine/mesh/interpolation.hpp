// ──────────────────────────────────────────────────────────────────────────
// Node-interpolation primitives shared by every sampler and channel adapter.
//
// The trilinear combination of the eight cell-corner values was previously
// copy-pasted (verbatim) into field_sampler.cpp, field_channels.cpp,
// surface_mapping.cpp and the thermal channel adapter, with only the
// bracket/clamp policy differing around it.  This header owns the arithmetic
// once; the callers keep their index/stride conventions (nodes are indexed
// flat, i + nx*(j + ny*k); `stride` scales the +i/+j/+k neighbor offsets —
// 1 for scalar grids, 3 for vector components).
// ──────────────────────────────────────────────────────────────────────────
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace exd::engine::mesh::interp {

/// Trilinear value at local cell fractions (fx, fy, fz) inside the cell whose
/// lower node sits at flat offset `base`.  `nx`/`ny` are the grid extents of
/// the two slow directions; `stride` is the per-component node spacing.
inline double trilinear(const std::vector<double>& values,
                        std::size_t base, std::size_t stride,
                        int nx, int ny,
                        double fx, double fy, double fz)
{
    const std::size_t sx = stride;                                            // +i neighbor
    const std::size_t sy = stride * static_cast<std::size_t>(nx);             // +j neighbor
    const std::size_t sz = stride * static_cast<std::size_t>(nx) *
                           static_cast<std::size_t>(ny);                      // +k neighbor

    const std::size_t i000 = base;
    const std::size_t i100 = base + sx;
    const std::size_t i010 = base + sy;
    const std::size_t i110 = base + sx + sy;
    const std::size_t i001 = base + sz;
    const std::size_t i101 = base + sx + sz;
    const std::size_t i011 = base + sy + sz;
    const std::size_t i111 = base + sx + sy + sz;

    const double c00 = values[i000] * (1.0 - fx) + values[i100] * fx;
    const double c10 = values[i010] * (1.0 - fx) + values[i110] * fx;
    const double c01 = values[i001] * (1.0 - fx) + values[i101] * fx;
    const double c11 = values[i011] * (1.0 - fx) + values[i111] * fx;

    const double c0 = c00 * (1.0 - fy) + c10 * fy;
    const double c1 = c01 * (1.0 - fy) + c11 * fy;

    return c0 * (1.0 - fz) + c1 * fz;
}

} // namespace exd::engine::mesh::interp
