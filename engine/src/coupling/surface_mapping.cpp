// Surface mapping: nearest-neighbor point transfer plus trilinear sampling of
// a structured scalar grid.  These are the raw mapping primitives used by the
// coupling manager to move field data between non-matching coupling surfaces
// (probe point sets and regular grids).

#include <exd/engine/coupling/surface_mapping.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace exd::engine::coupling
{

namespace
{

// Node-space tolerance for exact-boundary queries: a point that misses the
// bounding box by less than this (in normalized coordinates) is clamped onto
// the boundary instead of being reported out-of-bounds.
constexpr double kBoundaryEpsilon = 1e-12;

bool valid_scalar_grid(const StructuredScalarGrid& grid)
{
    for (int32_t d : grid.dims)
    {
        if (d < 2)
            return false;
    }
    for (double s : grid.spacing)
    {
        if (!(s > 0.0)) // also catches NaN
            return false;
    }
    const int64_t nx = grid.dims[0];
    const int64_t ny = grid.dims[1];
    const int64_t nz = grid.dims[2];
    const std::size_t node_count = static_cast<std::size_t>(nx * ny * nz);
    return grid.values.size() == node_count;
}

// Trilinear sampling of a validated scalar grid at a physical point.
// Returns false when the point lies outside the source bounding box (with the
// node-space epsilon above).  In-bounds points are clamped onto the last node
// so an exact upper-boundary query stays a pure node value.
bool sample_trilinear(const StructuredScalarGrid& grid,
                      const std::array<double, 3>& p,
                      double& value_out)
{
    const int nx = grid.dims[0];
    const int ny = grid.dims[1];
    const int nz = grid.dims[2];

    double n[3];
    for (int a = 0; a < 3; ++a)
    {
        n[a] = (p[a] - grid.origin[a]) / grid.spacing[a];
        const double upper = static_cast<double>(grid.dims[a] - 1);
        if (!std::isfinite(n[a]) || n[a] < -kBoundaryEpsilon ||
            n[a] > upper + kBoundaryEpsilon)
        {
            return false;
        }
        n[a] = std::clamp(n[a], 0.0, upper);
    }

    int i = static_cast<int>(n[0]);
    int j = static_cast<int>(n[1]);
    int k = static_cast<int>(n[2]);
    i = std::min(i, nx - 2);
    j = std::min(j, ny - 2);
    k = std::min(k, nz - 2);

    const double fx = n[0] - static_cast<double>(i);
    const double fy = n[1] - static_cast<double>(j);
    const double fz = n[2] - static_cast<double>(k);

    const std::size_t sx = 1u;
    const std::size_t sy = static_cast<std::size_t>(nx);
    const std::size_t sz = static_cast<std::size_t>(nx) *
                           static_cast<std::size_t>(ny);
    const std::size_t base = static_cast<std::size_t>(i) +
                             sy * static_cast<std::size_t>(j) +
                             sz * static_cast<std::size_t>(k);

    const std::size_t i000 = base;
    const std::size_t i100 = base + sx;
    const std::size_t i010 = base + sy;
    const std::size_t i110 = base + sx + sy;
    const std::size_t i001 = base + sz;
    const std::size_t i101 = base + sx + sz;
    const std::size_t i011 = base + sy + sz;
    const std::size_t i111 = base + sx + sy + sz;

    const double c00 = grid.values[i000] * (1.0 - fx) + grid.values[i100] * fx;
    const double c10 = grid.values[i010] * (1.0 - fx) + grid.values[i110] * fx;
    const double c01 = grid.values[i001] * (1.0 - fx) + grid.values[i101] * fx;
    const double c11 = grid.values[i011] * (1.0 - fx) + grid.values[i111] * fx;

    const double c0 = c00 * (1.0 - fy) + c10 * fy;
    const double c1 = c01 * (1.0 - fy) + c11 * fy;

    value_out = c0 * (1.0 - fz) + c1 * fz;
    return true;
}

} // anonymous namespace

bool transfer_nearest(std::span<const double> source_data,
                      std::span<const std::array<double, 3>> source_points,
                      std::span<const std::array<double, 3>> target_points,
                      std::span<double> target_data,
                      exd::engine::core::ModelStatus& status)
{
    status = exd::engine::core::ModelStatus{};

    if (source_data.size() != source_points.size())
    {
        status.ok = false;
        status.error = "transfer_nearest: source_data.size() (" +
                       std::to_string(source_data.size()) +
                       ") != source_points.size() (" +
                       std::to_string(source_points.size()) + ")";
        return false;
    }
    if (target_data.size() != target_points.size())
    {
        status.ok = false;
        status.error = "transfer_nearest: target_data.size() (" +
                       std::to_string(target_data.size()) +
                       ") != target_points.size() (" +
                       std::to_string(target_points.size()) + ")";
        return false;
    }
    if (source_points.empty())
    {
        status.ok = false;
        status.error = "transfer_nearest: source point set is empty";
        return false;
    }

    for (std::size_t i = 0; i < target_points.size(); ++i)
    {
        double best_dist2 = std::numeric_limits<double>::infinity();
        std::size_t best = 0;
        for (std::size_t s = 0; s < source_points.size(); ++s)
        {
            double d2 = 0.0;
            for (int a = 0; a < 3; ++a)
            {
                const double d = target_points[i][a] - source_points[s][a];
                d2 += d * d;
            }
            if (d2 < best_dist2)
            {
                best_dist2 = d2;
                best = s;
            }
        }
        target_data[i] = source_data[best];
    }

    return true;
}

bool transfer_trilinear(const StructuredScalarGrid& source,
                        const StructuredScalarGrid& target,
                        std::span<double> target_values_out,
                        exd::engine::core::ModelStatus& status)
{
    status = exd::engine::core::ModelStatus{};

    if (!valid_scalar_grid(source))
    {
        status.ok = false;
        status.error = "transfer_trilinear: source grid invalid "
                       "(need dims >= 2 per axis, spacing > 0, matching values size)";
        return false;
    }

    for (int32_t d : target.dims)
    {
        if (d < 1)
        {
            status.ok = false;
            status.error = "transfer_trilinear: target grid dims must be >= 1";
            return false;
        }
    }

    const int64_t t_nx = target.dims[0];
    const int64_t t_ny = target.dims[1];
    const int64_t t_nz = target.dims[2];
    const std::size_t target_nodes = static_cast<std::size_t>(t_nx * t_ny * t_nz);
    if (target_values_out.size() != target_nodes)
    {
        status.ok = false;
        status.error = "transfer_trilinear: target_values_out.size() (" +
                       std::to_string(target_values_out.size()) +
                       ") != product of target dims (" +
                       std::to_string(target_nodes) + ")";
        return false;
    }

    const std::size_t tnx = static_cast<std::size_t>(t_nx);
    const std::size_t tny = static_cast<std::size_t>(t_ny);

    std::size_t out_of_bounds = 0;
    for (std::size_t idx = 0; idx < target_nodes; ++idx)
    {
        const std::size_t k = idx / (tnx * tny);
        const std::size_t rem = idx % (tnx * tny);
        const std::size_t j = rem / tnx;
        const std::size_t i = rem % tnx;

        std::array<double, 3> p;
        p[0] = target.origin[0] + target.spacing[0] * static_cast<double>(i);
        p[1] = target.origin[1] + target.spacing[1] * static_cast<double>(j);
        p[2] = target.origin[2] + target.spacing[2] * static_cast<double>(k);

        double sample = 0.0;
        if (sample_trilinear(source, p, sample))
        {
            target_values_out[idx] = sample;
        }
        else
        {
            target_values_out[idx] = std::numeric_limits<double>::quiet_NaN();
            ++out_of_bounds;
        }
    }

    if (out_of_bounds > 0)
    {
        status.warnings.push_back("transfer_trilinear: " + std::to_string(out_of_bounds) +
                                  " target node(s) outside the source grid; set to NaN");
    }

    return true;
}

} // namespace exd::engine::coupling