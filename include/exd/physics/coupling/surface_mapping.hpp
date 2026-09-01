#pragma once

#include <exd/physics/coupling/field_channels.hpp>
#include <exd/physics/model_status.hpp>

#include <cstddef>
#include <cstdint>
#include <span>

namespace exd::physics::coupling {

/// Interpolation rule used when mapping values between coupling surfaces.
enum class InterpolationMode : uint8_t
{
    Nearest,     // copy the nearest source node value
    Trilinear,   // trilinear sampling of a structured grid at a physical point
};

/// Outcome of a surface-to-surface transfer.  The running transfer helpers
/// report through the bool return + ModelStatus channel; this value type is
/// kept for callers that want to carry the outcome as a record.
struct TransferResult
{
    bool ok = false;
    exd::physics::ModelStatus status;
    size_t transferred = 0;
};

/// Nearest-neighbor transfer between two point sets.
///
/// Brute-force O(n*m); documented for configuration-size coupling surfaces
/// (tens to hundreds of probe points), not for full surface meshes.  Each
/// target point receives the source_data entry of the nearest source point
/// (Euclidean distance; ties pick the first source point).
///
/// Sizes must match: source_data.size() == source_points.size() and
/// target_data.size() == target_points.size(), otherwise status is set to
/// error and the function returns false.
bool transfer_nearest(std::span<const double> source_data,
                      std::span<const std::array<double, 3>> source_points,
                      std::span<const std::array<double, 3>> target_points,
                      std::span<double> target_data,
                      exd::physics::ModelStatus& status);

/// Trilinear sampling of `source` at every node of `target`.
///
/// Target node physical coordinates are p = origin + spacing * node_index
/// elementwise.  A target node outside the source bounding box (beyond a
/// 1e-12 node-space epsilon) receives quiet_NaN and produces a warning.
/// target_values_out.size() must equal the product of the target dims.
bool transfer_trilinear(const StructuredScalarGrid& source,
                        const StructuredScalarGrid& target,
                        std::span<double> target_values_out,
                        exd::physics::ModelStatus& status);

} // namespace exd::physics::coupling