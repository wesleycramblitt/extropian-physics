#pragma once

// ─────────────────────────────────────────────────────
// Per-cell momentum sources for the fdm3 solver (W16).
//
// Two utilities that fill the N = nx·ny·nz flat force
// arrays consumed by FDM3Solver::set_body_force (units:
// ACCELERATION, m/s² per cell — the solver applies
// u += dt·f):
//
//   * IMMERSED MOVING-SOLID PENALTY ("frozen cells"):
//     inside a solid region the flow velocity is forced
//     toward the region's velocity:
//         f = −K·(u − u_solid),   K in 1/s
//     so a moving solid (a "vein", an impeller blade, a
//     piston, a fish…) is represented WITHOUT remeshing:
//     the mask is rebuilt each step from the moving
//     geometry, the grid stays fixed.
//
//   * BOUSSINESQ BUOYANCY:
//         f = −β·(T − T_ref)·g
//     the classic thermal→flow coupling (natural
//     convection) from a per-cell temperature field,
//     typically sampled from the thermal module's channel.
//
// Both are pure per-step sources: the caller composes
// them with any other body forces and sets the sum via
// set_body_force (the FSI-lite pattern).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/coupling/field_channels.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_config.hpp>

#include <array>
#include <span>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::fdm3 {

// ── cell geometry helpers ────────────────────────────

/// Number of fluid cells: N = nx·ny·nz (the force-array size).
inline size_t fdm3_cell_count(const FDM3Config& config)
{
    return static_cast<size_t>(config.nx) * config.ny * config.nz;
}

/// Flat 0-based cell index (matches FDM3Solver::set_body_force's layout).
inline size_t fdm3_cell_index(const FDM3Config& config, int i, int j, int k)
{
    return static_cast<size_t>(i) + static_cast<size_t>(config.nx) *
               (static_cast<size_t>(j) + static_cast<size_t>(config.ny) * k);
}

/// Cell-center coordinates (m).  The fdm3 domain is [0, lx]×[0, ly]×[0, lz]
/// with dx = lx/nx cells; cell (i,j,k) sits at ((i+½)·dx, …).
inline std::array<double, 3> fdm3_cell_center(const FDM3Config& config,
                                              int i, int j, int k)
{
    return {
        (static_cast<double>(i) + 0.5) * (config.lx / config.nx),
        (static_cast<double>(j) + 0.5) * (config.ly / config.ny),
        (static_cast<double>(k) + 0.5) * (config.lz / config.nz),
    };
}

// ── immersed moving-solid penalty ────────────────────

enum class ImmersedShape : uint8_t
{
    Sphere,
    Box,
};

/// One moving solid region inside the fdm3 domain.  The caller advances
/// `center`/`velocity` each step (from a rigid-body state, a prescribed
/// motion, or hand edits) — the mask follows the geometry without remeshing.
struct ImmersedSolid
{
    std::string name;
    ImmersedShape shape = ImmersedShape::Sphere;
    std::array<double, 3> center = {0, 0, 0};      // current centroid (m)
    std::array<double, 3> velocity = {0, 0, 0};    // region velocity (m/s)
    double radius = 0.1;                           // Sphere
    std::array<double, 3> half_extents = {0.1, 0.1, 0.1};  // Box
    /// Smoothing width of the solid-fraction transition (m); 0 = use the
    /// max cell size (one-cell smearing).  The smeared transition is what
    /// keeps the momentum source resolvable by the pressure projection.
    double blend_width = 0.0;
    double penalty = 100.0;                        // K (1/s): the relaxation
                                                   // rate toward u_solid.
                                                   // K·dt >= 2 is unstable for
                                                   // the explicit update —
                                                   // keep K·dt < 2 (e.g. K=100,
                                                   // dt=0.01 → one-step relax).
};

/// Mask of cells inside ANY solid region (diagnostics, rendering, editing).
bool immersed_solid_mask(const FDM3Config& config,
                         const std::vector<ImmersedSolid>& solids,
                         std::vector<bool>& mask);

/// Add the immersed-solid penalty to the force arrays:
///   f += −K·(u_cell − u_solid)  inside the regions
/// `u/v/w` are the CURRENT cell-centered velocities (nx·ny·nz); `fx/fy/fz`
/// are accumulated in place (start from 0 or from other forces).
///
/// NOTE: K·dt must stay below ~2 for the explicit update's stability; the
/// achievable freeze strength is therefore bounded by the neighbor-coupling
/// scale.  For hard frozen cells prefer `apply_kinematic_freeze` below.
bool add_immersed_solid_forces(const FDM3Config& config,
                               const std::vector<ImmersedSolid>& solids,
                               std::span<const double> u,
                               std::span<const double> v,
                               std::span<const double> w,
                               std::span<double> fx,
                               std::span<double> fy,
                               std::span<double> fz,
                               exd::engine::core::ModelStatus& status);

/// Kinematic freeze: after a solver step, blend the occupied cells'
/// velocity toward u_solid:  u ← (1−a)·u + a·u_solid.  Unconditionally
/// stable and effectively exact for a ≈ 1 — the standard "frozen cell"
/// immersed-solid treatment (the pressure projection's small divergence
/// leak is the accepted trade-off of the method).  `a` in [0, 1];
/// a = 0.9–0.99 recommended.
bool apply_kinematic_freeze(const FDM3Config& config,
                             const std::vector<ImmersedSolid>& solids,
                             double blend,
                             std::span<double> u,
                             std::span<double> v,
                             std::span<double> w,
                             exd::engine::core::ModelStatus& status);

// ── Boussinesq buoyancy ──────────────────────────────

/// Add the Boussinesq buoyancy to the force arrays:
///   f += −β·(T − T_ref)·g   (m/s²)
/// `temperature` is the per-cell temperature (K), typically sampled from the
/// thermal module's channel; `gravity` the world gravity (e.g. {0,0,-9.81}).
bool add_boussinesq_forces(const FDM3Config& config,
                           std::span<const double> temperature,
                           double beta, double t_ref,
                           const std::array<double, 3>& gravity,
                           std::span<double> fx,
                           std::span<double> fy,
                           std::span<double> fz,
                           exd::engine::core::ModelStatus& status);

/// Sample a temperature channel at the cell centers into a per-cell field.
/// The channel's lattice should match the fdm3 cell lattice (the CHT-lite
/// pattern); out-of-bounds samples are an error.
bool sample_temperature_field(const FDM3Config& config,
                              const coupling::IScalarField3D& channel,
                              std::vector<double>& temperature,
                              exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::fdm3
