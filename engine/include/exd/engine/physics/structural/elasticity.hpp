#pragma once

// ---------------------------------------------------------------------------
// Static linear elasticity (Navier-Cauchy) solver on a regular node-centered
// 3D grid, displacement form, small strain.
//
//   div sigma + rho*f = 0
//   sigma = lambda*tr(eps - eps_th)*I + 2*mu*(eps - eps_th)
//   eps   = 0.5*(grad u + grad u^T)
//
// `body_force` is a SPECIFIC force (acceleration, m/s^2, gravity-style).
// Optional isotropic thermal strain eps_th = alpha*(T - T_ref)*I is driven by
// the temperature channel passed to solve_elasticity(); T_ref is sampled at
// the grid origin node, so for vertical columns "bottom-anchored" behavior is
// obtained with the channel giving T = T_ref + g*z (delta T = g*z).
//
// Numerics: central second differences (19-point stencil) with the 4-point
// mixed-derivative stencil and one node-wide ghost layer.  Free surfaces use
// the standard central-difference Robin (ghost-mirror) treatment of the
// traction condition sigma*n = t.  Only the +z face carries a configurable
// (optionally regionally masked) traction; every other face is traction free.
// SOR (omega = 1.5) solves the three displacement components component-wise
// (Gauss-Seidel within each component block in a single sweep).
//
// Deterministic and exception-free: every failure surfaces through
// ModelStatus.  Dirichlet nodes are pinned through a per-node, per-component
// mask so roller supports (pin only u_x/u_y) are expressible.
// ---------------------------------------------------------------------------

#include <exd/engine/coupling/field_channels.hpp>
#include <exd/engine/core/model_status.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::structural {

/// Node-centered regular 3D grid for the displacement solve.
struct ElasticityGridConfig
{
    std::array<double, 3> origin = {0.0, 0.0, 0.0};     // coordinate of node (0,0,0) in m
    std::array<double, 3> spacing = {0.05, 0.05, 0.05}; // node spacing per axis in m, > 0
    std::array<int32_t, 3> dims = {0, 0, 0};            // node counts, >= 2 per axis
};

/// Isotropic linear-elastic material.
struct ElasticMaterialConfig
{
    double elastic_modulus = 200e9; // E in Pa, > 0
    double poisson_ratio = 0.3;     // nu in (-1, 0.5)
};

/// Full elastic-static configuration.
struct ElasticityConfig
{
    ElasticityGridConfig grid;
    ElasticMaterialConfig material;

    std::array<double, 3> body_force = {0.0, 0.0, 0.0}; // specific body force (m/s^2)

    // Per-node, per-component Dirichlet mask (index i + nx*(j + ny*k)).
    // fixed_mask[idx][c] == true pins component c of node idx at
    // fixed_displacement[idx][c] (zero displacement when fixed_displacement
    // is empty).  Empty = all nodes free (rigid modes must then be removed
    // by pinning at least one node).
    std::vector<std::array<bool, 3>> fixed_mask;

    // Per-node fixed displacement in m.  Empty = zeros for every pinned node.
    // Size must match the node count when non-empty.
    std::vector<std::array<double, 3>> fixed_displacement;

    std::array<double, 3> surface_traction = {0.0, 0.0, 0.0}; // traction on the +z face in Pa

    // Regional +z-face traction selector: traction_mask[idx] == true marks
    // nodes on the +z face that carry surface_traction.  Empty = the whole
    // +z face carries surface_traction.
    std::vector<bool> traction_mask;

    double thermal_expansion_coefficient = 0.0; // alpha in 1/K, >= 0 (0 disables thermal strain)

    double tolerance = 1e-10;         // SOR relative-residual target
    uint64_t max_iterations = 200000; // SOR sweep limit, > 0
};

/// Validate the elastic config.  Fatal problems return false and fill
/// `error`; non-fatal observations are appended to `warnings`.
bool validate_elasticity_config(const ElasticityConfig& config,
                                std::string& error,
                                std::vector<std::string>& warnings);

/// Result of an elastic solve.
struct ElasticityResult
{
    bool ok = false;
    ModelStatus status;

    exd::engine::coupling::StructuredVectorGrid displacement; // 3*N flat, values[3*idx + c], idx = i + nx*(j + ny*k)

    double max_residual = 0.0;          // best SOR relative residual achieved
    uint64_t iterations = 0;            // SOR sweeps performed
    double strain_energy = 0.0;         // sum over cells of 0.5*sigma:eps*cell_volume in J
    double max_effective_strain = 0.0;  // max von-Mises equivalent strain over nodes
};

/// Solve the static displacement field.  `temperature_channel` is optional:
/// when it is non-null and thermal_expansion_coefficient > 0, thermal strain
/// eps_th = alpha*(T - T_ref)*I is added; otherwise thermal terms are inert
/// (warned).
ElasticityResult solve_elasticity(const ElasticityConfig& config,
                                  ModelStatus& status,
                                  const exd::engine::coupling::IScalarField3D* temperature_channel = nullptr);

} // namespace exd::engine::physics::structural