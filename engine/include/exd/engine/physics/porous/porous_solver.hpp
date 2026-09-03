#pragma once

// ─────────────────────────────────────────────────────
// Porous media / Darcy flow (implementation_spec §11
// "Porous Media").
//
// Single-phase slightly-compressible pressure diffusion
// on a structured grid:
//
//   φ·c_t·∂p/∂t = ∇·((k/μ)·∇p) + q
//
// with hydraulic diffusivity K = k/(φ·μ·c_t), well
// sources q (1/s per unit volume), no-flow (mirror) faces
// by default and optional fixed-pressure (Dirichlet)
// faces.  Built on the core runtime: implicit diffusion
// via the shared DiffusionStepOperator + CG.
//
// Coupling: pressure channel (IScalarField3D) and Darcy
// velocity channel (IVectorField3D, v = −(k/μ)∇p) for
// particle/field coupling.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/coupling/field_channels.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>
#include <memory>
#include <vector>

namespace exd::engine::physics::porous {

struct PorousConfig
{
    mesh::StructuredGrid grid;

    // rock/fluid properties
    double permeability = 1e-12;        // k (m²)
    double viscosity = 1e-3;            // μ (Pa·s)
    double porosity = 0.2;              // φ (−)
    double compressibility = 1e-9;      // c_t (1/Pa) → K = k/(φ·μ·c_t) (m²/s)

    // well sources: injection rate per unit volume at node idx (1/s)
    std::vector<double> source_rate;    // node count; empty = none

    // fixed-pressure faces (Dirichlet pins); others are no-flow mirrors
    struct FacePressure
    {
        mesh::BoundaryId face = mesh::BoundaryId::XNeg;
        bool fixed = false;
        double value = 0.0;             // Pa
    };
    std::vector<FacePressure> boundary_faces;

    double initial_pressure = 0.0;      // Pa (uniform)
    std::vector<double> initial_pressure_field;  // per-node field; empty = uniform
    double dt = 0.01;                   // s (implicit → no stability limit)
    uint64_t max_steps = 10000;
    bool steady = false;                // iterate until max change < tolerance
    double steady_tolerance = 1e-10;
};

struct PorousResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;
    mesh::StructuredScalarGrid pressure;     // node-centered (Pa)
    uint64_t steps = 0;
    double time = 0.0;
    double max_change = 0.0;
    double total_mass = 0.0;                 // Σ p·(φ·c_t)·cell_vol (kg)
};

/// Validate the porous config.
bool validate_porous_config(const PorousConfig& config, exd::engine::core::ModelStatus& status);

/// Solve the pressure diffusion problem (deterministic, no exceptions).
PorousResult solve_porous(const PorousConfig& config);

/// Pressure sampling channel (spec §44).
std::unique_ptr<coupling::IScalarField3D> make_pressure_channel(const PorousResult& result);

/// Darcy velocity channel: v = −(k/μ)·∇p (node-centered, central differences).
std::unique_ptr<coupling::IVectorField3D> make_darcy_velocity_channel(const PorousConfig& config,
                                                                      const PorousResult& result);

/// Hydraulic diffusivity K = k/(φ·μ·c_t) (m²/s).
inline double hydraulic_diffusivity(const PorousConfig& config)
{
    return config.permeability /
           (config.porosity * config.viscosity * config.compressibility);
}

} // namespace exd::engine::physics::porous
