#pragma once

// ─────────────────────────────────────────────────────
// Species transport (implementation_spec §11 "Species
// Transport", §66 "Reaction-Diffusion").
//
// Multispecies advection-diffusion-reaction on a
// structured grid, BUILT ON THE CORE RUNTIME: mesh
// topology, FDM operators, matrix-free CG (§35) — the
// composable-operator pattern physics modules use from
// now on.  Per species s:
//
//   ∂c_s/∂t + u·∇c_s = ∇·(D_s ∇c_s) − k_s·c_s + Σ (conversions)
//
// Time stepping is operator-split per step:
//   1. explicit first-order upwind advection (fdm::upwind_advect)
//   2. implicit diffusion (I − dt·D_s·Δ) via CG — unconditionally
//      stable, SPD operator composed from FdmLaplacianOperator
//   3. exact first-order decay / pairwise A→B conversion
//
// Boundary faces: zero-flux mirror by default; optional fixed
// (Dirichlet) face values per species, re-pinned after each pass.
//
// Coupling: arbitrary velocity via a coupling::IVectorField3D
// channel (fdm3 flow, thermal…), concentration channels exported
// for downstream modules (particles, reaction, presets).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/coupling/field_channels.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace exd::engine::physics::species {

/// First-order pair conversion A → B (exact exponential pair treatment).
struct ConversionSpec
{
    size_t from = 0;
    size_t to = 1;
    double rate = 0.0;    // 1/s
};

/// Fixed (Dirichlet) face value per species.
struct FaceValueSpec
{
    mesh::BoundaryId face = mesh::BoundaryId::XNeg;
    bool fixed = false;
    std::vector<double> value;   // one per species (only first species used if shorter)
};

struct SpeciesConfig
{
    mesh::StructuredGrid grid;                      // node lattice (spec §9)
    std::vector<std::string> species = {"A"};       // names
    std::vector<double> diffusivity = {0.0};        // D_s (m²/s)
    std::vector<double> decay_rate = {0.0};         // k_s (1/s)
    std::vector<ConversionSpec> conversions;        // A → B reactions
    // 1 value per species (uniform) OR a full per-node field (spec §5 ICs)
    std::vector<std::vector<double>> initial_concentration = {{0.0}};

    // advection: uniform body velocity or sampled channel
    std::array<double, 3> body_velocity = {0.0, 0.0, 0.0};
    const coupling::IVectorField3D* velocity_channel = nullptr;

    std::vector<FaceValueSpec> boundary_faces;      // Dirichlet pins (others mirror)

    double dt = 0.001;                              // s (advective CFL clamped internally)
    uint64_t max_steps = 100000;
    bool steady = false;                            // iterate until max change < steady_tolerance
    double steady_tolerance = 1e-10;
};

struct SpeciesResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;
    std::vector<mesh::StructuredScalarGrid> concentration;  // per species, node-centered
    uint64_t steps = 0;
    double time = 0.0;
    std::vector<double> total_mass;                 // conservation ledger per species
    double max_change = 0.0;
};

/// Validate config; returns false with a diagnostic on fatal problems.
bool validate_species_config(const SpeciesConfig& config, exd::engine::core::ModelStatus& status);

/// Solve the species transport system (deterministic, no exceptions).
SpeciesResult solve_species(const SpeciesConfig& config);

/// Advance one time step (expert mode; state is provided by the caller).
/// `state` holds the per-species concentration grids (already configured).
bool step_species(const SpeciesConfig& config,
                  std::vector<mesh::StructuredScalarGrid>& state,
                  double dt, double& max_change_out,
                  exd::engine::core::ModelStatus& status);

/// Sampling channel for species `s` (spec §44 external consumer contract).
std::unique_ptr<coupling::IScalarField3D> make_concentration_channel(
    const std::vector<mesh::StructuredScalarGrid>& state, size_t s);

} // namespace exd::engine::physics::species
