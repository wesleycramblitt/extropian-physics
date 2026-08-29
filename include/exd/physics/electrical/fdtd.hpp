#pragma once

#include <exd/physics/model_status.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::electrical {

// ─────────────────────────────────────────────────────
// 3D FDTD (Yee lattice) — time-domain Maxwell.
// PEC box boundary, gaussian soft plane-wave source.
// Grid-friendly: same indexing machinery family as FDM.
// ─────────────────────────────────────────────────────

struct FdtdConfig
{
    // Yee cell counts (per axis)
    std::array<int32_t, 3> dims = {128, 16, 16};       // >= 3 per axis
    std::array<double, 3> spacing = {1e-3, 1e-3, 1e-3};// cell size (m)

    // Material (uniform)
    double eps_r = 1.0;
    double mu_r = 1.0;

    // Time step; 0 → auto: courant_factor × CFL (1/(c·sqrt(Σ1/dx²)))
    double dt = 0.0;
    double courant_factor = 0.99;

    // Soft plane-wave source: gaussian pulse injected into Ez at the
    // x-plane `source_plane_index`, uniform over y and z.
    int32_t source_plane_index = 8;
    double source_amplitude = 1.0;
    double source_t0 = 30.0;   // pulse center (in steps)
    double source_sigma = 6.0; // pulse width (in steps)

    int max_steps = 200;

    // Diagnostics
    bool record_energy = true; // per-step EM energy in history
};

/// Yee-lattice field data. Component arrays hold `dims` cells; convention
/// (half-index offsets implicit in the update loops):
///   Ex at (i+½, j,   k),   Ey at (i,   j+½, k),   Ez at (i,   j,   k+½)
///   Hx at (i,   j+½, k+½), Hy at (i+½, j,   k+½), Hz at (i+½, j+½, k)
struct FdtdField
{
    bool valid = false;
    std::array<int32_t, 3> dims = {0, 0, 0};
    std::vector<double> ex, ey, ez;   // V/m, each dims.x·dims.y·dims.z
    std::vector<double> hx, hy, hz;   // A/m, each dims.x·dims.y·dims.z
    uint64_t step = 0;
    double t = 0.0;
};

struct FdtdStepResult
{
    double max_e = 0.0;   // max |E| (V/m)
    double max_h = 0.0;   // max |H| (A/m)
    double energy = 0.0;  // Σ (ε|E|² + μ|H|²)·dV (J)
};

/// Allocate + zero a field for the config. Returns false and sets status on
/// invalid config (dims < 3, spacing <= 0, eps/mu <= 0).
bool init_fdtd_field(const FdtdConfig& config, FdtdField& field,
                     exd::physics::ModelStatus& status);

/// Advance one Yee step in place (PEC boundary, source injection, curl
/// updates). Returns false and sets status on invalid field/config.
bool step_fdtd(const FdtdConfig& config, FdtdField& field,
               FdtdStepResult& step_result, exd::physics::ModelStatus& status);

/// Run a full simulation: init then max_steps steps.
struct FdtdResult
{
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    FdtdField field;
    std::vector<FdtdStepResult> history;
};

FdtdResult run_fdtd(const FdtdConfig& config);

} // namespace exd::physics::electrical
