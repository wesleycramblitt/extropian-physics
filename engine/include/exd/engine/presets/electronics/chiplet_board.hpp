#pragma once

// ─────────────────────────────────────────────────────
// Chiplet-board preset (chiplet use case).
//
// A heterogeneous-conduction thermal board: multiple chip
// power sources, optional high-k heat-spreader regions,
// isothermal sink faces.  Built ON THE CORE RUNTIME: the
// steady heterogeneous heat equation
//   −∇·(k(x)·∇T) = q(x)
// with PER-NODE conductivity k (face value = harmonic
// mean for series conduction), per-node source q, and
// Dirichlet sink faces.  The operator is affine in the
// sink values; the linear part + CG solve it directly
// (the porous-media pattern, spec §35).
//
// Verified: bi-material strip reproduces the exact
// piecewise-linear profile; total chip power equals the
// sink flux (energy balance); the hottest point sits
// under the highest-power chip.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <vector>

namespace exd::engine::presets::electronics {

struct ChipletBoardConfig
{
    mesh::StructuredGrid grid;             // board lattice (x-y plane, thin z)

    double board_thickness = 0.001;        // (m) — source density = W/(A·t)
    double sink_temperature = 300.0;       // (K) on the board perimeter (all
                                           // six faces Dirichlet at this value)

    struct Chip
    {
        double x = 0.0;                    // center (m)
        double y = 0.0;
        int w_cells = 1;                   // footprint (cells)
        int h_cells = 1;
        double power_watts = 0.0;          // total dissipated (W)
    };
    std::vector<Chip> chips;

    struct Spreader
    {
        double x = 0.0;                    // center (m)
        double y = 0.0;
        int w_cells = 1;
        int h_cells = 1;
        double conductivity = 400.0;       // W/(m·K)
    };
    std::vector<Spreader> spreaders;

    double base_conductivity = 20.0;       // board material (W/(m·K))
    double tolerance = 1e-10;
    int max_iterations = 20000;
};

struct ChipletBoardResult
{
    bool ok = false;
    core::ModelStatus status;
    mesh::StructuredScalarGrid temperature;  // K, node-centered
    double total_power = 0.0;                // ∫q dV (W)
    double sink_flux = 0.0;                  // k·A·Σ|∂T/∂n| over the faces (W)
    double peak_temperature = 0.0;           // K
    double peak_x = 0.0, peak_y = 0.0;
    double max_residual = 0.0;
};

/// Solve the heterogeneous-conduction board (deterministic, no exceptions).
ChipletBoardResult solve_chiplet_board(const ChipletBoardConfig& config);

} // namespace exd::engine::presets::electronics
