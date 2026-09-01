#pragma once

// ---------------------------------------------------------------------
// Phase I acoustics domain: scalar wave equation on a
// regular structured grid.
//
//   d^2 p/dt^2 = c^2 * laplace(p)
//
// Boundary treatment: pressure-release (soft) walls p = 0
// on all boundary nodes for axes with 3+ nodes, which makes
// the box modes
//   p0 = A*sin(l*pi*x/Lx)*sin(m*pi*y/Ly)*sin(n*pi*z/Lz)
// exact eigenmodes of the discrete Dirichlet problem and
// enables the analytic box-mode verification in the tests.
// Axes with exactly 2 nodes have no interior; their two
// faces are treated as zero-normal-derivative symmetry
// planes (Neumann mirror) so a plane wave uniform across
// that thin dimension is supported (used by the 1D plane-
// wave test).  This keeps the operator deterministic.
//
// Time integration: explicit leapfrog
//   p^{n+1} = 2*p^n - p^{n-1} + (c*dt)^2 * laplace(p^n)
// started with p(-dt) = p(0) (zero initial velocity).
// When config.dt == 0 the step is CFL-adaptive:
//   dt = 0.8*dx_min/(c*sqrt(3))
// A requested dt above the exact stability limit
// c*dt*sqrt(1/dx^2 + 1/dy^2 + 1/dz^2) <= 1 is clamped to
// the limit and the caller gets a "CFL violation" warning.
//
// No exceptions, no I/O, pure deterministic free functions.
// ---------------------------------------------------------------------

#include <exd/physics/model_status.hpp>
#include <exd/physics/coupling/field_channels.hpp>   // StructuredScalarGrid (result payload only)

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::acoustics {

// ---------------------------------------------------------------------
// Grid configuration
// ---------------------------------------------------------------------

struct WaveGridConfig {
    std::array<double, 3> origin = {0.0, 0.0, 0.0};     // node (0,0,0) (m)
    std::array<double, 3> spacing = {0.05, 0.05, 0.05}; // node spacing per axis (m), > 0
    std::array<int32_t, 3> dims = {0, 0, 0};            // node counts (>= 2 per axis)
};

// ---------------------------------------------------------------------
// Solver configuration
// ---------------------------------------------------------------------

struct WaveConfig {
    WaveGridConfig grid;
    double sound_speed = 343.0;       // m/s
    // Uniform mean flow (m/s): solves the linearized convected wave equation
    // (dt + u·grad)^2 p = c^2 laplace p (aeroacoustic-lite). Zero = plain
    // wave equation. CFL is keyed on the effective speed c + |u|; a warning
    // fires when |u| >= c (no upstream propagation).
    std::array<double, 3> mean_flow = {0.0, 0.0, 0.0};
    double dt = 0.0;                  // s; 0 -> CFL-adaptive (0.8*dx_min/(c*sqrt(3)))
    uint64_t max_steps = 10000;       // leapfrog step cap
    int32_t probe_index = 0;          // flat node index, i + nx*(j + ny*k)
    double amplitude = 1.0;           // initial-pressure amplitude A

    // Box-mode indices of the seeded eigenmode, default the fundamental
    // (1,1,1).  Affects only the initial condition:
    //   p0 = A*sin(l*pi*x/Lx)*sin(m*pi*y/Ly)*sin(n*pi*z/Lz)
    // with the constant 1 factor on any 2-node (Neumann) axis.
    std::array<int32_t, 3> initial_mode = {1, 1, 1};
    // Optional arbitrary initial pressure field (flat, node count).  Empty =
    // the box-mode seed (initial_mode).  Zero initial velocity either way.  A
    // localized pulse here is what the mean-flow tests use.
    std::vector<double> initial_pressure = {};
};

// ---------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------

struct WaveResult {
    bool ok = false;
    ModelStatus status;
    coupling::StructuredScalarGrid pressure_final;  // node-centered pressure after the run
    std::vector<double> probe_history;              // per-step probe samples (element 0 = t0)
    double dt_used = 0.0;                           // step actually used after adapt/clamp
    double max_pressure = 0.0;                      // max |p| over all nodes across the run
    double energy_initial = 0.0;                    // sum p^2 * cell_volume at t = 0
    double energy_final = 0.0;                      // sum p^2 * cell_volume at the final step
};

// ---------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------

/// Validate a wave config.  Fatal problems return false and set `error`;
/// non-fatal problems (e.g. a CFL-violating requested dt) are appended to
/// `warnings`.
bool validate_wave_config(const WaveConfig& config,
                          std::string& error,
                          std::vector<std::string>& warnings);

/// Solve the scalar wave equation.  On failure (`status.ok == false`) the
/// returned result has `ok == false` and an empty field.  Deterministic.
WaveResult solve_wave(const WaveConfig& config, ModelStatus& status);

} // namespace exd::physics::acoustics