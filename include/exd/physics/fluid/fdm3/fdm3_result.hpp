#pragma once

// Result and field types for the 3D incompressible FDM CFD solver (fdm3).

#include <cstddef>
#include <string>
#include <vector>

namespace exd::physics::fluid::fdm3 {

// ─────────────────────────────────────────────────────
// Field data (flat arrays at cell centers, row-major)
// ─────────────────────────────────────────────────────

struct FDM3FieldData {
    int nx = 0;
    int ny = 0;
    int nz = 0;

    // Velocity and pressure at cell centers.
    // Each array has size nx*ny*nz, flat index i + nx*(j + ny*k).
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> p;

    // Grid coordinates at cell centers, e.g. x[i] = (i + 0.5) * dx.
    std::vector<double> x;      // size nx
    std::vector<double> y;      // size ny
    std::vector<double> z;      // size nz

    size_t index(int i, int j, int k) const {
        return static_cast<size_t>(i) + static_cast<size_t>(nx) *
               (static_cast<size_t>(j) + static_cast<size_t>(ny) * static_cast<size_t>(k));
    }
};

// ─────────────────────────────────────────────────────
// Per-step snapshot
// ─────────────────────────────────────────────────────

struct FDM3StepResult {
    double time = 0.0;          // Simulated time after this step (s)
    int step = 0;               // Step number (1-based)
    double residual_u = 0.0;    // Max |u - u_old| during the step (m/s)
    double residual_v = 0.0;    // Max |v - v_old| during the step (m/s)
    double residual_w = 0.0;    // Max |w - w_old| during the step (m/s)
    double residual_p = 0.0;    // Max |p - p_old| during the step (Pa)
    double max_velocity = 0.0;  // Max |(u,v,w)| over interior cells (m/s)
    double cfl = 0.0;           // CFL = max_velocity * dt / min(dx,dy,dz)
    double divergence = 0.0;    // Max |div u| over interior cells (1/s)
};

// ─────────────────────────────────────────────────────
// Solver result
// ─────────────────────────────────────────────────────

struct FDM3Result {
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;

    // Final state
    FDM3FieldData field;

    // Time history
    std::vector<FDM3StepResult> history;

    // Summary
    int steps_taken = 0;
    double final_time = 0.0;
    bool converged = false;
};

} // namespace exd::physics::fluid::fdm3