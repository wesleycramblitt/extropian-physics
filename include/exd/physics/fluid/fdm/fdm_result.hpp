#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::fluid::fdm {

// ─────────────────────────────────────────────────────
// Per-step snapshot
// ─────────────────────────────────────────────────────

struct FDMStepResult {
    double time = 0.0;
    int step = 0;
    double residual_u = 0.0;
    double residual_v = 0.0;
    double residual_p = 0.0;
    double max_velocity = 0.0;
    double cfl = 0.0;
    double drag_coefficient = 0.0;
    double lift_coefficient = 0.0;
};

// ─────────────────────────────────────────────────────
// Field data (flat arrays, row-major)
// ─────────────────────────────────────────────────────

struct FDMFieldData {
    int nx = 0;
    int ny = 0;

    // Velocity at cell centers (interpolated from staggered locations)
    std::vector<double> u;      // size nx*ny, row-major [i + j*nx]
    std::vector<double> v;      // size nx*ny

    // Pressure at cell centers
    std::vector<double> p;      // size nx*ny

    // Vorticity at cell centers (optional)
    std::vector<double> omega;  // size nx*ny, empty if not computed

    // Grid coordinates at cell centers
    std::vector<double> x;      // size nx
    std::vector<double> y;      // size ny

    size_t index(int i, int j) const { return static_cast<size_t>(i + j * nx); }
};

// ─────────────────────────────────────────────────────
// Solver result
// ─────────────────────────────────────────────────────

struct FDMResult {
    bool valid = false;
    bool converged = false;
    std::string error;
    std::vector<std::string> warnings;

    // Final state
    FDMFieldData field;

    // Time history
    std::vector<FDMStepResult> history;

    // Summary
    int total_steps = 0;
    double final_time = 0.0;
    double final_residual = 0.0;
    double max_velocity = 0.0;
    double reynolds_number = 0.0;  // Based on domain length and max velocity
};

} // namespace exd::physics::fluid::fdm
