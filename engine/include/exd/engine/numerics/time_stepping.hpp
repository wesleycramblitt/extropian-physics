#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace exd::engine::numerics {

// ─────────────────────────────────────────────────────
// Shared time-stepping infrastructure, promoted from the
// private implementation so every solver and domain uses
// one adaptive-stepping policy instead of its own dt logic.
// ─────────────────────────────────────────────────────

/// Adaptive, CFL-based time stepper (advection + diffusion).
struct TimeStepper
{
    double t = 0.0;           // current simulation time (s)
    double dt = 0.001;        // current time step (s)
    double dt_min = 1e-8;     // step clamp lower bound (s)
    double dt_max = 1.0;      // step clamp upper bound (s)
    double cfl_target = 0.8;  // target CFL number
    double cfl_current = 0.0; // CFL of the last adapt() call
    uint64_t step_count = 0;  // steps taken

    /// CFL for advection on a uniform grid.
    [[nodiscard]] double compute_cfl_u(double u_max, double dx) const
    {
        return u_max * dt / dx;
    }

    /// CFL for diffusion on a uniform grid.
    [[nodiscard]] double compute_cfl_nu(double nu, double dx) const
    {
        return nu * dt / (dx * dx);
    }

    /// Adjust dt to hit `cfl_target` (clamped to [dt_min, dt_max]).
    void adapt(double u_max, double nu, double dx)
    {
        const double cfl_adv = compute_cfl_u(u_max, dx);
        const double cfl_dif = compute_cfl_nu(nu, dx);
        const double cfl = std::max(cfl_adv, cfl_dif);
        cfl_current = cfl;
        if (cfl > 0.0 && cfl_target > 0.0)
        {
            const double ratio = cfl_target / cfl;
            dt = std::clamp(dt * ratio, dt_min, dt_max);
        }
    }

    /// Advance one step.
    void advance()
    {
        t += dt;
        ++step_count;
    }
};

/// Windowed residual-based convergence monitor.
struct ConvergenceMonitor
{
    double tolerance = 1e-6;   // residual threshold
    double residual = 1.0;     // latest residual
    bool converged = false;    // residual < tolerance
    uint64_t max_iterations = 1000000;

    /// Feed the latest residual; returns whether convergence was reached.
    bool check(double current_residual)
    {
        residual = current_residual;
        converged = residual < tolerance;
        return converged;
    }
};

} // namespace exd::engine::numerics
