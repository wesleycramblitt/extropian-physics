#include <cstdint>
// Time stepping: adaptive, CFL-based, and multi-rate stepping.

#include <cstdio>
#include <algorithm>

namespace exd::physics {

struct TimeStepper {
    double t = 0.0;           // current simulation time
    double dt = 0.001;        // current time step
    double dt_min = 1e-8;
    double dt_max = 1.0;
    double cfl_target = 0.8;  // target CFL number
    double cfl_current = 0.0; // CFL from last step
    uint64_t step_count = 0;

    /// Compute CFL number for uniform grid
    double compute_cfl_u(double u_max, double dx) const {
        return u_max * dt / dx;
    }

    /// Compute CFL number for diffusion
    double compute_cfl_nu(double nu, double dx) const {
        return nu * dt / (dx * dx);
    }

    /// Adjust dt to hit target CFL
    void adapt(double u_max, double nu, double dx) {
        double cfl_adv = compute_cfl_u(u_max, dx);
        double cfl_dif = compute_cfl_nu(nu, dx);
        double cfl = std::max(cfl_adv, cfl_dif);
        cfl_current = cfl;

        if (cfl > 0 && cfl_target > 0) {
            double ratio = cfl_target / cfl;
            dt = std::clamp(dt * ratio, dt_min, dt_max);
        }
    }

    void advance() { t += dt; step_count++; }
};

struct ConvergenceMonitor {
    double tolerance = 1e-6;
    double residual = 1.0;
    bool converged = false;
    uint64_t max_iterations = 1000000;

    bool check(double current_residual) {
        residual = current_residual;
        converged = residual < tolerance;
        return converged;
    }
};

} // namespace exd::physics
