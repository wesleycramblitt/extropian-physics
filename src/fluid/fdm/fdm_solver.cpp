#include <exd/physics/fluid/fdm/fdm_solver.hpp>
#include "fdm_internal.hpp"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace exd::physics::fluid::fdm {

// ─────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────

double max_velocity(const FDMGrid& g) {
    double umax = 0.0;
    for (int j = 1; j <= g.ny; ++j) {
        for (int i = 1; i <= g.nx; ++i) {
            double mag = std::sqrt(g.u[g.idx(i, j)]*g.u[g.idx(i, j)] +
                                   g.v[g.idx(i, j)]*g.v[g.idx(i, j)]);
            if (mag > umax) umax = mag;
        }
    }
    return umax;
}

double max_divergence(const FDMGrid& g) {
    double dmax = 0.0;
    for (int j = 1; j <= g.ny; ++j) {
        for (int i = 1; i <= g.nx; ++i) {
            double div = std::abs(spatial::divergence(g, i, j));
            if (div > dmax) dmax = div;
        }
    }
    return dmax;
}

// ─────────────────────────────────────────────────────
// Field extraction
// ─────────────────────────────────────────────────────

static FDMFieldData extract_field(const FDMGrid& g, const FDMConfig& config,
                                  bool compute_vorticity) {
    FDMFieldData field;
    field.nx = g.nx;
    field.ny = g.ny;

    const double dx = config.dx();
    const double dy = config.dy();
    field.x.resize(g.nx);
    field.y.resize(g.ny);
    for (int i = 0; i < g.nx; ++i)
        field.x[i] = (i + 0.5) * dx;
    for (int j = 0; j < g.ny; ++j)
        field.y[j] = (j + 0.5) * dy;

    field.u.resize(g.nx * g.ny);
    field.v.resize(g.nx * g.ny);
    field.p.resize(g.nx * g.ny);
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            size_t id = field.index(i, j);
            field.u[id] = g.u[g.idx(i+1, j+1)];
            field.v[id] = g.v[g.idx(i+1, j+1)];
            field.p[id] = g.p[g.idx(i+1, j+1)];
        }
    }

    if (compute_vorticity) {
        field.omega.resize(g.nx * g.ny);
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i)
                field.omega[field.index(i, j)] = spatial::vorticity(g, i+1, j+1);
    }

    return field;
}

static void apply_time_integration(FDMGrid& g, double dt, double nu,
                                   AdvectionScheme scheme, TimeIntegration method) {
    switch (method) {
        case TimeIntegration::ForwardEuler: integrate_forward_euler(g, dt, nu, scheme); break;
        case TimeIntegration::Heun:         integrate_heun(g, dt, nu, scheme); break;
        case TimeIntegration::RK4:          integrate_rk4(g, dt, nu, scheme); break;
        case TimeIntegration::CrankNicolson:integrate_crank_nicolson(g, dt, nu, scheme); break;
    }
}

// ─────────────────────────────────────────────────────
// Validation
// ─────────────────────────────────────────────────────

static bool validate_config(const FDMConfig& config, std::string& error) {
    if (config.nx < 2 || config.ny < 2) { error = "Grid must be at least 2x2"; return false; }
    if (config.lx <= 0.0 || config.ly <= 0.0) { error = "Domain dimensions must be positive"; return false; }
    if (config.rho <= 0.0) { error = "Density must be positive"; return false; }
    if (config.mu < 0.0) { error = "Viscosity must be non-negative"; return false; }
    if (config.dt <= 0.0) { error = "Time step must be positive"; return false; }
    if (config.max_steps < 1) { error = "Maximum steps must be >= 1"; return false; }
    if (config.sor_omega < 1.0 || config.sor_omega >= 2.0) { error = "SOR omega must be in [1, 2)"; return false; }
    if (config.boundary_conditions.empty()) { error = "At least one BC required"; return false; }
    return true;
}

// ─────────────────────────────────────────────────────
// Public API: solve_fdm
// ─────────────────────────────────────────────────────

FDMResult solve_fdm(const FDMConfig& config) {
    FDMResult result;

    std::string error;
    if (!validate_config(config, error)) {
        result.valid = false;
        result.error = error;
        return result;
    }

    const double nu = config.kinematic_viscosity();
    const double dx = config.dx();
    const double dy = config.dy();

    FDMGrid grid;
    grid.allocate(config.nx, config.ny);
    grid.dx = dx;
    grid.dy = dy;
    grid.initialize(config.initial_u, config.initial_v, config.initial_p);
    apply_boundary_conditions(grid, config);

    double t = 0.0;
    bool converged = false;

    for (int step = 0; step < config.max_steps; ++step) {
        double dt = config.dt;

        // Save state for SIMPLE
        std::copy(grid.u.begin(), grid.u.end(), grid.u_old.begin());
        std::copy(grid.v.begin(), grid.v.end(), grid.v_old.begin());
        std::copy(grid.p.begin(), grid.p.end(), grid.p_old.begin());

        // 1. Time integration (predictor step)
        apply_time_integration(grid, dt, nu, config.advection_scheme, config.time_integration);

        // 2. Pressure-velocity coupling (SIMPLE)
        // Compute pressure correction from divergence of predicted velocity
        compute_pressure_rhs(grid, dt, grid.rhs);

        // Solve for p' with zero initial guess
        std::fill(grid.p.begin(), grid.p.end(), 0.0);
        solve_pressure_poisson(grid, grid.rhs, config.pressure_max_iterations,
                               config.pressure_tolerance, config.sor_omega);
        // Now grid.p contains p'

        // Copy p' to p_prime, restore old pressure
        std::copy(grid.p.begin(), grid.p.end(), grid.p_prime.begin());
        std::copy(grid.p_old.begin(), grid.p_old.end(), grid.p.begin());

        // Correct velocity: u = u* - alpha_u * dt * grad(p')
        for (int j = 1; j <= grid.ny; ++j) {
            for (int i = 1; i <= grid.nx; ++i) {
                size_t id = grid.idx(i, j);
                double dpdx = (grid.p_prime[grid.idx(i+1, j)] - grid.p_prime[grid.idx(i-1, j)]) / (2.0 * dx);
                double dpdy = (grid.p_prime[grid.idx(i, j+1)] - grid.p_prime[grid.idx(i, j-1)]) / (2.0 * dy);
                grid.u[id] = grid.u_old[id] + config.velocity_under_relaxation *
                    (grid.u[id] - grid.u_old[id] - dt * dpdx);
                grid.v[id] = grid.v_old[id] + config.velocity_under_relaxation *
                    (grid.v[id] - grid.v_old[id] - dt * dpdy);
            }
        }

        // Update pressure: p = p_old + alpha_p * p'
        for (int j = 1; j <= grid.ny; ++j)
            for (int i = 1; i <= grid.nx; ++i) {
                size_t id = grid.idx(i, j);
                grid.p[id] = grid.p_old[id] + config.pressure_under_relaxation * grid.p_prime[id];
            }

        // 3. Apply boundary conditions
        apply_boundary_conditions(grid, config);

        t += dt;

        // Diagnostics
        double umax = max_velocity(grid);
        double div_max = max_divergence(grid);

        FDMStepResult step_result;
        step_result.time = t;
        step_result.step = step + 1;
        step_result.max_velocity = umax;
        step_result.cfl = umax * dt / dx;
        step_result.residual_p = div_max;
        result.history.push_back(step_result);

        if (step > 0 && step % config.convergence_window == 0) {
            if (div_max < config.convergence_tolerance) {
                converged = true;
                break;
            }
        }
    }

    result.field = extract_field(grid, config, config.output_vorticity);
    result.valid = true;
    result.converged = converged;
    result.total_steps = static_cast<int>(result.history.size());
    result.final_time = t;
    result.max_velocity = max_velocity(grid);
    result.reynolds_number = result.max_velocity * config.lx / nu;
    result.final_residual = max_divergence(grid);

    if (!converged)
        result.warnings.push_back("Solver did not converge within max_steps");

    return result;
}

FDMStepResult step_fdm(const FDMConfig& config, double& t, int& step, FDMFieldData& field) {
    FDMStepResult result;
    const double nu = config.kinematic_viscosity();
    const double dx = config.dx();

    FDMGrid grid;
    grid.allocate(config.nx, config.ny);
    grid.dx = dx;
    grid.dy = config.dy();
    grid.initialize(config.initial_u, config.initial_v, config.initial_p);

    // Time integration + SIMPLE
    apply_time_integration(grid, config.dt, nu, config.advection_scheme, config.time_integration);

    compute_pressure_rhs(grid, config.dt, grid.rhs);
    std::fill(grid.p.begin(), grid.p.end(), 0.0);
    solve_pressure_poisson(grid, grid.rhs, config.pressure_max_iterations,
                           config.pressure_tolerance, config.sor_omega);
    std::copy(grid.p.begin(), grid.p.end(), grid.p_prime.begin());

    for (int j = 1; j <= grid.ny; ++j)
        for (int i = 1; i <= grid.nx; ++i) {
            size_t id = grid.idx(i, j);
            double dpdx = (grid.p_prime[grid.idx(i+1, j)] - grid.p_prime[grid.idx(i-1, j)]) / (2.0 * dx);
            double dpdy = (grid.p_prime[grid.idx(i, j+1)] - grid.p_prime[grid.idx(i, j-1)]) / (2.0 * config.dy());
            grid.u[id] -= config.velocity_under_relaxation * config.dt * dpdx;
            grid.v[id] -= config.velocity_under_relaxation * config.dt * dpdy;
        }

    apply_boundary_conditions(grid, config);
    t += config.dt;
    step++;
    field = extract_field(grid, config, config.output_vorticity);

    result.time = t;
    result.step = step;
    result.max_velocity = max_velocity(grid);
    result.cfl = result.max_velocity * config.dt / dx;
    result.residual_p = max_divergence(grid);
    return result;
}

} // namespace exd::physics::fluid::fdm
