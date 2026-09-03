#pragma once

// Private header for FDM CFD solver internals.
// NOT visible to downstream consumers (PRIVATE include path in CMake).

#include <exd/engine/physics/fluid/fdm/fdm_config.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_result.hpp>

#include <vector>
#include <memory>

namespace exd::engine::physics::fluid::fdm {

// ─────────────────────────────────────────────────────
// Collocated grid fields
// ─────────────────────────────────────────────────────
//
// All variables (u, v, p) live at cell centers.
// Grid: (nx+2) x (ny+2) including one layer of ghost cells.
// Physical cells: i=1..nx, j=1..ny
// Ghost cells: i=0, i=nx+1, j=0, j=ny+1
//
// Access: data[i + j * stride]

struct FDMGrid {
    int nx = 0;
    int ny = 0;
    double dx = 0.0;
    double dy = 0.0;

    // Current state
    std::vector<double> u;  // x-velocity at cell centers
    std::vector<double> v;  // y-velocity at cell centers
    std::vector<double> p;  // pressure at cell centers

    // Temporary storage
    std::vector<double> u_old, v_old, p_old;
    std::vector<double> u_tmp, v_tmp;
    std::vector<double> p_prime;  // pressure correction
    std::vector<double> rhs;      // RHS for pressure Poisson

    int stride() const { return nx + 2; }
    size_t idx(int i, int j) const { return static_cast<size_t>(i + j * stride()); }

    void allocate(int nx_, int ny_);
    void initialize(double u0, double v0, double p0);
};

// ─────────────────────────────────────────────────────
// Spatial discretization (central differences on collocated grid)
// ─────────────────────────────────────────────────────

namespace spatial {

inline double dudx(const FDMGrid& g, int i, int j) {
    return (g.u[g.idx(i+1, j)] - g.u[g.idx(i-1, j)]) / (2.0 * g.dx);
}

inline double dudy(const FDMGrid& g, int i, int j) {
    return (g.u[g.idx(i, j+1)] - g.u[g.idx(i, j-1)]) / (2.0 * g.dy);
}

inline double dvdx(const FDMGrid& g, int i, int j) {
    return (g.v[g.idx(i+1, j)] - g.v[g.idx(i-1, j)]) / (2.0 * g.dx);
}

inline double dvdy(const FDMGrid& g, int i, int j) {
    return (g.v[g.idx(i, j+1)] - g.v[g.idx(i, j-1)]) / (2.0 * g.dy);
}

inline double d2udx2(const FDMGrid& g, int i, int j) {
    return (g.u[g.idx(i+1, j)] - 2.0*g.u[g.idx(i, j)] + g.u[g.idx(i-1, j)]) / (g.dx*g.dx);
}

inline double d2udy2(const FDMGrid& g, int i, int j) {
    return (g.u[g.idx(i, j+1)] - 2.0*g.u[g.idx(i, j)] + g.u[g.idx(i, j-1)]) / (g.dy*g.dy);
}

inline double d2vdx2(const FDMGrid& g, int i, int j) {
    return (g.v[g.idx(i+1, j)] - 2.0*g.v[g.idx(i, j)] + g.v[g.idx(i-1, j)]) / (g.dx*g.dx);
}

inline double d2vdy2(const FDMGrid& g, int i, int j) {
    return (g.v[g.idx(i, j+1)] - 2.0*g.v[g.idx(i, j)] + g.v[g.idx(i, j-1)]) / (g.dy*g.dy);
}

inline double dpdx(const FDMGrid& g, int i, int j) {
    return (g.p[g.idx(i+1, j)] - g.p[g.idx(i-1, j)]) / (2.0 * g.dx);
}

inline double dpdy(const FDMGrid& g, int i, int j) {
    return (g.p[g.idx(i, j+1)] - g.p[g.idx(i, j-1)]) / (2.0 * g.dy);
}

inline double divergence(const FDMGrid& g, int i, int j) {
    return dudx(g, i, j) + dvdy(g, i, j);
}

inline double vorticity(const FDMGrid& g, int i, int j) {
    return dvdx(g, i, j) - dudy(g, i, j);
}

} // namespace spatial

// ─────────────────────────────────────────────────────
// Advection operators
// ─────────────────────────────────────────────────────

namespace advection {

double convective_u(const FDMGrid& g, int i, int j, AdvectionScheme scheme);
double convective_v(const FDMGrid& g, int i, int j, AdvectionScheme scheme);

} // namespace advection

// ─────────────────────────────────────────────────────
// Time integration
// ─────────────────────────────────────────────────────

void compute_rhs(const FDMGrid& g, double nu, AdvectionScheme scheme,
                 std::vector<double>& du_out, std::vector<double>& dv_out);

void integrate_forward_euler(FDMGrid& g, double dt, double nu, AdvectionScheme scheme);
void integrate_heun(FDMGrid& g, double dt, double nu, AdvectionScheme scheme);
void integrate_rk4(FDMGrid& g, double dt, double nu, AdvectionScheme scheme);
void integrate_crank_nicolson(FDMGrid& g, double dt, double nu, AdvectionScheme scheme);

// ─────────────────────────────────────────────────────
// Pressure-velocity coupling
// ─────────────────────────────────────────────────────

void compute_pressure_rhs(const FDMGrid& g, double dt, std::vector<double>& rhs);
int solve_pressure_poisson(FDMGrid& g, const std::vector<double>& rhs,
                           int max_iterations, double tolerance, double omega);

// ─────────────────────────────────────────────────────
// Boundary conditions
// ─────────────────────────────────────────────────────

void apply_boundary_conditions(FDMGrid& g, const FDMConfig& config);

// ─────────────────────────────────────────────────────
// Utility
// ─────────────────────────────────────────────────────

double max_velocity(const FDMGrid& g);
double max_divergence(const FDMGrid& g);

} // namespace exd::engine::physics::fluid::fdm
