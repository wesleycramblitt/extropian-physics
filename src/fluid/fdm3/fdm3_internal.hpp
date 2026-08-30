#pragma once

// Private header for the 3D FDM CFD solver (fdm3) internals.
// NOT visible to downstream consumers (PRIVATE include path in CMake).

#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>

#include <vector>

namespace exd::physics::fluid::fdm3 {

// ─────────────────────────────────────────────────────
// Collocated 3D grid fields
// ─────────────────────────────────────────────────────
//
// All variables (u, v, w, p) live at cell centers.
// Grid: (nx+2) x (ny+2) x (nz+2) including one layer of ghost cells.
// Physical cells: i=1..nx, j=1..ny, k=1..nz
// Ghost cells:    i=0, i=nx+1, j=0, j=ny+1, k=0, k=nz+1
//
// Flat, row-major access: data[idx(i,j,k)] with
//   idx(i,j,k) = i + sx * (j + sy * k),   sx = nx+2, sy = ny+2

struct FDM3Grid {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    // Current state
    std::vector<double> u;  // x-velocity at cell centers
    std::vector<double> v;  // y-velocity at cell centers
    std::vector<double> w;  // z-velocity at cell centers
    std::vector<double> p;  // pressure at cell centers

    // Temporary storage
    std::vector<double> u_old, v_old, w_old, p_old;
    std::vector<double> p_prime;  // pressure correction
    std::vector<double> rhs;      // RHS for pressure Poisson

    int sx() const { return nx + 2; }   // stride in x
    int sy() const { return ny + 2; }   // stride in y
    size_t total() const { return static_cast<size_t>(sx()) * sy() * (nz + 2); }
    size_t idx(int i, int j, int k) const {
        return static_cast<size_t>(i + sx() * (j + sy() * k));
    }

    void allocate(int nx_, int ny_, int nz_);
    void initialize(double u0, double v0, double w0, double p0);
};

// ─────────────────────────────────────────────────────
// Spatial discretization (central differences on the collocated grid)
// ─────────────────────────────────────────────────────

namespace spatial {

inline double dudx(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i + 1, j, k)] - g.u[g.idx(i - 1, j, k)]) / (2.0 * g.dx);
}
inline double dudy(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i, j + 1, k)] - g.u[g.idx(i, j - 1, k)]) / (2.0 * g.dy);
}
inline double dudz(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i, j, k + 1)] - g.u[g.idx(i, j, k - 1)]) / (2.0 * g.dz);
}
inline double dvdx(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i + 1, j, k)] - g.v[g.idx(i - 1, j, k)]) / (2.0 * g.dx);
}
inline double dvdy(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i, j + 1, k)] - g.v[g.idx(i, j - 1, k)]) / (2.0 * g.dy);
}
inline double dvdz(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i, j, k + 1)] - g.v[g.idx(i, j, k - 1)]) / (2.0 * g.dz);
}
inline double dwdx(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i + 1, j, k)] - g.w[g.idx(i - 1, j, k)]) / (2.0 * g.dx);
}
inline double dwdy(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i, j + 1, k)] - g.w[g.idx(i, j - 1, k)]) / (2.0 * g.dy);
}
inline double dwdz(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i, j, k + 1)] - g.w[g.idx(i, j, k - 1)]) / (2.0 * g.dz);
}

inline double d2udx2(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i + 1, j, k)] - 2.0 * g.u[g.idx(i, j, k)] + g.u[g.idx(i - 1, j, k)]) /
           (g.dx * g.dx);
}
inline double d2udy2(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i, j + 1, k)] - 2.0 * g.u[g.idx(i, j, k)] + g.u[g.idx(i, j - 1, k)]) /
           (g.dy * g.dy);
}
inline double d2udz2(const FDM3Grid& g, int i, int j, int k) {
    return (g.u[g.idx(i, j, k + 1)] - 2.0 * g.u[g.idx(i, j, k)] + g.u[g.idx(i, j, k - 1)]) /
           (g.dz * g.dz);
}
inline double d2vdx2(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i + 1, j, k)] - 2.0 * g.v[g.idx(i, j, k)] + g.v[g.idx(i - 1, j, k)]) /
           (g.dx * g.dx);
}
inline double d2vdy2(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i, j + 1, k)] - 2.0 * g.v[g.idx(i, j, k)] + g.v[g.idx(i, j - 1, k)]) /
           (g.dy * g.dy);
}
inline double d2vdz2(const FDM3Grid& g, int i, int j, int k) {
    return (g.v[g.idx(i, j, k + 1)] - 2.0 * g.v[g.idx(i, j, k)] + g.v[g.idx(i, j, k - 1)]) /
           (g.dz * g.dz);
}
inline double d2wdx2(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i + 1, j, k)] - 2.0 * g.w[g.idx(i, j, k)] + g.w[g.idx(i - 1, j, k)]) /
           (g.dx * g.dx);
}
inline double d2wdy2(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i, j + 1, k)] - 2.0 * g.w[g.idx(i, j, k)] + g.w[g.idx(i, j - 1, k)]) /
           (g.dy * g.dy);
}
inline double d2wdz2(const FDM3Grid& g, int i, int j, int k) {
    return (g.w[g.idx(i, j, k + 1)] - 2.0 * g.w[g.idx(i, j, k)] + g.w[g.idx(i, j, k - 1)]) /
           (g.dz * g.dz);
}

inline double dpdx(const FDM3Grid& g, int i, int j, int k) {
    return (g.p[g.idx(i + 1, j, k)] - g.p[g.idx(i - 1, j, k)]) / (2.0 * g.dx);
}
inline double dpdy(const FDM3Grid& g, int i, int j, int k) {
    return (g.p[g.idx(i, j + 1, k)] - g.p[g.idx(i, j - 1, k)]) / (2.0 * g.dy);
}
inline double dpdz(const FDM3Grid& g, int i, int j, int k) {
    return (g.p[g.idx(i, j, k + 1)] - g.p[g.idx(i, j, k - 1)]) / (2.0 * g.dz);
}

inline double divergence(const FDM3Grid& g, int i, int j, int k) {
    return dudx(g, i, j, k) + dvdy(g, i, j, k) + dwdz(g, i, j, k);
}

} // namespace spatial

// ─────────────────────────────────────────────────────
// Advection operators
// ─────────────────────────────────────────────────────

namespace advection {

double convective_u(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme);
double convective_v(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme);
double convective_w(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme);

} // namespace advection

// ─────────────────────────────────────────────────────
// Time integration
// ─────────────────────────────────────────────────────

void compute_rhs(const FDM3Grid& g, const FDM3Config& config,
                 std::vector<double>& du_out, std::vector<double>& dv_out,
                 std::vector<double>& dw_out);

void integrate_forward_euler(FDM3Grid& g, const FDM3Config& config, double dt);
void integrate_heun(FDM3Grid& g, const FDM3Config& config, double dt);
void integrate_rk4(FDM3Grid& g, const FDM3Config& config, double dt);
void integrate_crank_nicolson(FDM3Grid& g, const FDM3Config& config, double dt);
void apply_time_integration(FDM3Grid& g, const FDM3Config& config, double dt);

// ─────────────────────────────────────────────────────
// Grid setup / field extraction
// ─────────────────────────────────────────────────────

void initialize_grid(FDM3Grid& g, const FDM3Config& config);
void extract_field(const FDM3Grid& g, const FDM3Config& config, FDM3FieldData& field);

// ─────────────────────────────────────────────────────
// Pressure-velocity coupling
// ─────────────────────────────────────────────────────

void compute_pressure_rhs(FDM3Grid& g, const FDM3Config& config, double dt);
void solve_pressure_poisson(FDM3Grid& g, const FDM3Config& config);
void correct_velocity(FDM3Grid& g, const FDM3Config& config, double dt);

// ─────────────────────────────────────────────────────
// Boundary conditions
// ─────────────────────────────────────────────────────

void apply_boundary_conditions(FDM3Grid& g, const FDM3Config& config);

// Refresh the ghost cells of an arbitrary field at Periodic faces from the
// paired interior values.  Used by the pressure Poisson relaxation so the
// in-place (point-relaxed) field sees a consistent periodic operator.
void update_periodic_field_ghosts(FDM3Grid& g, const FDM3Config& config,
                                  std::vector<double>& field);

// ─────────────────────────────────────────────────────
// Body force helper
// ─────────────────────────────────────────────────────
//
// Adds dt * f to the predictor velocities (f in m/s², force per unit mass).
// The arrays hold one value per interior cell (nx*ny*nz) in the order
// i + nx*(j + ny*k).  Because the force is constant within a step, applying
// it once after the explicit predictor is equivalent to including it in every
// Runge-Kutta stage (for Heun/RK4 the stage weight sums to 1).

void apply_body_force_to_predictor(FDM3Grid& g,
                                   const std::vector<double>& fx,
                                   const std::vector<double>& fy,
                                   const std::vector<double>& fz,
                                   double dt);

// ─────────────────────────────────────────────────────
// Utility / diagnostics
// ─────────────────────────────────────────────────────

double max_velocity(const FDM3Grid& g);
double max_divergence(const FDM3Grid& g);
void compute_diagnostics(const FDM3Grid& g, const FDM3Config& config,
                         FDM3StepResult& step_result);

} // namespace exd::physics::fluid::fdm3