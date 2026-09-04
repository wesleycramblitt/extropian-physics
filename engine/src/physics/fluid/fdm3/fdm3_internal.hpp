#pragma once

// Private header for the 3D FDM CFD solver (fdm3) internals.
// NOT visible to downstream consumers (PRIVATE include path in CMake).

#include <exd/engine/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_result.hpp>

#include <vector>

namespace exd::engine::physics::fluid::fdm3 {

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

// The body-force arrays (nx·ny·nz, 0-based flat) are optional per-step
// accelerations (m/s²).  When non-null they enter the RHS of EVERY stage of
// the time integrator (constant during the step) — the correct treatment
// for all integrators; appending them AFTER the step (the old Heun bug)
// mismatched the time levels and destabilized Heun with sustained forces.
void integrate_forward_euler(FDM3Grid& g, const FDM3Config& config, double dt,
                             const std::vector<double>* fx = nullptr,
                             const std::vector<double>* fy = nullptr,
                             const std::vector<double>* fz = nullptr);
void integrate_heun(FDM3Grid& g, const FDM3Config& config, double dt,
                    const std::vector<double>* fx = nullptr,
                    const std::vector<double>* fy = nullptr,
                    const std::vector<double>* fz = nullptr);
void integrate_rk4(FDM3Grid& g, const FDM3Config& config, double dt,
                   const std::vector<double>* fx = nullptr,
                   const std::vector<double>* fy = nullptr,
                   const std::vector<double>* fz = nullptr);
void integrate_crank_nicolson(FDM3Grid& g, const FDM3Config& config, double dt,
                              const std::vector<double>* fx = nullptr,
                              const std::vector<double>* fy = nullptr,
                              const std::vector<double>* fz = nullptr);
void apply_time_integration(FDM3Grid& g, const FDM3Config& config, double dt,
                            const std::vector<double>* fx = nullptr,
                            const std::vector<double>* fy = nullptr,
                            const std::vector<double>* fz = nullptr);

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

/// How `project_velocity` applies the pressure correction (the ONE
/// fractional-step operator shared by the solver's SIMPLE block and the
/// multi-stage integrators — no more copied Poisson blocks).
enum class ProjectionMode {
    OuterSIMPLE,   // the step's correction: stash p', restore the real
                   // pressure, under-relaxed velocity/pressure update with
                   // u_old/p_old as the reference state
    InnerStage,    // a stage projection: fully remove the predictor's
                   // divergence; the caller's real pressure is restored so
                   // subsequent stages see the same pressure level
};
void project_velocity(FDM3Grid& g, const FDM3Config& config, double dt,
                      ProjectionMode mode);

/// Apply  u -= relax * dt * grad(p_src)  to the velocity field.
void apply_pressure_gradient(FDM3Grid& g, double dt, double relax,
                             const std::vector<double>& p_src);

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
// ─────────────────────────────────────────────────────
// Utility / diagnostics
// ─────────────────────────────────────────────────────

double max_velocity(const FDM3Grid& g);
double max_divergence(const FDM3Grid& g);
void compute_diagnostics(const FDM3Grid& g, const FDM3Config& config,
                         FDM3StepResult& step_result);

} // namespace exd::engine::physics::fluid::fdm3