#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::fdm {

// ─────────────────────────────────────────────────────
// Time integration method
// ─────────────────────────────────────────────────────

enum class TimeIntegration : uint8_t {
    ForwardEuler,   // First-order explicit (CFL ≤ 1)
    Heun,           // Second-order explicit Runge-Kutta (CFL ≤ 1)
    RK4,            // Fourth-order explicit Runge-Kutta (CFL ≤ 1)
    CrankNicolson,  // Second-order implicit trapezoidal (unconditionally stable)
};

// ─────────────────────────────────────────────────────
// Advection scheme
// ─────────────────────────────────────────────────────

enum class AdvectionScheme : uint8_t {
    Central,        // Second-order central (stable for Re > ~100)
    Upwind,         // First-order upwind (stable, diffusive)
    Hybrid,         // Central where Pe < 2, upwind otherwise
};

// ─────────────────────────────────────────────────────
// Boundary condition types for FDM
// ─────────────────────────────────────────────────────

enum class FDMBoundaryType : uint8_t {
    Inlet,          // Specified velocity
    Outlet,         // Zero-gradient (convective outflow)
    Wall,           // No-slip (u = 0)
    Symmetry,       // Normal gradient = 0, tangential free
    Periodic,       // Periodic pair
    FixedPressure,  // Specified pressure (Dirichlet)
};

// ─────────────────────────────────────────────────────
// Boundary condition definition
// ─────────────────────────────────────────────────────

enum class BoundaryEdge : uint8_t {
    Left, Right, Bottom, Top
};

struct FDMMaterialBoundaryCondition {
    BoundaryEdge edge = BoundaryEdge::Left;
    FDMBoundaryType type = FDMBoundaryType::Wall;

    // For Inlet: specified velocity components
    double u_value = 0.0;
    double v_value = 0.0;

    // For FixedPressure: specified pressure
    double p_value = 0.0;

    // For Periodic: index of the paired edge (0-3)
    BoundaryEdge paired_edge = BoundaryEdge::Right;
};

// ─────────────────────────────────────────────────────
// Solver configuration
// ─────────────────────────────────────────────────────

struct FDMConfig {
    // ── Grid ──────────────────────────────────────────
    int nx = 64;                // Cells in x-direction
    int ny = 64;                // Cells in y-direction
    double lx = 1.0;            // Domain length in x (m)
    double ly = 1.0;            // Domain length in y (m)

    // ── Physical parameters ───────────────────────────
    double rho = 1.225;         // Density (kg/m³)
    double mu = 1.81e-5;        // Dynamic viscosity (Pa·s)
    double nu = 0.0;            // Kinematic viscosity (m²/s). If 0, computed from mu/rho.

    // ── Time stepping ─────────────────────────────────
    TimeIntegration time_integration = TimeIntegration::ForwardEuler;
    double dt = 0.001;          // Time step (s)
    int max_steps = 1000;       // Maximum time steps
    double cfl_target = 0.5;    // Target CFL (for adaptive dt)
    bool adaptive_dt = false;   // Enable CFL-based adaptive dt

    // ── Spatial discretization ────────────────────────
    AdvectionScheme advection_scheme = AdvectionScheme::Hybrid;

    // ── Pressure solver ───────────────────────────────
    int pressure_max_iterations = 200;
    double pressure_tolerance = 1e-6;
    double sor_omega = 1.2;     // SOR relaxation factor (1.0 = Gauss-Seidel)

    // ── SIMPLE coupling ───────────────────────────────
    double velocity_under_relaxation = 0.7;
    double pressure_under_relaxation = 0.3;
    int simple_max_iterations = 50;
    double simple_tolerance = 1e-5;

    // ── Convergence ───────────────────────────────────
    double convergence_tolerance = 1e-6;
    int convergence_window = 100;   // Steps over which to check convergence

    // ── Output control ────────────────────────────────
    bool output_velocity_field = true;
    bool output_pressure_field = true;
    bool output_vorticity = false;
    bool output_streamfunction = false;
    int output_interval = 100;      // Write results every N steps

    // ── Boundary conditions ───────────────────────────
    std::vector<FDMMaterialBoundaryCondition> boundary_conditions;

    // ── Initial conditions ────────────────────────────
    double initial_u = 0.0;     // Initial u-velocity (m/s)
    double initial_v = 0.0;     // Initial v-velocity (m/s)
    double initial_p = 0.0;     // Initial pressure (Pa)

    // ── Derived (computed in setup) ───────────────────
    double dx() const { return lx / static_cast<double>(nx); }
    double dy() const { return ly / static_cast<double>(ny); }
    double kinematic_viscosity() const { return nu > 0.0 ? nu : mu / rho; }
};

} // namespace exd::engine::physics::fluid::fdm
