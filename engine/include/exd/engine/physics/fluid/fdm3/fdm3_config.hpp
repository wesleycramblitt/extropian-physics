#pragma once

// Configuration types for the 3D incompressible FDM CFD solver (fdm3 module).
//
// The fdm3 module is the 3D counterpart of the 2D FDM solver in
// exd::engine::physics::fluid::fdm.  It reuses the 2D time-integration, advection and
// boundary-condition enums (via using-declarations) so that consumers do not
// need to define or switch on duplicate enum types.

#include <exd/engine/physics/fluid/fdm/fdm_config.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::fdm3 {

// Reused 2D solver enums (see exd::engine::physics::fluid::fdm::fdm_config.hpp).
using exd::engine::physics::fluid::fdm::TimeIntegration;
using exd::engine::physics::fluid::fdm::AdvectionScheme;
using exd::engine::physics::fluid::fdm::FDMBoundaryType;

// Convenience alias for the 2D fdm namespace (used for enum qualification).
namespace fdm = exd::engine::physics::fluid::fdm;

// ─────────────────────────────────────────────────────
// Boundary faces of the 3D box
// ─────────────────────────────────────────────────────

enum class BoundaryFace : uint8_t {
    XMin = 0,   // x = 0 face
    XMax,       // x = lx face
    YMin,       // y = 0 face
    YMax,       // y = ly face
    ZMin,       // z = 0 face
    ZMax,       // z = lz face
};

/// Natural opposite face of a given face (XMin <-> XMax, etc.).
inline BoundaryFace natural_opposite(BoundaryFace f) {
    switch (f) {
    case BoundaryFace::XMin: return BoundaryFace::XMax;
    case BoundaryFace::XMax: return BoundaryFace::XMin;
    case BoundaryFace::YMin: return BoundaryFace::YMax;
    case BoundaryFace::YMax: return BoundaryFace::YMin;
    case BoundaryFace::ZMin: return BoundaryFace::ZMax;
    case BoundaryFace::ZMax: return BoundaryFace::ZMin;
    }
    return BoundaryFace::XMax;
}

// ─────────────────────────────────────────────────────
// Boundary condition definition
// ─────────────────────────────────────────────────────

struct FDM3BoundaryCondition {
    BoundaryFace face = BoundaryFace::XMin;
    FDMBoundaryType type = FDMBoundaryType::Wall;

    // For Inlet: specified velocity components.
    double u_value = 0.0;
    double v_value = 0.0;
    double w_value = 0.0;

    // For FixedPressure: specified pressure.
    double p_value = 0.0;

    // For Periodic: paired face.  Only the natural opposite face is
    // supported in v1; a different pairing emits a validation warning.
    BoundaryFace paired_face = BoundaryFace::XMax;
};

// ─────────────────────────────────────────────────────
// Solver configuration
// ─────────────────────────────────────────────────────

struct FDM3Config {
    // ── Grid ────────────────────────────────────
    int nx = 32;                // Cells in x-direction
    int ny = 32;                // Cells in y-direction
    int nz = 32;                // Cells in z-direction
    double lx = 1.0;            // Domain length in x (m)
    double ly = 1.0;            // Domain length in y (m)
    double lz = 1.0;            // Domain length in z (m)

    // ── Physical parameters ─────────────────────
    double rho = 1.225;         // Density (kg/m³)
    double mu = 1.81e-5;        // Dynamic viscosity (Pa·s)
    double nu = 0.0;            // Kinematic viscosity (m²/s). If 0, mu/rho.

    // ── Time stepping ───────────────────────────
    TimeIntegration time_integration = TimeIntegration::ForwardEuler;
    double dt = 0.001;          // Base time step (s); clamped by CFL when adaptive_dt
    int max_steps = 1000;       // Maximum time steps
    double cfl_target = 0.5;    // Target CFL (for adaptive dt)
    bool adaptive_dt = false;   // Enable CFL-based adaptive dt

    // ── Spatial discretization ──────────────────
    AdvectionScheme advection_scheme = AdvectionScheme::Hybrid;

    // ── Pressure solver ─────────────────────────
    int pressure_max_iterations = 200;
    double pressure_tolerance = 1e-6;
    double sor_omega = 1.2;     // SOR relaxation factor (1.0 = Gauss-Seidel)

    // ── SIMPLE coupling ─────────────────────────
    double velocity_under_relaxation = 0.7;
    double pressure_under_relaxation = 0.3;

    // ── Convergence ─────────────────────────────
    double convergence_tolerance = 1e-6;
    int convergence_window = 100;   // Steps over which to check convergence

    // ── Boundary conditions ─────────────────────
    std::vector<FDM3BoundaryCondition> boundary_conditions;

    // ── Initial conditions ──────────────────────
    double initial_u = 0.0;     // Initial u-velocity (m/s)
    double initial_v = 0.0;     // Initial v-velocity (m/s)
    double initial_w = 0.0;     // Initial w-velocity (m/s)
    double initial_p = 0.0;

    uint32_t field_stamp_interval = 100; // run_fdm3_simulation cadence when no
                                         // OutputScheduler is supplied     // Initial pressure (Pa)

    // ── Derived ─────────────────────────────────
    double dx() const { return lx / static_cast<double>(nx); }
    double dy() const { return ly / static_cast<double>(ny); }
    double dz() const { return lz / static_cast<double>(nz); }
    double kinematic_viscosity() const { return nu > 0.0 ? nu : mu / rho; }

    /// Validate the configuration.  Returns false with a fatal `error` for
    /// invalid settings; non-fatal problems are appended to `warnings`.
    bool validate(std::string& error, std::vector<std::string>& warnings) const;
};

} // namespace exd::engine::physics::fluid::fdm3