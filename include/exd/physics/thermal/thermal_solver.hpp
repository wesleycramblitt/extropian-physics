#pragma once

// ---------------------------------------------------------------------
// Phase I thermal domain: steady-state (and transient, implicit CN-style)
// heat equation on a regular structured grid.
//
//   rho*cp dT/dt = k*laplace(T) - rho*cp*(u . grad T) + q_source
//
// Steady solve sets dT/dt = 0; the transient solve uses one implicit
// backward-Euler-ish step per advance (unconditionally stable).  The
// advecting velocity is either uniform (body_velocity) or sampled per
// node from an optional IVectorField3D channel (CHT-lite coupling):
//   u(x) = velocity_channel ? sample(x) : body_velocity.
//
// Discretization is a cell-centered-style FDM reduced to
// the shared node grid exposed to the coupling layer:
// central 7-point Laplacian, first-order upwind advection
// (this solver's stable default; the Peclet guard below
// warns when upwinding is numerically diffusive), SOR
// relaxation (omega = 1.5) to the requested tolerance.
//
// The temperature is returned as a plain
// coupling::StructuredScalarGrid (same origin/spacing/dims
// as the config).  This module does NOT build channels or
// call coupling factories - the caller wraps the grid with
// make_scalar_grid_field() when a channel is needed.
//
// No exceptions, no I/O, pure deterministic free functions.
// ---------------------------------------------------------------------

#include <exd/physics/model_status.hpp>
#include <exd/physics/coupling/field_channels.hpp>   // StructuredScalarGrid (result payload only)

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::thermal {

// ---------------------------------------------------------------------
// Grid configuration
// ---------------------------------------------------------------------

struct ThermalGridConfig {
    std::array<double, 3> origin = {0.0, 0.0, 0.0};     // node (0,0,0) (m)
    std::array<double, 3> spacing = {0.05, 0.05, 0.05}; // node spacing per axis (m), > 0
    std::array<int32_t, 3> dims = {0, 0, 0};            // node counts (>= 2 per axis)
};

// ---------------------------------------------------------------------
// Material configuration
// ---------------------------------------------------------------------

struct ThermalMaterialConfig {
    double conductivity = 50.0;    // W/(m*K)
    double density = 7800.0;       // kg/m^3
    double specific_heat = 500.0;  // J/(kg*K)
};

// ---------------------------------------------------------------------
// Boundary kinds
// ---------------------------------------------------------------------

enum class ThermalBoundaryKind : uint8_t {
    FixedValue,  // prescribed temperature on the face (Dirichlet)
    Insulated,   // zero normal flux, dT/dn = 0 (ghost-cell mirror)
};

// ---------------------------------------------------------------------
// Solver configuration
//
// Face order for boundary_kind / boundary_values:
//   [0] = +x, [1] = -x, [2] = +y, [3] = -y, [4] = +z, [5] = -z
//
// Boundary priority for grid nodes that sit on a corner or
// edge: any FixedValue face pins the node; a node is only
// relaxed (with mirrored ghosts) when every face it touches
// is Insulated.
//
// `dt` is present for interface symmetry with the acoustics
// module and is validated (> 0) but is INERT for the v1
// steady-state solve (an unsteady thermal solve would use
// it in a later phase).
// ---------------------------------------------------------------------

struct ThermalConfig {
    ThermalGridConfig grid;
    ThermalMaterialConfig material;

    std::array<double, 6> boundary_values = {300.0, 300.0, 300.0, 300.0, 300.0, 300.0};
    std::array<ThermalBoundaryKind, 6> boundary_kind = {
        ThermalBoundaryKind::FixedValue,
        ThermalBoundaryKind::FixedValue,
        ThermalBoundaryKind::FixedValue,
        ThermalBoundaryKind::FixedValue,
        ThermalBoundaryKind::FixedValue,
        ThermalBoundaryKind::FixedValue,
    };

    double source_density = 0.0;                 // uniform volumetric heating (W/m^3)
    std::array<double, 3> body_velocity = {0.0, 0.0, 0.0}; // uniform advecting velocity (m/s)
    // Optional per-node advection: when non-null, sampled at every node and
    // overrides body_velocity (out-of-bounds samples fall back to
    // body_velocity with one warning). Caller keeps the channel alive.
    const coupling::IVectorField3D* velocity_channel = nullptr;

    // Transient mode (simulate_thermal): dt is the time step, end_time the
    // horizon, initial_temperature the uniform seed, max_time_steps a cap.
    bool transient = false;
    double end_time = 1.0;             // s, > 0 (transient only)
    double initial_temperature = 300.0; // K, > 0 (transient seed; fixed nodes override)
    uint64_t max_time_steps = 1000000;  // transient step cap, > 0

    double dt = 1e-3;            // transient time step (s); steady solve ignores it
    uint64_t max_steps = 100000; // SOR sweep cap
    double tolerance = 1e-8;     // SOR convergence: max |dT| per sweep <= tolerance
};

// ---------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------

struct ThermalResult {
    bool ok = false;
    ModelStatus status;
    coupling::StructuredScalarGrid temperature;  // node-centered, same grid as config
    double max_residual = 0.0;                   // final SOR residual (max |dT|)
    uint64_t iterations = 0;                     // SOR sweeps executed
    double max_temperature = 0.0;
    double min_temperature = 0.0;
    double total_power = 0.0;                    // sum source_density * cell_volume (W)
};

// ---------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------

/// Validate a thermal config.  Fatal problems return false and set `error`;
/// non-fatal problems (e.g. advection at Pe > 1) are appended to `warnings`.
bool validate_thermal_config(const ThermalConfig& config,
                             std::string& error,
                             std::vector<std::string>& warnings);

/// Solve the steady-state heat equation.
///
/// On failure (`status.ok == false`) the returned result has `ok == false`
/// and an empty temperature grid.  On success, `status` is ok and carries any
/// warnings.  Deterministic given the same config.
ThermalResult solve_thermal(const ThermalConfig& config, ModelStatus& status);

// ---------------------------------------------------------------------
// Transient stepping + coupling-writable state (W11)
// ---------------------------------------------------------------------

/// Mutable thermal state for transient runs and link-driven coupling.
/// `temperature` is node-centered on the config grid; `fixed` marks the
/// nodes pinned by config boundaries or external writes.
struct ThermalState {
    double time = 0.0;                          // s
    coupling::StructuredScalarGrid temperature; // node-centered
    std::vector<bool> fixed;                    // node count, pinned mask
};

/// Build the initial state: uniform seed temperature, config boundaries
/// applied and pinned.  Returns false (with status error) on bad config.
bool init_thermal_state(ThermalState& state, const ThermalConfig& config,
                        ModelStatus& status);

/// Advance the state by one implicit step of `dt` (unconditionally
/// stable; SOR to config.tolerance).  Pinned nodes (config + external
/// writes) hold their values.  Deterministic.
bool advance_thermal(ThermalState& state, double dt, const ThermalConfig& config,
                     ModelStatus& status);

/// Write a temperature at `point`, pinning the nearest node for all
/// subsequent steps.  Used by coupling sinks to impose interface
/// temperatures.  Fails on out-of-bounds points.
bool set_temperature_point(ThermalState& state, const std::array<double, 3>& point,
                           double value, ModelStatus& status);

/// Trilinear scalar channel over a thermal state's temperature field
/// (read-only).  The caller must keep `state` alive.
std::unique_ptr<coupling::IScalarField3D> make_temperature_channel(const ThermalState& state);

/// Time-march the transient solve (transient == false -> the steady solve).
/// Records the final state; on failure result.ok is false.
ThermalResult simulate_thermal(const ThermalConfig& config, ModelStatus& status);

} // namespace exd::physics::thermal