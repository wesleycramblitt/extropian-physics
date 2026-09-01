#pragma once

// ---------------------------------------------------------------------------
// Lagrangian particle tracking over sampled background fields.
//
// A cloud of particles is placed on a deterministic lattice inside a box and
// integrated through
//
//   dx/dt = v
//   dv/dt = g - k*(v - v_f(x))     (v_f from the optional flow channel)
//
// with the shared solver::integrate_step() integrator over one big flat state
// vector [x0,y0,z0,vx0,vy0,vz0, x1, ...].  The flow channel is sampled per
// particle per derivative call; out-of-bounds samples yield v_f = 0 plus a
// single status warning.
//
// Lumped/Lagrangian modules are exempt from the per-domain channel rule
// (same as the engine and compressor modules): the flow channel is an
// optional consumption-only input, not a required output.
// ---------------------------------------------------------------------------

#include <exd/physics/coupling/field_channels.hpp>
#include <exd/physics/model_status.hpp>
#include <exd/physics/solver/integrators.hpp>

#include <array>
#include <cstdint>
#include <vector>

namespace exd::physics::particles {

/// Configuration for a particle-cloud solve.
struct ParticleConfig
{
    uint64_t particle_count = 100;           // requested count, >= 1 (lattice is the closest
                                             // integer decomposition <= requested; the actual
                                             // count is final_positions.size())
    std::array<double, 3> origin = {0.0, 0.0, 0.0};       // spawn box corner (m)
    std::array<double, 3> spawn_extent = {1.0, 1.0, 1.0}; // spawn box size (m), >= 0

    std::array<double, 3> initial_velocity = {0.0, 0.0, 0.0};
    std::array<double, 3> gravity = {0.0, 0.0, -9.81};    // m/s^2

    double drag_coefficient = 0.0;             // specific drag k (1/s), >= 0
    const coupling::IVectorField3D* flow_channel = nullptr; // optional background v_f (m/s)

    double dt = 1e-3;                          // time step (s), > 0
    uint64_t max_steps = 10000;                // step limit, > 0
    uint64_t history_interval = 10;            // probe cadence (steps), >= 1

    solver::IntegratorConfig integration;      // default RK4
};

/// Validate the particle config.  Fatal problems return false and fill
/// `error`; non-fatal observations are appended to `warnings`.
bool validate_particle_config(const ParticleConfig& config,
                              std::string& error,
                              std::vector<std::string>& warnings);

/// Result of a particle solve.
struct ParticleResult
{
    bool ok = false;
    ModelStatus status;

    std::vector<std::array<double, 3>> final_positions;  // one entry per spawned particle
    std::vector<std::array<double, 3>> final_velocities;

    std::vector<std::array<double, 3>> trajectory_probe; // position history of particle 0
    std::vector<std::array<double, 3>> velocity_probe;   // velocity history of particle 0
    std::vector<double> time_history;                    // times matching both probes

    double mean_speed = 0.0; // mean |v| over final_velocities (m/s)
};

/// Integrate the cloud.  Deterministic: identical configs give bit-identical
/// results.  Failures surface through `status` and the result.ok flag.
ParticleResult solve_particles(const ParticleConfig& config, ModelStatus& status);

} // namespace exd::physics::particles