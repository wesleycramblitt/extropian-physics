#pragma once

// ─────────────────────────────────────────────────────
// Incompressible CFD preset (implementation_spec §14,
// §52 "Incompressible CFD").
//
// Assembly, not a solver: builds a default fdm3
// incompressible-flow configuration (channel: inlet →
// no-slip walls → fixed-pressure outlet) with default
// numerical settings and boundary conventions.  All
// numerical work stays in physics/fluid/fdm3.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_solver.hpp>

namespace exd::engine::presets::cfd {

struct IncompressibleCfdConfig
{
    // domain
    double length = 1.0, width = 0.4, height = 0.4;
    int nx = 20, ny = 8, nz = 8;
    // flow
    double inlet_velocity = 0.2;
    double rho = 1.225;
    double mu = 0.05;              // ν·dt·Σ1/h² < 0.5 at dt = 0.01
    double dt = 0.01;
    int max_steps = 900;
    // numerics defaults (fdm3 SIMPLE + SOR)
    double pressure_tolerance = 5e-6;
    int pressure_max_iterations = 60;
    double sor_omega = 1.45;
    double velocity_under_relaxation = 0.8;
    double pressure_under_relaxation = 0.35;
    double convergence_tolerance = 1e-6;
};

/// Build the default channel-flow fdm3 configuration (inlet ×4 walls ×
/// fixed-pressure outlet).  The caller can override any field afterwards —
/// presets provide defaults, not hardcoded modes (§45).
inline physics::fluid::fdm3::FDM3Config make_channel_flow(
    const IncompressibleCfdConfig& cfg, core::ModelStatus& status)
{
    using namespace exd::engine::physics::fluid::fdm3;
    FDM3Config flow;
    flow.nx = cfg.nx; flow.ny = cfg.ny; flow.nz = cfg.nz;
    flow.lx = cfg.length; flow.ly = cfg.width; flow.lz = cfg.height;
    flow.rho = cfg.rho;
    flow.mu = cfg.mu;
    flow.dt = cfg.dt;
    flow.max_steps = cfg.max_steps;
    flow.pressure_max_iterations = cfg.pressure_max_iterations;
    flow.pressure_tolerance = cfg.pressure_tolerance;
    flow.sor_omega = cfg.sor_omega;
    flow.velocity_under_relaxation = cfg.velocity_under_relaxation;
    flow.pressure_under_relaxation = cfg.pressure_under_relaxation;
    flow.convergence_tolerance = cfg.convergence_tolerance;

    FDM3BoundaryCondition inlet;
    inlet.face = BoundaryFace::XMin;
    inlet.type = FDMBoundaryType::Inlet;
    inlet.u_value = cfg.inlet_velocity;
    flow.boundary_conditions.push_back(inlet);
    FDM3BoundaryCondition outlet;
    outlet.face = BoundaryFace::XMax;
    outlet.type = FDMBoundaryType::FixedPressure;
    outlet.p_value = 0.0;
    flow.boundary_conditions.push_back(outlet);
    for (auto face : {BoundaryFace::YMin, BoundaryFace::YMax,
                      BoundaryFace::ZMin, BoundaryFace::ZMax})
    {
        FDM3BoundaryCondition bc;
        bc.face = face;
        bc.type = FDMBoundaryType::Wall;
        flow.boundary_conditions.push_back(bc);
    }
    std::string error;
    std::vector<std::string> warnings;
    if (!flow.validate(error, warnings))
    {
        status.ok = false;
        status.error = "preset incompressible_cfd: " + error;
        for (auto& w : warnings) status.warnings.push_back(w);
        return flow;
    }
    return flow;
}

} // namespace exd::engine::presets::cfd
