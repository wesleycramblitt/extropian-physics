#pragma once

// ─────────────────────────────────────────────────────
// Species-in-flow preset (implementation_spec §52 CFD +
// §66 reaction-diffusion; combustion-lite).
//
// Assembly: an fdm3 channel flow supplies the velocity
// channel; the species module transports a decaying scalar
// (first-order kinetics).  Steady anchor:
// c_out/c_in = exp(−k·L/u_mean) (advection-dominated).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/engine/physics/species/species_solver.hpp>
#include <exd/engine/presets/cfd/incompressible.hpp>

#include <memory>

namespace exd::engine::presets::multiphysics {

struct SpeciesInFlowConfig
{
    presets::cfd::IncompressibleCfdConfig flow;   // duct + velocity
    double decay_rate = 0.1;                      // 1/s (first-order)
    double inlet_concentration = 1.0;
    double species_dt = 0.01;
    uint64_t max_steps = 400000;
    double steady_tolerance = 1e-11;
};

struct SpeciesInFlowResult
{
    bool ok = false;
    core::ModelStatus status;
    exd::engine::physics::fluid::fdm3::FDM3Result flow_result;
    exd::engine::physics::species::SpeciesResult species;
    double outlet_concentration = 0.0;
    double analytic_outlet = 0.0;                 // c_in·exp(−k·L/u)
};

/// Run the species-in-flow preset.
SpeciesInFlowResult run_species_in_flow(const SpeciesInFlowConfig& config);

} // namespace exd::engine::presets::multiphysics
