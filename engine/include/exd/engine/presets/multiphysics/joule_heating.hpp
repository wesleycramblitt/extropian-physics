#pragma once

// ─────────────────────────────────────────────────────
// Joule heating preset (implementation_spec §52 "Joule
// Heating" EM preset).
//
// Assembly: the electromagnetics static-field module solves
// the potential between two electrodes; the volumetric heat
// q = σ·|E|² feeds the thermal module through the SOURCE
// CHANNEL (spatially varying q).  Power anchor:
// P = ∫σ|E|²dV ≈ V²·σ·A/L (parallel-plate); the steady
// thermal field with fixed ends is the parabola
// T(x) = T0 + (q̄/2k)·x·(L−x).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/electromagnetics/static_fields.hpp>
#include <exd/engine/physics/thermal/thermal_solver.hpp>

#include <memory>

namespace exd::engine::presets::multiphysics {

struct JouleHeatingConfig
{
    // electrodes
    int nx = 41, ny = 9, nz = 9;
    double spacing = 0.025;          // m
    double voltage = 10.0;           // V across the x faces
    // conductor
    double conductivity = 5.96e7;    // copper S/m (rate the heat source)
    // thermal
    double thermal_conductivity = 400.0;   // W/(m·K)
    double t_wall = 300.0;           // fixed temperature at both x faces
    double tolerance = 1e-6;
    uint64_t max_iterations = 200000;
};

struct JouleHeatingResult
{
    bool ok = false;
    core::ModelStatus status;
    exd::engine::physics::electromagnetics::StaticFieldResult electric;
    exd::engine::physics::thermal::ThermalResult thermal;
    double total_power = 0.0;        // ∫σ|E|²dV (W)
    double analytic_power = 0.0;     // V²σA/d (W)
    double center_source = 0.0;      // q at the domain center (W/m³)
    double mid_temperature = 0.0;
};

/// Run the Joule heating preset.
JouleHeatingResult run_joule_heating(const JouleHeatingConfig& config);

} // namespace exd::engine::presets::multiphysics
