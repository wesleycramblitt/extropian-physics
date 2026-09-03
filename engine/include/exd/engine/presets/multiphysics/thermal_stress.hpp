#pragma once

// ─────────────────────────────────────────────────────
// Thermal-stress preset (implementation_spec §52 "Thermal
// Stress" FEA preset, FDM-flavored).
//
// Assembly: the thermal module produces a temperature
// field; the structural module solves the Navier–Cauchy
// displacement with the thermal strain ε_th = α·(T−T_ref)
// driven by the temperature CHANNEL — the existing
// coupling point of the two modules.  Free-bar case:
// u_x(x) = α·g·x²/2 for T(x) = T_ref + g·x (exact).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/structural/elasticity.hpp>
#include <exd/engine/physics/thermal/thermal_solver.hpp>

#include <memory>

namespace exd::engine::presets::multiphysics {

struct ThermalStressConfig
{
    // thermal slab: linear T from left to right
    int nx = 31, ny = 4, nz = 4;
    double length = 1.0;
    double t_left = 300.0;
    double t_right = 400.0;
    double conductivity = 50.0;
    // structure (same lattice)
    double elastic_modulus = 200e9;
    double poisson_ratio = 0.3;
    double thermal_expansion = 1e-5;
    double t_ref = 300.0;
    double tolerance = 1e-10;
    uint64_t max_iterations = 200000;
};

struct ThermalStressResult
{
    bool ok = false;
    core::ModelStatus status;
    exd::engine::physics::thermal::ThermalResult thermal;
    exd::engine::physics::structural::ElasticityResult structural;
    double measured_tip_displacement = 0.0;   // u_x at the bar end (center line)
};

/// Run the thermal-stress preset.
ThermalStressResult run_thermal_stress(const ThermalStressConfig& config);

} // namespace exd::engine::presets::multiphysics
