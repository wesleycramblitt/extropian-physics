#pragma once

// ─────────────────────────────────────────────────────
// Aeroacoustic-lite preset (implementation_spec §52
// "Natural Convection / CFD + wave" coupling family).
//
// Assembly: an fdm3 incompressible duct field (presets/
// cfd/incompressible) supplies the cross-section MEAN
// velocity; the acoustics module solves the linearized
// convected wave equation with that mean flow — the
// documented "acoustics ↔ fdm3" demo.  All numerics stay
// in the physics modules.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/acoustics/wave_solver.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/engine/presets/cfd/incompressible.hpp>

namespace exd::engine::presets::multiphysics {

struct AeroacousticsConfig
{
    // duct/flow (fdm3)
    presets::cfd::IncompressibleCfdConfig flow;
    // acoustics
    double sound_speed = 343.0;
    std::array<double, 3> pulse_center = {0.3, 0.2, 0.2};
    double pulse_width = 0.02;
    double amplitude = 100.0;
    int32_t probe_index = -1;            // flat node; -1 = midpoint of the domain
    uint64_t max_steps = 10000;
    bool include_flow = true;            // false → plain wave (control)
};

struct AeroacousticsResult
{
    bool ok = false;
    core::ModelStatus status;
    exd::engine::physics::fluid::fdm3::FDM3Result flow_result;
    exd::engine::physics::acoustics::WaveResult wave_result;
    std::array<double, 3> mean_velocity = {0, 0, 0};
};

/// Run the aeroacoustic-lite preset.
AeroacousticsResult run_aeroacoustics(const AeroacousticsConfig& config);

} // namespace exd::engine::presets::multiphysics
