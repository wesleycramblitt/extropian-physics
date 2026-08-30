#pragma once

// ─────────────────────────────────────────────────────
// Engine application: slider-crank + cycle + governor +
// loads over the shared integrator module.
//
//   θ̇ = ω
//   J(θ)·ω̇ = M_gas(θ) − M_load(ω) − ½·(dJ/dθ)·ω²
//
// (the ½·(dJ/dθ)·ω² inertia-torque term is required for
// energy consistency — J varies with the crank angle).
// ─────────────────────────────────────────────────────

#include <exd/physics/control/controller.hpp>
#include <exd/physics/engine/engine_config.hpp>
#include <exd/physics/engine/engine_result.hpp>
#include <exd/physics/model_status.hpp>

namespace exd::physics::engine {

/// Advance the engine by one step (config.integration
/// method; default RK4). `t` is the time at the START of
/// the step; `state` updates in place. `governor` may be
/// null (no governor: full heat release / fixed admission);
/// when non-null it is updated exactly once per step (the
/// throttle is held constant inside the step so RK4 stage
/// evaluations never mutate controller state).
EngineStepResult step_engine(EngineState& state,
                             double t,
                             const EngineConfig& config,
                             control::IController* governor,
                             ModelStatus& status);

/// Time-march the engine over `max_steps`, recording
/// history and (optionally) streaming a CSV machine-state
/// series to `config.csv_path` via the io::CsvSeriesWriter
/// (real-time compatible: one row per step).
EngineSimResult simulate_engine(const EngineConfig& config,
                                ModelStatus& status);

} // namespace exd::physics::engine
