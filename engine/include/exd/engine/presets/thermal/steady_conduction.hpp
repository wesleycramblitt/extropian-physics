#pragma once

// ─────────────────────────────────────────────────────
// Steady heat conduction preset (implementation_spec
// §52 "Heat Conduction" CLASS — NOT a product).
//
// Assembly: the heterogeneous-conduction module with
// sensible defaults (all faces at the sink temperature;
// region builders for materials/sources).  APPLICABLE TO
// AN ENTIRE CLASS: chiplet boards, heat sinks, cooling
// channels, reactor walls, building envelopes — the user
// imports geometry/materials/BCs/parameters into the
// config and the same solver runs:
//
//   config.grid       ← CAD/mesh geometry
//   region_materials  ← CAD bodies (conductivity per body)
//   region_sources    ← loads (W/m³ per volume)
//   face kind/value   ← BCs
//   per-node fields   ← data-driven override (any CAD
//                       tool can emit them)
//
// The module owns the numerics; the preset owns the
// defaults and the entry point.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/thermal/heterogeneous.hpp>

namespace exd::engine::presets::thermal {

using SteadyConductionConfig = exd::engine::physics::thermal::HeterogeneousConductionConfig;
using SteadyConductionResult = exd::engine::physics::thermal::HeterogeneousConductionResult;

/// Run the steady-conduction class with the given (user-supplied) data.
inline SteadyConductionResult run_steady_conduction(const SteadyConductionConfig& config)
{
    return exd::engine::physics::thermal::solve_heterogeneous_conduction(config);
}

} // namespace exd::engine::presets::thermal
