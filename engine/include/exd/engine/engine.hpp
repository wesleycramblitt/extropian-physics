#pragma once

// ─────────────────────────────────────────────────────
// extropian-physics engine — umbrella header.
//
// One composable numerical runtime (docs/implementation_spec.md):
//   DATA + OPERATORS + PHYSICS + DISCRETIZATION + COUPLING
//   + EXECUTION + HARDWARE
//
// Layer order (spec §56): core → numerics → discretization →
// physics → coupling → presets; backends are orthogonal.
//
// Include the umbrella for the core runtime pieces; modules
// are included on demand (include what you use).
// ─────────────────────────────────────────────────────

// core runtime (spec Phase 1)
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/core/units.hpp>
#include <exd/engine/core/memory.hpp>
#include <exd/engine/core/field.hpp>
#include <exd/engine/core/entity_set.hpp>
#include <exd/engine/core/state.hpp>
#include <exd/engine/core/operator.hpp>
#include <exd/engine/core/execution.hpp>

// backends (orthogonal; CPU today, CUDA per spec Phase 11)
#include <exd/engine/backends/cpu.hpp>

// mesh (spec Phase 2)
#include <exd/engine/mesh/structured.hpp>
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/mesh/boundary.hpp>
#include <exd/engine/mesh/validation.hpp>

// numerics (spec Phase 5): ODE integrators + time stepping
#include <exd/engine/numerics/integrators.hpp>
#include <exd/engine/numerics/time_stepping.hpp>

// coupling (spec Phase 9)
#include <exd/engine/coupling/plugin_interface.hpp>
#include <exd/engine/coupling/field_channels.hpp>
