#pragma once

#include <exd/engine/core/model_status.hpp>

namespace exd::engine::physics::rigid_body {

/// Backward-compatible alias: the canonical ModelStatus lives at
/// `exd::engine` (umbrella level) so infrastructure modules (solver,
/// coupling) can use it without depending on a domain.
using ModelStatus = exd::engine::core::ModelStatus;

} // namespace exd::engine::physics::rigid_body
