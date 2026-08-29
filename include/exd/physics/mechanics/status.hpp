#pragma once

#include <exd/physics/model_status.hpp>

namespace exd::physics::mechanics {

/// Backward-compatible alias: the canonical ModelStatus lives at
/// `exd::physics` (umbrella level) so infrastructure modules (solver,
/// coupling) can use it without depending on a domain.
using ModelStatus = exd::physics::ModelStatus;

} // namespace exd::physics::mechanics
