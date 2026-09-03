#pragma once

#include <string>
#include <vector>

namespace exd::engine::core {

/// Status channel shared by all solver modules (no-exception doctrine).
/// `ok == false` is a fatal error; `warnings` are non-fatal.
struct ModelStatus
{
    bool ok = true;
    std::string error;
    std::vector<std::string> warnings;
};

} // namespace exd::engine::core

// The engine root re-exports ModelStatus so every layer (`exd::engine::*`)
// can use the bare name exactly as `exd::physics::ModelStatus` was used
// under the old root namespace (spec: single no-exception status channel).
namespace exd::engine {
using core::ModelStatus;
} // namespace exd::engine
