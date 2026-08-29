#pragma once

#include <string>
#include <vector>

namespace exd::physics {

/// Status channel shared by all solver modules (no-exception doctrine).
/// `ok == false` is a fatal error; `warnings` are non-fatal.
struct ModelStatus
{
    bool ok = true;
    std::string error;
    std::vector<std::string> warnings;
};

} // namespace exd::physics
