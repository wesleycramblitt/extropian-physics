// Benchmark case registry (W18 plan, Sections 2 / 2A / 2B).
#pragma once

#include "benchmark_common.hpp"

#include <string>
#include <vector>

namespace bench {

struct CaseSpec {
    const char* name;        // e.g. "mms"
    const char* family;      // e.g. "B1 MMS"
    const char* anchor;      // e.g. "exact manufactured solution"
    void (*run)(const RunSpec&);
};

std::vector<CaseSpec> all_cases();

/// Module-family registry (Section 2B C-series).
std::vector<CaseSpec> module_cases();

} // namespace bench
