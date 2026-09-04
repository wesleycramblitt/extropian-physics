// ────────────────────────────────────────────────────────────────────────
// Shared harness helpers for the benchmark demos (W18 plan; see
// docs/benchmark_plan.md for the family definitions and tiers).
//
// Every case: configure → solve → compare against the exact/public anchor →
// print a metric row (name, family, anchor, value, reference, tier state,
// wall time).  No assertions: the demos are meant to be READ, run in CI
// smoke mode (tiny grids, seconds) and re-run at full size on a capable
// machine via `--full`.
// ────────────────────────────────────────────────────────────────────────
#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace bench {

/// Wall-clock seconds since program start (steady).
struct Stopwatch {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    double elapsed_s() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
};

/// Relative L2 error of `got` vs `ref` (same length; 0 for an all-zero ref).
inline double rel_l2(const std::vector<double>& got, const std::vector<double>& ref) {
    double num = 0.0, den = 0.0;
    const size_t n = std::min(got.size(), ref.size());
    for (size_t i = 0; i < n; ++i) {
        const double d = got[i] - ref[i];
        num += d * d;
        den += ref[i] * ref[i];
    }
    return (den > 0.0) ? std::sqrt(num / den) : 0.0;
}

inline double rel_l2(const double* got, const double* ref, size_t n) {
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double d = got[i] - ref[i];
        num += d * d;
        den += ref[i] * ref[i];
    }
    return (den > 0.0) ? std::sqrt(num / den) : 0.0;
}

/// One result row (CSV-style single line for later phases to collect).
inline void row(const std::string& family, const std::string& metric,
                double value, double reference, const std::string& tier,
                double seconds) {
    std::cout << family << "," << metric << "," << value << "," << reference << ","
              << tier << "," << seconds << "\n";
}

/// Human-readable verdict line with the tier label.
inline void verdict(const std::string& family, const std::string& anchor,
                    const std::string& tier, const std::string& detail) {
    std::cout << "[bench] " << family << " — " << anchor << " [" << tier << "] "
              << detail << "\n";
}

/// Configuration of the current run (grid size + tier) so each case can
/// size its workload: smoke = seconds, full = the planned sweeps.
struct RunSpec {
    int grid = 0;          // case-specific refinement knob (0 = case default)
    bool full = false;     // full sweep vs smoke tier
};

} // namespace bench
