#pragma once

// ─────────────────────────────────────────────────────
// CPU backend (implementation_spec §40–§41, §58).
//
// The backend provides the universal execution primitives
// (parallel_for, reduction, synchronize, allocate, copy).
// Physics modules never name a backend — they call the
// `exec` entry points, and the execution layer routes to
// the active backend (CPU today, CUDA per spec Phase 11).
//
// Determinism: the default CPU backend is serial; a
// thread_count > 1 selects std::execution::par for
// element-wise loops (results remain deterministic for
// independent element kernels).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/memory.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace exd::engine::backends {

class CpuBackend
{
public:
    std::string_view name() const { return "cpu"; }

    /// Number of worker threads for parallel_for (1 = serial, deterministic).
    int thread_count = 1;

    /// Execute fn(i) for i in [0, n). Element order is unspecified when
    /// thread_count > 1 — kernels must be data-independent per element.
    /// (std::execution::par pulls in TBB on libstdc++; std::thread keeps the
    /// engine dependency-free.)
    void parallel_for(size_t n, const std::function<void(size_t)>& fn) const
    {
        if (n == 0) return;
        if (thread_count > 1)
        {
            const size_t nt = static_cast<size_t>(thread_count);
            const size_t block = (n + nt - 1) / nt;
            std::vector<std::thread> workers;
            workers.reserve(nt);
            for (size_t t = 0; t < nt; ++t)
            {
                const size_t lo = t * block;
                const size_t hi = std::min(lo + block, n);
                workers.emplace_back([&, lo, hi]() {
                    for (size_t i = lo; i < hi; ++i) fn(i);
                });
            }
            for (auto& w : workers) w.join();
        }
        else
        {
            for (size_t i = 0; i < n; ++i) fn(i);
        }
    }

    /// Reduce values[0..n) with binary op (deterministic serial order by
    /// default; parallel reduction for associative ops when thread_count > 1).
    double reduce(const double* values, size_t n, double identity,
                  const std::function<double(double, double)>& op) const
    {
        if (thread_count > 1 && n > 1024)
        {
            const size_t nt = static_cast<size_t>(thread_count);
            std::vector<double> partial(nt, identity);
            const size_t block = (n + nt - 1) / nt;
            std::vector<std::thread> workers;
            workers.reserve(nt);
            for (size_t t = 0; t < nt; ++t)
            {
                const size_t lo = t * block;
                const size_t hi = std::min(lo + block, n);
                workers.emplace_back([&, t, lo, hi]() {
                    for (size_t i = lo; i < hi; ++i) partial[t] = op(partial[t], values[i]);
                });
            }
            for (auto& w : workers) w.join();
            double acc = identity;
            for (double p : partial) acc = op(acc, p);
            return acc;
        }
        double acc = identity;
        for (size_t i = 0; i < n; ++i) acc = op(acc, values[i]);
        return acc;
    }

    double sum(const double* values, size_t n) const
    {
        return reduce(values, n, 0.0, [](double a, double b) { return a + b; });
    }

    void synchronize() const {}

    core::Buffer allocate(core::MemorySpace space, size_t bytes) const
    {
        return core::Buffer(space, bytes);
    }
};

/// The process-wide default CPU backend.
inline const CpuBackend& cpu_backend()
{
    static const CpuBackend instance;
    return instance;
}

} // namespace exd::engine::backends
