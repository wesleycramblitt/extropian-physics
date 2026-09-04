// ──────────────────────────────────────────────────────────────────────────
// Generic SOR relaxation skeleton (the "one relaxation engine" every Poisson
// domain shares).  What is shared: the sweep loop structure, the omega
// under-relaxation, the residual accounting (RMS or max-change), the
// convergence check, the iteration cap and the diagnostics.  What stays with
// each domain: the point stencil (weights, sign, per-node diagonal), the
// ghost policy (pre-filled padded ghosts, mirror substitution, axis
// neighbor hooks), per-node skipping (Dirichlet/fixed nodes), and any
// between-sweep refresh (e.g. periodic ghosts).
//
// Consumers: fdm3 pressure, fdm 2D pressure, static fields, thermal
// conductivity/advection.  Previously each re-implemented this skeleton —
// and the copies had already drifted (fdm 2D carries +rhs solving
// Lap(p) = -rhs while fdm3/static use the standard sign; the residual
// semantics differ: RMS vs max-change).  Now ONE skeleton with the
// caller-owned point update keeps the sign difference visible at the call
// sites instead of hiding it in duplicated loops.
// ──────────────────────────────────────────────────────────────────────────
#pragma once

#include <cmath>
#include <cstdint>
#include <functional>

namespace exd::engine::numerics {

/// Iteration policy shared by every relaxation consumer.
struct SorConfig {
    double omega = 1.5;
    double tolerance = 1e-6;
    uint64_t max_iterations = 100000;
};

/// How a sweep's residual is measured (callers keep their historical
/// semantics: the RMS forms used by the pressure solvers, the max-change
/// forms used by static fields and thermal).
enum class SorResidualMode {
    RootMeanSquare,   // sqrt( sum(dp^2) / swept_count )
    MaxChange,        // max |dp| over the swept nodes
};

struct SorResult {
    uint64_t iterations = 0;
    bool converged = false;
    double residual = 1.0;
};

/// Run the SOR skeleton over the interior [0..nx) x [0..ny) x [0..nz).
///
/// The caller provides:
///   next(i, j, k)      -> the UNRELAXED new value at the interior index
///                          (the caller resolves its own stencil weights,
///                          sign, ghosts and per-node diagonal)
///   cur(i, j, k)       -> the current value at the interior index
///   set(i, j, k, v)    -> write the relaxed value
///   skip(i, j, k)      -> true = Dirichlet/fixed node: never updated and
///                          excluded from the residual
///   between_sweeps()   -> optional hook run after every sweep (e.g. the
///                          periodic-ghost refresh)
/// The engine applies  v <- cur + omega*(next - cur)  and returns the
/// iteration count, `converged`, and the final residual.  Never throws; a
/// non-converged run is reported via `converged == false` for the caller to
/// warn on.
/// One relaxation pass over the interior; returns the residual per `mode`.
/// The caller owns the convergence loop and iteration accounting (used by
/// the thermal solver, whose point update mixes steady/transient modes).
template <typename NextFn, typename CurFn, typename SetFn, typename SkipFn>
double sor_one_pass(int nx, int ny, int nz,
                    NextFn&& next, CurFn&& cur, SetFn&& set, SkipFn&& skip,
                    double omega, SorResidualMode mode)
{
    double sum_dp2 = 0.0;
    double max_dp = 0.0;
    uint64_t swept = 0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                if (skip(i, j, k))
                    continue;
                const double v_new = next(i, j, k);
                const double v_cur = cur(i, j, k);
                const double dp = v_new - v_cur;
                set(i, j, k, v_cur + omega * dp);
                ++swept;
                sum_dp2 += dp * dp;
                const double adp = std::fabs(dp);
                max_dp = (adp > max_dp) ? adp : max_dp;
            }
    return (mode == SorResidualMode::RootMeanSquare)
        ? ((swept > 0u) ? std::sqrt(sum_dp2 / static_cast<double>(swept)) : 0.0)
        : max_dp;
}

template <typename NextFn, typename CurFn, typename SetFn, typename SkipFn>
SorResult sor_solve(int nx, int ny, int nz,
                    NextFn&& next, CurFn&& cur, SetFn&& set, SkipFn&& skip,
                    const SorConfig& cfg, SorResidualMode mode,
                    const std::function<void()>& between_sweeps = {})
{
    SorResult result;
    for (uint64_t it = 0; it < cfg.max_iterations; ++it) {
        result.residual = sor_one_pass(nx, ny, nz, next, cur, set, skip,
                                       cfg.omega, mode);
        if (between_sweeps)
            between_sweeps();
        result.iterations = it + 1;
        if (result.residual < cfg.tolerance) {
            result.converged = true;
            break;
        }
    }
    return result;
}

} // namespace exd::engine::numerics
