#include "fdm_internal.hpp"

#include <exd/engine/numerics/sor.hpp>

#include <cmath>
#include <algorithm>

namespace exd::engine::physics::fluid::fdm {

void compute_pressure_rhs(const FDMGrid& g, double dt, std::vector<double>& rhs) {
    const int s = g.stride();
    rhs.assign(s * (g.ny + 2), 0.0);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            rhs[g.idx(i, j)] = spatial::divergence(g, i, j) / dt;
        }
}

// 2D point stencil (weights + the historical +rhs convention solving
// Lap(p) = -rhs — kept explicit HERE rather than dormant in a duplicated
// loop; the shared skeleton is numerics::sor_solve, nz = 1 gives the
// 5-point form).
int solve_pressure_poisson(FDMGrid& g, const std::vector<double>& rhs,
                           int max_iterations, double tolerance, double omega) {
    const int s = g.stride();
    const double dx2 = g.dx * g.dx;
    const double dy2 = g.dy * g.dy;
    const double denom = 2.0 * (dx2 + dy2);

    // the engine iterates 0-based [0..nx); the padded interior is [1..nx]
    auto next = [&](int i, int j, int) {
        const int I = i + 1, J = j + 1;
        const size_t id = g.idx(I, J);
        return (rhs[id] * dx2 * dy2
                + dy2 * (g.p[g.idx(I + 1, J)] + g.p[g.idx(I - 1, J)])
                + dx2 * (g.p[g.idx(I, J + 1)] + g.p[g.idx(I, J - 1)])) / denom;
    };
    auto cur = [&](int i, int j, int) { return g.p[g.idx(i + 1, j + 1)]; };
    auto setp = [&](int i, int j, int, double v) { g.p[g.idx(i + 1, j + 1)] = v; };
    auto skip_none = [](int, int, int) { return false; };

    numerics::SorConfig cfg;
    cfg.omega = omega;
    cfg.tolerance = tolerance;
    cfg.max_iterations = static_cast<uint64_t>(max_iterations);
    auto r = numerics::sor_solve(g.nx, g.ny, 1, next, cur, setp, skip_none,
                                 cfg, numerics::SorResidualMode::RootMeanSquare);
    return static_cast<int>(r.iterations);
}

} // namespace exd::engine::physics::fluid::fdm
