#include "fdm_internal.hpp"
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

int solve_pressure_poisson(FDMGrid& g, const std::vector<double>& rhs,
                           int max_iterations, double tolerance, double omega) {
    const int s = g.stride();
    const double dx2 = g.dx * g.dx;
    const double dy2 = g.dy * g.dy;
    const double denom = 2.0 * (dx2 + dy2);

    double residual = 1.0;
    int iter = 0;

    for (iter = 0; iter < max_iterations && residual > tolerance; ++iter) {
        residual = 0.0;
        for (int j = 1; j <= g.ny; ++j) {
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j);
                double p_new = (rhs[id] * dx2 * dy2
                    + dy2 * (g.p[g.idx(i+1, j)] + g.p[g.idx(i-1, j)])
                    + dx2 * (g.p[g.idx(i, j+1)] + g.p[g.idx(i, j-1)])
                ) / denom;
                double dp = p_new - g.p[id];
                g.p[id] += omega * dp;
                residual += dp * dp;
            }
        }
        residual = std::sqrt(residual / (g.nx * g.ny));
    }

    return iter;
}

} // namespace exd::engine::physics::fluid::fdm
