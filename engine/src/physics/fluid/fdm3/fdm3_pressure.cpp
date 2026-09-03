#include "fdm3_internal.hpp"
#include <cmath>

namespace exd::engine::physics::fluid::fdm3 {

// ─────────────────────────────────────────────────────
// Pressure-velocity coupling (SIMPLE)
// ─────────────────────────────────────────────────────

void compute_pressure_rhs(FDM3Grid& g, const FDM3Config&, double dt) {
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i)
                g.rhs[g.idx(i, j, k)] = spatial::divergence(g, i, j, k) / dt;
}

// 7-point SOR point relaxation for the discrete 3D pressure Poisson
// equation.  Second-order centering on the collocated grid:
//
//   (p_E - 2p + p_W)/dx² + (p_N - 2p + p_S)/dy² + (p_F - 2p + p_B)/dz² = rhs
//
// which yields the update
//
//   p_new = (dy²dz²(p_E+p_W) + dx²dz²(p_N+p_S) + dx²dy²(p_F+p_B)
//                               - rhs·dx²dy²dz²) / (2(dy²dz² + dx²dz² + dx²dy²))
//
// NOTE ON SIGN: the 2D fdm solver's relaxation carries +rhs (solving
// Laplacian(p) = -rhs); combined with a u -= dt·grad(p') correction that
// sign never actually projects divergence away.  Here we use the standard
// Poisson sign (Laplacian(p) = +rhs), which together with the
// u = u* - dt·grad(p') correction in correct_velocity() enforces the
// incompressibility constraint.
//
// Ghost cells are treated as Dirichlet values (they are preset by the caller
// before the solve and never updated by the relaxation).
void solve_pressure_poisson(FDM3Grid& g, const FDM3Config& config) {
    const double dx2 = g.dx * g.dx;
    const double dy2 = g.dy * g.dy;
    const double dz2 = g.dz * g.dz;
    const double dy2dz2 = dy2 * dz2;
    const double dx2dz2 = dx2 * dz2;
    const double dx2dy2 = dx2 * dy2;
    const double denom = 2.0 * (dy2dz2 + dx2dz2 + dx2dy2);

    const double omega = config.sor_omega;
    const double tolerance = config.pressure_tolerance;
    const double scale = 1.0 / static_cast<double>(g.nx * g.ny * g.nz);

    double residual = 1.0;
    int iter = 0;
    for (iter = 0; iter < config.pressure_max_iterations && residual > tolerance; ++iter) {
        residual = 0.0;
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    size_t id = g.idx(i, j, k);
                    double p_new = (dy2dz2 * (g.p[g.idx(i + 1, j, k)] + g.p[g.idx(i - 1, j, k)])
                        + dx2dz2 * (g.p[g.idx(i, j + 1, k)] + g.p[g.idx(i, j - 1, k)])
                        + dx2dy2 * (g.p[g.idx(i, j, k + 1)] + g.p[g.idx(i, j, k - 1)])
                        - g.rhs[id] * dx2dy2 * dz2) / denom;
                    double dp = p_new - g.p[id];
                    g.p[id] += omega * dp;
                    residual += dp * dp;
                }
        // Keep periodic ghost values consistent with the evolving field so
        // the periodic direction sees a periodic operator.
        update_periodic_field_ghosts(g, config, g.p);
        residual = std::sqrt(residual * scale);
    }
}

// ─────────────────────────────────────────────────────
// Velocity correction + pressure update
// ─────────────────────────────────────────────────────
//
// Same form as the 2D solver: with u* the predicted velocity and p' the
// pressure-correction (in grid.p_prime),
//
//   u = u_old + alpha_u * (u* - u_old - dt * dp'/dx)   (and v, w)
//   p = p_old + alpha_p * p'
void correct_velocity(FDM3Grid& g, const FDM3Config& config, double dt) {
    const double alpha_u = config.velocity_under_relaxation;
    const double alpha_p = config.pressure_under_relaxation;

    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                double dpdx = (g.p_prime[g.idx(i + 1, j, k)] - g.p_prime[g.idx(i - 1, j, k)]) /
                              (2.0 * g.dx);
                double dpdy = (g.p_prime[g.idx(i, j + 1, k)] - g.p_prime[g.idx(i, j - 1, k)]) /
                              (2.0 * g.dy);
                double dpdz = (g.p_prime[g.idx(i, j, k + 1)] - g.p_prime[g.idx(i, j, k - 1)]) /
                              (2.0 * g.dz);

                g.u[id] = g.u_old[id] + alpha_u * (g.u[id] - g.u_old[id] - dt * dpdx);
                g.v[id] = g.v_old[id] + alpha_u * (g.v[id] - g.v_old[id] - dt * dpdy);
                g.w[id] = g.w_old[id] + alpha_u * (g.w[id] - g.w_old[id] - dt * dpdz);
                g.p[id] = g.p_old[id] + alpha_p * g.p_prime[id];
            }
}

} // namespace exd::engine::physics::fluid::fdm3