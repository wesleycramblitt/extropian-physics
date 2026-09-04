#include "fdm3_internal.hpp"

#include <exd/engine/numerics/sor.hpp>

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

    // The point stencil (weights, standard Poisson sign, padded-Dirichlet
    // ghosts) stays here; the sweep/relaxation/residual skeleton is shared
    // (numerics::sor_solve).
    // NOTE: the engine iterates 0-based indices over [0..nx) while the
    // padded grid's interior sits at [1..nx] — the wrappers must shift.
    auto next = [&](int i, int j, int k) {
        const int I = i + 1, J = j + 1, K = k + 1;
        const size_t id = g.idx(I, J, K);
        return (dy2dz2 * (g.p[g.idx(I + 1, J, K)] + g.p[g.idx(I - 1, J, K)])
                + dx2dz2 * (g.p[g.idx(I, J + 1, K)] + g.p[g.idx(I, J - 1, K)])
                + dx2dy2 * (g.p[g.idx(I, J, K + 1)] + g.p[g.idx(I, J, K - 1)])
                - g.rhs[id] * dx2dy2 * dz2) / denom;
    };
    auto cur = [&](int i, int j, int k) { return g.p[g.idx(i + 1, j + 1, k + 1)]; };
    auto setp = [&](int i, int j, int k, double v) { g.p[g.idx(i + 1, j + 1, k + 1)] = v; };
    auto skip_none = [](int, int, int) { return false; };

    numerics::SorConfig cfg;
    cfg.omega = config.sor_omega;
    cfg.tolerance = config.pressure_tolerance;
    cfg.max_iterations = config.pressure_max_iterations;

    // Keep periodic ghost values consistent with the evolving field so the
    // periodic direction sees a periodic operator (no-op without periodic
    // faces — the refresh iterates the face list).
    auto refresh = [&]() { update_periodic_field_ghosts(g, config, g.p); };
    auto r = numerics::sor_solve(g.nx, g.ny, g.nz, next, cur, setp, skip_none,
                                 cfg, numerics::SorResidualMode::RootMeanSquare, refresh);
    (void)r;   // the caller policy treats a capped run as acceptable
}

// ─────────────────────────────────────────────────────
// Shared fractional-step projection operator
// ─────────────────────────────────────────────────────
void apply_pressure_gradient(FDM3Grid& g, double dt, double relax,
                             const std::vector<double>& p_src) {
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const size_t id = g.idx(i, j, k);
                const double dpdx = (p_src[g.idx(i + 1, j, k)] - p_src[g.idx(i - 1, j, k)]) /
                                    (2.0 * g.dx);
                const double dpdy = (p_src[g.idx(i, j + 1, k)] - p_src[g.idx(i, j - 1, k)]) /
                                    (2.0 * g.dy);
                const double dpdz = (p_src[g.idx(i, j, k + 1)] - p_src[g.idx(i, j, k - 1)]) /
                                    (2.0 * g.dz);
                g.u[id] -= relax * dt * dpdx;
                g.v[id] -= relax * dt * dpdy;
                g.w[id] -= relax * dt * dpdz;
            }
}

void project_velocity(FDM3Grid& g, const FDM3Config& config, double dt,
                      ProjectionMode mode) {
    if (mode == ProjectionMode::InnerStage) {
        // The stage projection: solve p' from div(u*)/dt with a zero initial
        // guess (zero-Dirichlet ghosts, same policy as the outer solve),
        // fully correct the stage state, and restore the caller's real
        // pressure so the next stage's RHS sees the same pressure as k1.
        const std::vector<double> p_save = g.p;
        compute_pressure_rhs(g, config, dt);
        std::fill(g.p.begin(), g.p.end(), 0.0);
        solve_pressure_poisson(g, config);
        apply_pressure_gradient(g, dt, 1.0, g.p);
        g.p = p_save;
        return;
    }
    // OuterSIMPLE: the step's correction — p' stays in the stash, the real
    // pressure is restored, and the under-relaxed correction is applied
    // against the pre-step reference state.
    compute_pressure_rhs(g, config, dt);
    std::fill(g.p.begin(), g.p.end(), 0.0);
    solve_pressure_poisson(g, config);
    std::copy(g.p.begin(), g.p.end(), g.p_prime.begin());
    std::copy(g.p_old.begin(), g.p_old.end(), g.p.begin());
    correct_velocity(g, config, dt);
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

    // u'' = u_old + alpha_u*(u - u_old) - alpha_u*dt*grad(p')
    //      = blend (keeping 1-alpha of the OLD state) THEN the shared
    //        gradient at relax = alpha_u
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const size_t id = g.idx(i, j, k);
                g.u[id] = (1.0 - alpha_u) * g.u_old[id] + alpha_u * g.u[id];
                g.v[id] = (1.0 - alpha_u) * g.v_old[id] + alpha_u * g.v[id];
                g.w[id] = (1.0 - alpha_u) * g.w_old[id] + alpha_u * g.w[id];
                g.p[id] = g.p_old[id] + alpha_p * g.p_prime[id];
            }
    apply_pressure_gradient(g, dt, alpha_u, g.p_prime);
}

} // namespace exd::engine::physics::fluid::fdm3