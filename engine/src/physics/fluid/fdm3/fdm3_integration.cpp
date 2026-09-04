#include "fdm3_internal.hpp"
#include <algorithm>

namespace exd::engine::physics::fluid::fdm3 {

// ─────────────────────────────────────────────────────
// Momentum RHS
// ─────────────────────────────────────────────────────
//
// du/dt = -(u·grad)u + nu*Laplacian(u) - grad(p)
//
// Pressure is treated as kinematic pressure (no 1/rho factor), exactly like
// the 2D FDM solver; rho feeds kinematic_viscosity() only.

/// Add the body force (0-based nx·ny·nz accelerations) to a padded RHS.
/// The force is CONSTANT during the step and enters every integrator stage
/// (the old post-step u += dt·f application was a latent time-level
/// inconsistency; the observed Heun instability had a different root —
/// see the HEUN STABILITY FIX note in integrate_heun).
void add_body_force_to_rhs(const FDM3Grid& g,
                           const std::vector<double>* fx,
                           const std::vector<double>* fy,
                           const std::vector<double>* fz,
                           std::vector<double>& du, std::vector<double>& dv,
                           std::vector<double>& dw) {
    if (!fx && !fy && !fz) return;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const size_t id = g.idx(i, j, k);
                const size_t fid = static_cast<size_t>((i - 1) + g.nx * ((j - 1) + g.ny * (k - 1)));
                if (fx) du[id] += (*fx)[fid];
                if (fy) dv[id] += (*fy)[fid];
                if (fz) dw[id] += (*fz)[fid];
            }
}

void compute_rhs(const FDM3Grid& g, const FDM3Config& config,
                 std::vector<double>& du_out, std::vector<double>& dv_out,
                 std::vector<double>& dw_out) {
    const double nu = config.kinematic_viscosity();
    const AdvectionScheme scheme = config.advection_scheme;

    du_out.assign(g.total(), 0.0);
    dv_out.assign(g.total(), 0.0);
    dw_out.assign(g.total(), 0.0);

    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);

                double conv_u = advection::convective_u(g, i, j, k, scheme);
                double diff_u = nu * (spatial::d2udx2(g, i, j, k) +
                                      spatial::d2udy2(g, i, j, k) +
                                      spatial::d2udz2(g, i, j, k));
                du_out[id] = -conv_u + diff_u - spatial::dpdx(g, i, j, k);

                double conv_v = advection::convective_v(g, i, j, k, scheme);
                double diff_v = nu * (spatial::d2vdx2(g, i, j, k) +
                                      spatial::d2vdy2(g, i, j, k) +
                                      spatial::d2vdz2(g, i, j, k));
                dv_out[id] = -conv_v + diff_v - spatial::dpdy(g, i, j, k);

                double conv_w = advection::convective_w(g, i, j, k, scheme);
                double diff_w = nu * (spatial::d2wdx2(g, i, j, k) +
                                      spatial::d2wdy2(g, i, j, k) +
                                      spatial::d2wdz2(g, i, j, k));
                dw_out[id] = -conv_w + diff_w - spatial::dpdz(g, i, j, k);
            }
}

void integrate_forward_euler(FDM3Grid& g, const FDM3Config& config, double dt,
                             const std::vector<double>* fx,
                             const std::vector<double>* fy,
                             const std::vector<double>* fz) {
    std::vector<double> du, dv, dw;
    compute_rhs(g, config, du, dv, dw);
    add_body_force_to_rhs(g, fx, fy, fz, du, dv, dw);

    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] += dt * du[id];
                g.v[id] += dt * dv[id];
                g.w[id] += dt * dw[id];
            }
}

void integrate_heun(FDM3Grid& g, const FDM3Config& config, double dt,
                    const std::vector<double>* fx,
                    const std::vector<double>* fy,
                    const std::vector<double>* fz) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    // k1 (with the constant body force)
    std::vector<double> du1, dv1, dw1;
    compute_rhs(g, config, du1, dv1, dw1);
    add_body_force_to_rhs(g, fx, fy, fz, du1, dv1, dw1);

    // HEUN STABILITY FIX (was: "the Heun body-force bug", forces innocent):
    // the plain trapezoid's second stage samples the FULL-step predictor,
    // and that stage state couples with the stale pressure through the
    // under-relaxed SIMPLE correction into a growing mode in wall+inlet
    // channels at dt*u >= ~0.002 (verified: plain = blowup at step 442,
    // dt*0.5 or Symmetry BCs = stable; RK4's damped tableau never sees it).
    // Two-layer fix (validated: stable 3000+ steps at dt*u = 0.002):
    //   1) the second stage samples the predictor UNDER-RELAXED by the
    //      config's velocity_under_relaxation (the same conservativeness
    //      knob the SIMPLE correction already uses);
    //   2) the predictor is projected divergence-free (inner pressure
    //      correction, u** = u* - dt*grad(p'), Laplacian(p') = div(u*)/dt)
    //      before the k2 evaluation, with the stale real pressure restored
    //      afterwards so k2 sees the same pressure level as k1.  The outer
    //      SIMPLE correction completes the step as before.
    {
        // Predictor: u* = u^n + s*dt*k1  (s = velocity_under_relaxation)
        const double s = config.velocity_under_relaxation;
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    size_t id = g.idx(i, j, k);
                    g.u[id] = g.u_old[id] + s * dt * du1[id];
                    g.v[id] = g.v_old[id] + s * dt * dv1[id];
                    g.w[id] = g.w_old[id] + s * dt * dw1[id];
                }
        std::vector<double> p_save = g.p;
        compute_pressure_rhs(g, config, dt);
        std::fill(g.p.begin(), g.p.end(), 0.0); // zero initial guess + zero-Dirichlet ghosts
        solve_pressure_poisson(g, config);
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    size_t id = g.idx(i, j, k);
                    double dpdx = (g.p[g.idx(i + 1, j, k)] - g.p[g.idx(i - 1, j, k)]) /
                                  (2.0 * g.dx);
                    double dpdy = (g.p[g.idx(i, j + 1, k)] - g.p[g.idx(i, j - 1, k)]) /
                                  (2.0 * g.dy);
                    double dpdz = (g.p[g.idx(i, j, k + 1)] - g.p[g.idx(i, j, k - 1)]) /
                                  (2.0 * g.dz);
                    g.u[id] -= dt * dpdx;
                    g.v[id] -= dt * dpdy;
                    g.w[id] -= dt * dpdz;
                }
        g.p = std::move(p_save);   // restore the stale real pressure for k2
    }

    // k2 (with the constant body force)
    std::vector<double> du2, dv2, dw2;
    compute_rhs(g, config, du2, dv2, dw2);
    add_body_force_to_rhs(g, fx, fy, fz, du2, dv2, dw2);

    // Corrector
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + 0.5 * dt * (du1[id] + du2[id]);
                g.v[id] = g.v_old[id] + 0.5 * dt * (dv1[id] + dv2[id]);
                g.w[id] = g.w_old[id] + 0.5 * dt * (dw1[id] + dw2[id]);
            }
}

void integrate_rk4(FDM3Grid& g, const FDM3Config& config, double dt,
                   const std::vector<double>* fx,
                   const std::vector<double>* fy,
                   const std::vector<double>* fz) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    std::vector<double> k1_u, k1_v, k1_w;
    compute_rhs(g, config, k1_u, k1_v, k1_w);
    add_body_force_to_rhs(g, fx, fy, fz, k1_u, k1_v, k1_w);
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + 0.5 * dt * k1_u[id];
                g.v[id] = g.v_old[id] + 0.5 * dt * k1_v[id];
                g.w[id] = g.w_old[id] + 0.5 * dt * k1_w[id];
            }

    std::vector<double> k2_u, k2_v, k2_w;
    compute_rhs(g, config, k2_u, k2_v, k2_w);
    add_body_force_to_rhs(g, fx, fy, fz, k2_u, k2_v, k2_w);
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + 0.5 * dt * k2_u[id];
                g.v[id] = g.v_old[id] + 0.5 * dt * k2_v[id];
                g.w[id] = g.w_old[id] + 0.5 * dt * k2_w[id];
            }

    std::vector<double> k3_u, k3_v, k3_w;
    compute_rhs(g, config, k3_u, k3_v, k3_w);
    add_body_force_to_rhs(g, fx, fy, fz, k3_u, k3_v, k3_w);
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + dt * k3_u[id];
                g.v[id] = g.v_old[id] + dt * k3_v[id];
                g.w[id] = g.w_old[id] + dt * k3_w[id];
            }

    std::vector<double> k4_u, k4_v, k4_w;
    compute_rhs(g, config, k4_u, k4_v, k4_w);
    add_body_force_to_rhs(g, fx, fy, fz, k4_u, k4_v, k4_w);
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + (dt / 6.0) * (k1_u[id] + 2.0 * k2_u[id] +
                                                      2.0 * k3_u[id] + k4_u[id]);
                g.v[id] = g.v_old[id] + (dt / 6.0) * (k1_v[id] + 2.0 * k2_v[id] +
                                                      2.0 * k3_v[id] + k4_v[id]);
                g.w[id] = g.w_old[id] + (dt / 6.0) * (k1_w[id] + 2.0 * k2_w[id] +
                                                      2.0 * k3_w[id] + k4_w[id]);
            }
}

// NOTE: "Crank-Nicolson" here is implemented as the second-order explicit
// trapezoidal rule (predictor-corrector pair), i.e. Heun.  A true implicit
// Crank-Nicolson requires a linear solve of the coupled system; that is
// future work.  The explicit trapezoid matches the 2D FDM solver's behavior.
void integrate_crank_nicolson(FDM3Grid& g, const FDM3Config& config, double dt,
                              const std::vector<double>* fx,
                              const std::vector<double>* fy,
                              const std::vector<double>* fz) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    // RHS at time n (with the constant body force)
    std::vector<double> du_n, dv_n, dw_n;
    compute_rhs(g, config, du_n, dv_n, dw_n);
    add_body_force_to_rhs(g, fx, fy, fz, du_n, dv_n, dw_n);

    // PROJECTED PREDICTOR (same two-layer fix as integrate_heun): the
    // trapezoid's second stage samples the UNDER-RELAXED predictor,
    // projected divergence-free before the second RHS evaluation.
    {
        // Predictor: u* = u^n + s*dt*du_n  (s = velocity_under_relaxation)
        const double s = config.velocity_under_relaxation;
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    size_t id = g.idx(i, j, k);
                    g.u[id] = g.u_old[id] + s * dt * du_n[id];
                    g.v[id] = g.v_old[id] + s * dt * dv_n[id];
                    g.w[id] = g.w_old[id] + s * dt * dw_n[id];
                }
        std::vector<double> p_save = g.p;
        compute_pressure_rhs(g, config, dt);
        std::fill(g.p.begin(), g.p.end(), 0.0);
        solve_pressure_poisson(g, config);
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    size_t id = g.idx(i, j, k);
                    double dpdx = (g.p[g.idx(i + 1, j, k)] - g.p[g.idx(i - 1, j, k)]) /
                                  (2.0 * g.dx);
                    double dpdy = (g.p[g.idx(i, j + 1, k)] - g.p[g.idx(i, j - 1, k)]) /
                                  (2.0 * g.dy);
                    double dpdz = (g.p[g.idx(i, j, k + 1)] - g.p[g.idx(i, j, k - 1)]) /
                                  (2.0 * g.dz);
                    g.u[id] -= dt * dpdx;
                    g.v[id] -= dt * dpdy;
                    g.w[id] -= dt * dpdz;
                }
        g.p = std::move(p_save);
    }

    // RHS at predicted state (with the constant body force)
    std::vector<double> du_star, dv_star, dw_star;
    compute_rhs(g, config, du_star, dv_star, dw_star);
    add_body_force_to_rhs(g, fx, fy, fz, du_star, dv_star, dw_star);

    // Trapezoidal average
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + 0.5 * dt * (du_n[id] + du_star[id]);
                g.v[id] = g.v_old[id] + 0.5 * dt * (dv_n[id] + dv_star[id]);
                g.w[id] = g.w_old[id] + 0.5 * dt * (dw_n[id] + dw_star[id]);
            }
}

void apply_time_integration(FDM3Grid& g, const FDM3Config& config, double dt,
                            const std::vector<double>* fx,
                            const std::vector<double>* fy,
                            const std::vector<double>* fz) {
    switch (config.time_integration) {
    case TimeIntegration::ForwardEuler: integrate_forward_euler(g, config, dt, fx, fy, fz); break;
    case TimeIntegration::Heun:         integrate_heun(g, config, dt, fx, fy, fz); break;
    case TimeIntegration::RK4:          integrate_rk4(g, config, dt, fx, fy, fz); break;
    case TimeIntegration::CrankNicolson:integrate_crank_nicolson(g, config, dt, fx, fy, fz); break;
    }
}

} // namespace exd::engine::physics::fluid::fdm3