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

void integrate_forward_euler(FDM3Grid& g, const FDM3Config& config, double dt) {
    std::vector<double> du, dv, dw;
    compute_rhs(g, config, du, dv, dw);

    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] += dt * du[id];
                g.v[id] += dt * dv[id];
                g.w[id] += dt * dw[id];
            }
}

void integrate_heun(FDM3Grid& g, const FDM3Config& config, double dt) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    // k1
    std::vector<double> du1, dv1, dw1;
    compute_rhs(g, config, du1, dv1, dw1);

    // Predictor: u* = u^n + dt*k1
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + dt * du1[id];
                g.v[id] = g.v_old[id] + dt * dv1[id];
                g.w[id] = g.w_old[id] + dt * dw1[id];
            }

    // k2
    std::vector<double> du2, dv2, dw2;
    compute_rhs(g, config, du2, dv2, dw2);

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

void integrate_rk4(FDM3Grid& g, const FDM3Config& config, double dt) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    std::vector<double> k1_u, k1_v, k1_w;
    compute_rhs(g, config, k1_u, k1_v, k1_w);
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
void integrate_crank_nicolson(FDM3Grid& g, const FDM3Config& config, double dt) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());
    std::copy(g.w.begin(), g.w.end(), g.w_old.begin());

    // RHS at time n
    std::vector<double> du_n, dv_n, dw_n;
    compute_rhs(g, config, du_n, dv_n, dw_n);

    // Predictor (explicit Euler step)
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                g.u[id] = g.u_old[id] + dt * du_n[id];
                g.v[id] = g.v_old[id] + dt * dv_n[id];
                g.w[id] = g.w_old[id] + dt * dw_n[id];
            }

    // RHS at predicted state
    std::vector<double> du_star, dv_star, dw_star;
    compute_rhs(g, config, du_star, dv_star, dw_star);

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

void apply_time_integration(FDM3Grid& g, const FDM3Config& config, double dt) {
    switch (config.time_integration) {
    case TimeIntegration::ForwardEuler: integrate_forward_euler(g, config, dt); break;
    case TimeIntegration::Heun:         integrate_heun(g, config, dt); break;
    case TimeIntegration::RK4:          integrate_rk4(g, config, dt); break;
    case TimeIntegration::CrankNicolson:integrate_crank_nicolson(g, config, dt); break;
    }
}

void apply_body_force_to_predictor(FDM3Grid& g,
                                   const std::vector<double>& fx,
                                   const std::vector<double>& fy,
                                   const std::vector<double>& fz,
                                   double dt) {
    const size_t cells = static_cast<size_t>(g.nx) * g.ny * g.nz;
    if (fx.size() < cells || fy.size() < cells || fz.size() < cells)
        return;

    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                size_t fid = static_cast<size_t>((i - 1) + g.nx * ((j - 1) + g.ny * (k - 1)));
                g.u[id] += dt * fx[fid];
                g.v[id] += dt * fy[fid];
                g.w[id] += dt * fz[fid];
            }
}

} // namespace exd::engine::physics::fluid::fdm3