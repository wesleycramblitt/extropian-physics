#include "fdm_internal.hpp"
#include <algorithm>

namespace exd::engine::physics::fluid::fdm {

// RHS: du/dt = -(u·∇)u + ν∇²u - (1/ρ)∇p
void compute_rhs(const FDMGrid& g, double nu, AdvectionScheme scheme,
                 std::vector<double>& du_out, std::vector<double>& dv_out) {
    const int s = g.stride();
    du_out.assign(s * (g.ny + 2), 0.0);
    dv_out.assign(s * (g.ny + 2), 0.0);

    for (int j = 1; j <= g.ny; ++j) {
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            double conv_u = advection::convective_u(g, i, j, scheme);
            double diff_u = nu * (spatial::d2udx2(g, i, j) + spatial::d2udy2(g, i, j));
            double grad_px = spatial::dpdx(g, i, j);
            du_out[id] = -conv_u + diff_u - grad_px;

            double conv_v = advection::convective_v(g, i, j, scheme);
            double diff_v = nu * (spatial::d2vdx2(g, i, j) + spatial::d2vdy2(g, i, j));
            double grad_py = spatial::dpdy(g, i, j);
            dv_out[id] = -conv_v + diff_v - grad_py;
        }
    }
}

void integrate_forward_euler(FDMGrid& g, double dt, double nu, AdvectionScheme scheme) {
    std::vector<double> du, dv;
    compute_rhs(g, nu, scheme, du, dv);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] += dt * du[id];
            g.v[id] += dt * dv[id];
        }
}

void integrate_heun(FDMGrid& g, double dt, double nu, AdvectionScheme scheme) {
    const int s = g.stride();
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());

    // k1
    std::vector<double> k1_u, k1_v;
    compute_rhs(g, nu, scheme, k1_u, k1_v);

    // Predictor: u* = u^n + dt*k1
    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + dt * k1_u[id];
            g.v[id] = g.v_old[id] + dt * k1_v[id];
        }

    // k2
    std::vector<double> k2_u, k2_v;
    compute_rhs(g, nu, scheme, k2_u, k2_v);

    // Corrector
    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + 0.5 * dt * (k1_u[id] + k2_u[id]);
            g.v[id] = g.v_old[id] + 0.5 * dt * (k1_v[id] + k2_v[id]);
        }
}

void integrate_rk4(FDMGrid& g, double dt, double nu, AdvectionScheme scheme) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());

    std::vector<double> k1_u, k1_v;
    compute_rhs(g, nu, scheme, k1_u, k1_v);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + 0.5 * dt * k1_u[id];
            g.v[id] = g.v_old[id] + 0.5 * dt * k1_v[id];
        }

    std::vector<double> k2_u, k2_v;
    compute_rhs(g, nu, scheme, k2_u, k2_v);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + 0.5 * dt * k2_u[id];
            g.v[id] = g.v_old[id] + 0.5 * dt * k2_v[id];
        }

    std::vector<double> k3_u, k3_v;
    compute_rhs(g, nu, scheme, k3_u, k3_v);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + dt * k3_u[id];
            g.v[id] = g.v_old[id] + dt * k3_v[id];
        }

    std::vector<double> k4_u, k4_v;
    compute_rhs(g, nu, scheme, k4_u, k4_v);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + (dt / 6.0) * (k1_u[id] + 2*k2_u[id] + 2*k3_u[id] + k4_u[id]);
            g.v[id] = g.v_old[id] + (dt / 6.0) * (k1_v[id] + 2*k2_v[id] + 2*k3_v[id] + k4_v[id]);
        }
}

void integrate_crank_nicolson(FDMGrid& g, double dt, double nu, AdvectionScheme scheme) {
    std::copy(g.u.begin(), g.u.end(), g.u_old.begin());
    std::copy(g.v.begin(), g.v.end(), g.v_old.begin());

    // RHS at time n
    std::vector<double> du_n, dv_n;
    compute_rhs(g, nu, scheme, du_n, dv_n);

    // Predictor
    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + dt * du_n[id];
            g.v[id] = g.v_old[id] + dt * dv_n[id];
        }

    // RHS at predicted state
    std::vector<double> du_star, dv_star;
    compute_rhs(g, nu, scheme, du_star, dv_star);

    // Trapezoidal average
    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            size_t id = g.idx(i, j);
            g.u[id] = g.u_old[id] + 0.5 * dt * (du_n[id] + du_star[id]);
            g.v[id] = g.v_old[id] + 0.5 * dt * (dv_n[id] + dv_star[id]);
        }
}

} // namespace exd::engine::physics::fluid::fdm
