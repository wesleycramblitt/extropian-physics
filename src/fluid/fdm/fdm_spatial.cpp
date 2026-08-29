#include "fdm_internal.hpp"
#include <cmath>
#include <algorithm>

namespace exd::physics::fluid::fdm {
namespace advection {

// Upwind interpolation helper
static double upwind_interp(double vm, double vc, double vp, double vel) {
    if (vel > 0.0) return vc;
    return vp;
}

double convective_u(const FDMGrid& g, int i, int j, AdvectionScheme scheme) {
    double u = g.u[g.idx(i, j)];
    double v_at_u = 0.25 * (g.v[g.idx(i-1, j)] + g.v[g.idx(i, j)] +
                             g.v[g.idx(i-1, j+1)] + g.v[g.idx(i, j+1)]);

    double dudx_adv, dudy_adv;

    switch (scheme) {
    case AdvectionScheme::Central:
        dudx_adv = u * (g.u[g.idx(i+1, j)] - g.u[g.idx(i-1, j)]) / (2.0 * g.dx);
        dudy_adv = v_at_u * (g.u[g.idx(i, j+1)] - g.u[g.idx(i, j-1)]) / (2.0 * g.dy);
        break;

    case AdvectionScheme::Upwind:
        dudx_adv = (u > 0.0)
            ? u * (u - g.u[g.idx(i-1, j)]) / g.dx
            : u * (g.u[g.idx(i+1, j)] - u) / g.dx;
        dudy_adv = (v_at_u > 0.0)
            ? v_at_u * (u - g.u[g.idx(i, j-1)]) / g.dy
            : v_at_u * (g.u[g.idx(i, j+1)] - u) / g.dy;
        break;

    case AdvectionScheme::Hybrid: {
        double peclet_x = std::abs(u) * g.dx / 1e-6;
        double peclet_y = std::abs(v_at_u) * g.dy / 1e-6;
        if (peclet_x < 2.0) {
            dudx_adv = u * (g.u[g.idx(i+1, j)] - g.u[g.idx(i-1, j)]) / (2.0 * g.dx);
        } else {
            dudx_adv = (u > 0.0)
                ? u * (u - g.u[g.idx(i-1, j)]) / g.dx
                : u * (g.u[g.idx(i+1, j)] - u) / g.dx;
        }
        if (peclet_y < 2.0) {
            dudy_adv = v_at_u * (g.u[g.idx(i, j+1)] - g.u[g.idx(i, j-1)]) / (2.0 * g.dy);
        } else {
            dudy_adv = (v_at_u > 0.0)
                ? v_at_u * (u - g.u[g.idx(i, j-1)]) / g.dy
                : v_at_u * (g.u[g.idx(i, j+1)] - u) / g.dy;
        }
        break;
    }
    }

    return dudx_adv + dudy_adv;
}

double convective_v(const FDMGrid& g, int i, int j, AdvectionScheme scheme) {
    double v = g.v[g.idx(i, j)];
    double u_at_v = 0.25 * (g.u[g.idx(i, j-1)] + g.u[g.idx(i+1, j-1)] +
                             g.u[g.idx(i, j)] + g.u[g.idx(i+1, j)]);

    double dvdx_adv, dvdy_adv;

    switch (scheme) {
    case AdvectionScheme::Central:
        dvdx_adv = u_at_v * (g.v[g.idx(i+1, j)] - g.v[g.idx(i-1, j)]) / (2.0 * g.dx);
        dvdy_adv = v * (g.v[g.idx(i, j+1)] - g.v[g.idx(i, j-1)]) / (2.0 * g.dy);
        break;

    case AdvectionScheme::Upwind:
        dvdx_adv = (u_at_v > 0.0)
            ? u_at_v * (v - g.v[g.idx(i-1, j)]) / g.dx
            : u_at_v * (g.v[g.idx(i+1, j)] - v) / g.dx;
        dvdy_adv = (v > 0.0)
            ? v * (v - g.v[g.idx(i, j-1)]) / g.dy
            : v * (g.v[g.idx(i, j+1)] - v) / g.dy;
        break;

    case AdvectionScheme::Hybrid: {
        double peclet_x = std::abs(u_at_v) * g.dx / 1e-6;
        double peclet_y = std::abs(v) * g.dy / 1e-6;
        if (peclet_x < 2.0) {
            dvdx_adv = u_at_v * (g.v[g.idx(i+1, j)] - g.v[g.idx(i-1, j)]) / (2.0 * g.dx);
        } else {
            dvdx_adv = (u_at_v > 0.0)
                ? u_at_v * (v - g.v[g.idx(i-1, j)]) / g.dx
                : u_at_v * (g.v[g.idx(i+1, j)] - v) / g.dx;
        }
        if (peclet_y < 2.0) {
            dvdy_adv = v * (g.v[g.idx(i, j+1)] - g.v[g.idx(i, j-1)]) / (2.0 * g.dy);
        } else {
            dvdy_adv = (v > 0.0)
                ? v * (v - g.v[g.idx(i, j-1)]) / g.dy
                : v * (g.v[g.idx(i, j+1)] - v) / g.dy;
        }
        break;
    }
    }

    return dvdx_adv + dvdy_adv;
}

} // namespace advection
} // namespace exd::physics::fluid::fdm
