#include "fdm3_internal.hpp"
#include <cmath>
#include <algorithm>

namespace exd::engine::physics::fluid::fdm3 {
namespace advection {

// ─── helpers ─────────────────────────────────────────────────────────

// 4-point average of component values at the corners of the interpolation
// stencil around a cell: for the two tangential axes, offsets {0,+1} in the
// +/- tangential directions (keeps the stencil symmetric, mirroring the 2D
// collocated scheme's cross-flux interpolation).
static double interp4(const FDM3Grid& g, const std::vector<double>& f,
                      int i, int j, int k, int i_off, int j_off, int k_off) {
    return 0.25 * (f[g.idx(i, j, k)] + f[g.idx(i + i_off, j, k)] +
                   f[g.idx(i, j + j_off, k)] + f[g.idx(i, j, k + k_off)]);
}

// One-dimensional convective derivative of scalar field `f` at (i,j,k) in
// the direction `axis` (0=x, 1=y, 2=z) transported by velocity `vel`.
static double convective_axis(const FDM3Grid& g, const std::vector<double>& f,
                              int i, int j, int k, int axis, double vel,
                              AdvectionScheme scheme) {
    double fc = f[g.idx(i, j, k)];
    double d_central, d_upwind;
    if (axis == 0) {
        d_central = (f[g.idx(i + 1, j, k)] - f[g.idx(i - 1, j, k)]) / (2.0 * g.dx);
        d_upwind = (vel > 0.0) ? (fc - f[g.idx(i - 1, j, k)]) / g.dx
                               : (f[g.idx(i + 1, j, k)] - fc) / g.dx;
    } else if (axis == 1) {
        d_central = (f[g.idx(i, j + 1, k)] - f[g.idx(i, j - 1, k)]) / (2.0 * g.dy);
        d_upwind = (vel > 0.0) ? (fc - f[g.idx(i, j - 1, k)]) / g.dy
                               : (f[g.idx(i, j + 1, k)] - fc) / g.dy;
    } else {
        d_central = (f[g.idx(i, j, k + 1)] - f[g.idx(i, j, k - 1)]) / (2.0 * g.dz);
        d_upwind = (vel > 0.0) ? (fc - f[g.idx(i, j, k - 1)]) / g.dz
                               : (f[g.idx(i, j, k + 1)] - fc) / g.dz;
    }

    switch (scheme) {
    case AdvectionScheme::Central:
        return vel * d_central;
    case AdvectionScheme::Upwind:
        return vel * d_upwind;
    case AdvectionScheme::Hybrid: {
        double spacing = (axis == 0) ? g.dx : (axis == 1) ? g.dy : g.dz;
        double peclet = std::abs(vel) * spacing / 1e-6;
        return vel * ((peclet < 2.0) ? d_central : d_upwind);
    }
    }
    return vel * d_central;
}

// ─── convective terms ────────────────────────────────────────────────

double convective_u(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme) {
    double u = g.u[g.idx(i, j, k)];
    // Advecting transverse velocities at the u-cell center.
    double v_at_u = interp4(g, g.v, i, j, k, 1, 0, 1);  // avg over (x,z)
    double w_at_u = interp4(g, g.w, i, j, k, 1, 1, 0);  // avg over (x,y)
    return convective_axis(g, g.u, i, j, k, 0, u, scheme) +
           convective_axis(g, g.u, i, j, k, 1, v_at_u, scheme) +
           convective_axis(g, g.u, i, j, k, 2, w_at_u, scheme);
}

double convective_v(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme) {
    double v = g.v[g.idx(i, j, k)];
    double u_at_v = interp4(g, g.u, i, j, k, 0, 1, 1);  // avg over (y,z)
    double w_at_v = interp4(g, g.w, i, j, k, 1, 1, 0);  // avg over (x,y)
    return convective_axis(g, g.v, i, j, k, 0, u_at_v, scheme) +
           convective_axis(g, g.v, i, j, k, 1, v, scheme) +
           convective_axis(g, g.v, i, j, k, 2, w_at_v, scheme);
}

double convective_w(const FDM3Grid& g, int i, int j, int k, AdvectionScheme scheme) {
    double w = g.w[g.idx(i, j, k)];
    double u_at_w = interp4(g, g.u, i, j, k, 0, 1, 1);  // avg over (y,z)
    double v_at_w = interp4(g, g.v, i, j, k, 1, 0, 1);  // avg over (x,z)
    return convective_axis(g, g.w, i, j, k, 0, u_at_w, scheme) +
           convective_axis(g, g.w, i, j, k, 1, v_at_w, scheme) +
           convective_axis(g, g.w, i, j, k, 2, w, scheme);
}

} // namespace advection
} // namespace exd::engine::physics::fluid::fdm3