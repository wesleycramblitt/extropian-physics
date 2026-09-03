#include "fdm_internal.hpp"

namespace exd::engine::physics::fluid::fdm {

void FDMGrid::allocate(int nx_, int ny_) {
    nx = nx_;
    ny = ny_;
    const int sz = (nx + 2) * (ny + 2);
    u.assign(sz, 0.0);
    v.assign(sz, 0.0);
    p.assign(sz, 0.0);
    u_old.assign(sz, 0.0);
    v_old.assign(sz, 0.0);
    p_old.assign(sz, 0.0);
    u_tmp.assign(sz, 0.0);
    v_tmp.assign(sz, 0.0);
    p_prime.assign(sz, 0.0);
    rhs.assign(sz, 0.0);
}

void FDMGrid::initialize(double u0, double v0, double p0) {
    for (int j = 0; j <= ny + 1; ++j)
        for (int i = 0; i <= nx + 1; ++i) {
            size_t id = idx(i, j);
            u[id] = u0;
            v[id] = v0;
            p[id] = p0;
        }
}

} // namespace exd::engine::physics::fluid::fdm
