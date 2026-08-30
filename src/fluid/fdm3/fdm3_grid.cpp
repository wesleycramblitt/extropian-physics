#include "fdm3_internal.hpp"

namespace exd::physics::fluid::fdm3 {

void FDM3Grid::allocate(int nx_, int ny_, int nz_) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    const size_t sz = total();
    u.assign(sz, 0.0);
    v.assign(sz, 0.0);
    w.assign(sz, 0.0);
    p.assign(sz, 0.0);
    u_old.assign(sz, 0.0);
    v_old.assign(sz, 0.0);
    w_old.assign(sz, 0.0);
    p_old.assign(sz, 0.0);
    p_prime.assign(sz, 0.0);
    rhs.assign(sz, 0.0);
}

void FDM3Grid::initialize(double u0, double v0, double w0, double p0) {
    // Interior cells only; ghost cells stay zero (boundary conditions are
    // applied on top of the initial state by the caller).
    for (int k = 1; k <= nz; ++k)
        for (int j = 1; j <= ny; ++j)
            for (int i = 1; i <= nx; ++i) {
                size_t id = idx(i, j, k);
                u[id] = u0;
                v[id] = v0;
                w[id] = w0;
                p[id] = p0;
            }
}

void initialize_grid(FDM3Grid& g, const FDM3Config& config) {
    g.allocate(config.nx, config.ny, config.nz);
    g.dx = config.dx();
    g.dy = config.dy();
    g.dz = config.dz();
    g.initialize(config.initial_u, config.initial_v, config.initial_w, config.initial_p);
}

} // namespace exd::physics::fluid::fdm3