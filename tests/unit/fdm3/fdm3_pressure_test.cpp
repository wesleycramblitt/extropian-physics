#include <doctest/doctest.h>
#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>
#include "../../../../src/fluid/fdm3/fdm3_internal.hpp"

#include <cmath>
#include <vector>

using namespace exd::physics::fluid::fdm3;

namespace {

double p_exact(double x, double y, double z) {
    return x * x + y * y + z * z;  // Laplacian = 6 everywhere
}

// Build a hand-tuned grid for the Poisson test:
//   * interior RHS = +6 (the relaxation solves Laplacian(p) = +rhs)
//   * Dirichlet ghosts set to the analytic solution at ghost cell centers
//   * interior p starts from an arbitrary initial guess
FDM3Grid make_poisson_grid(int n, double omega_tol_epoch = 0.0) {
    (void)omega_tol_epoch;
    FDM3Grid g;
    g.allocate(n, n, n);
    g.dx = 1.0 / n;
    g.dy = 1.0 / n;
    g.dz = 1.0 / n;

    for (int k = 1; k <= n; ++k)
        for (int j = 1; j <= n; ++j)
            for (int i = 1; i <= n; ++i) {
                const double x = (i - 0.5) * g.dx;
                const double y = (j - 0.5) * g.dy;
                const double z = (k - 0.5) * g.dz;
                const size_t id = g.idx(i, j, k);
                g.rhs[id] = 6.0;
                g.p[id] = 0.0;  // zero initial guess
            }

    // Dirichlet ghosts: analytic solution at the ghost cell centers.
    for (int k = 1; k <= n; ++k)
        for (int j = 1; j <= n; ++j) {
            const double y = (j - 0.5) * g.dy;
            const double z = (k - 0.5) * g.dz;
            g.p[g.idx(0, j, k)] = p_exact(-0.5 * g.dx, y, z);
            g.p[g.idx(n + 1, j, k)] = p_exact((n + 0.5) * g.dx, y, z);
        }
    for (int k = 1; k <= n; ++k)
        for (int i = 1; i <= n; ++i) {
            const double x = (i - 0.5) * g.dx;
            const double z = (k - 0.5) * g.dz;
            g.p[g.idx(i, 0, k)] = p_exact(x, -0.5 * g.dy, z);
            g.p[g.idx(i, n + 1, k)] = p_exact(x, (n + 0.5) * g.dy, z);
        }
    for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= n; ++i) {
            const double x = (i - 0.5) * g.dx;
            const double y = (j - 0.5) * g.dy;
            g.p[g.idx(i, j, 0)] = p_exact(x, y, -0.5 * g.dz);
            g.p[g.idx(i, j, n + 1)] = p_exact(x, y, (n + 0.5) * g.dz);
        }
    return g;
}

double linf_error(const FDM3Grid& g, int n) {
    double err = 0.0;
    for (int k = 1; k <= n; ++k)
        for (int j = 1; j <= n; ++j)
            for (int i = 1; i <= n; ++i) {
                const double x = (i - 0.5) * g.dx;
                const double y = (j - 0.5) * g.dy;
                const double z = (k - 0.5) * g.dz;
                err = std::max(err, std::abs(g.p[g.idx(i, j, k)] - p_exact(x, y, z)));
            }
    return err;
}

} // anonymous namespace

TEST_CASE("7-point SOR converges to the analytic quadratic on a 12^3 grid") {
    const int n = 12;

    FDM3Config config;
    config.pressure_max_iterations = 4000;
    config.pressure_tolerance = 1e-10;
    config.sor_omega = 1.5;

    FDM3Grid g = make_poisson_grid(n);
    solve_pressure_poisson(g, config);

    // The discrete second derivative of a quadratic is exact, so with exact
    // Dirichlet ghosts the converged solution matches the analytic values at
    // cell centers to solver tolerance.
    CHECK(linf_error(g, n) < 1e-4);
}

TEST_CASE("SOR converges from a nonzero initial guess and with 8^3 grid") {
    const int n = 8;

    FDM3Config config;
    config.pressure_max_iterations = 4000;
    config.pressure_tolerance = 1e-8;
    config.sor_omega = 1.2;

    FDM3Grid g = make_poisson_grid(n);
    for (int k = 1; k <= n; ++k)
        for (int j = 1; j <= n; ++j)
            for (int i = 1; i <= n; ++i)
                g.p[g.idx(i, j, k)] = 123.0;  // nonzero initial guess

    solve_pressure_poisson(g, config);
    CHECK(linf_error(g, n) < 1e-6);
}

TEST_CASE("divergence of uniform flow is zero on the collocated grid") {
    FDM3Grid g;
    g.allocate(8, 8, 8);
    g.dx = 0.125;
    g.dy = 0.125;
    g.dz = 0.125;
    // Fill every cell (ghosts included) with the uniform flow so the
    // centered differences vanish exactly.
    for (int k = 0; k <= 8 + 1; ++k)
        for (int j = 0; j <= 8 + 1; ++j)
            for (int i = 0; i <= 8 + 1; ++i) {
                const size_t id = g.idx(i, j, k);
                g.u[id] = 1.0;
                g.v[id] = 0.5;
                g.w[id] = -0.25;
            }

    double max_div = 0.0;
    for (int k = 1; k <= 8; ++k)
        for (int j = 1; j <= 8; ++j)
            for (int i = 1; i <= 8; ++i)
                max_div = std::max(max_div, std::abs(spatial::divergence(g, i, j, k)));
    CHECK(max_div < 1e-12);
}