#include <doctest/doctest.h>
#include <exd/engine/physics/fluid/fdm/fdm_config.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_result.hpp>
#include "../../../../../src/physics/fluid/fdm/fdm_internal.hpp"
#include <cmath>
#include <vector>

using namespace exd::engine::physics::fluid::fdm;

TEST_CASE("Pressure Poisson solver: uniform source field") {
    FDMGrid g;
    g.allocate(8, 8);
    g.dx = 1.0 / 8.0;
    g.dy = 1.0 / 8.0;

    std::vector<double> rhs(g.stride() * (g.ny + 2), 1.0);

    int iters = solve_pressure_poisson(g, rhs, 500, 1e-8, 1.5);
    CHECK(iters < 500);

    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i)
            CHECK(g.p[g.idx(i, j)] > 0.0);
}

TEST_CASE("Pressure Poisson solver: convergence with SOR vs Gauss-Seidel") {
    FDMGrid base;
    base.allocate(8, 8);
    base.dx = 1.0 / 8.0;
    base.dy = 1.0 / 8.0;

    std::vector<double> rhs(base.stride() * (base.ny + 2), 1.0);

    FDMGrid g_gs = base;
    int iters_gs = solve_pressure_poisson(g_gs, rhs, 1000, 1e-6, 1.0);

    FDMGrid g_sor = base;
    int iters_sor = solve_pressure_poisson(g_sor, rhs, 1000, 1e-6, 1.5);

    CHECK(iters_sor < iters_gs);
}

TEST_CASE("Divergence: zero for uniform flow") {
    FDMGrid g;
    g.allocate(8, 8);
    g.dx = 1.0 / 8.0;
    g.dy = 1.0 / 8.0;
    g.initialize(1.0, 0.0, 0.0);

    double max_div = 0.0;
    for (int j = 1; j <= g.ny; ++j)
        for (int i = 1; i <= g.nx; ++i) {
            double div = std::abs(spatial::divergence(g, i, j));
            if (div > max_div) max_div = div;
        }
    CHECK(max_div < 1e-10);
}