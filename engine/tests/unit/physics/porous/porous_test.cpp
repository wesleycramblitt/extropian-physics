// porous_test.cpp — Darcy pressure diffusion (W14);
// verification anchors (exact discrete invariants):
//   * steady linear pressure with fixed ends (Δ of linear = 0 exactly)
//   * pressure bump variance grows at 2·K·t (discrete Green's function)
//   * no-flow reservoir with a uniform source: total mass grows linearly
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/physics/porous/porous_solver.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::mesh;
using namespace exd::engine::physics::porous;

TEST_CASE("Porous: fixed-pressure reservoir reaches the exact linear profile")
{
    PorousConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {2, 0.2, 0.2}, {41, 3, 3});
    cfg.permeability = 1e-12;
    cfg.viscosity = 1e-3;
    cfg.porosity = 0.2;
    cfg.compressibility = 1e-9;
    // implicit diffusion has no stability limit, so dt is chosen for fast
    // pseudo-transient convergence while keeping the CG well conditioned:
    // dt·K·λ₁ = 0.62 (rate 0.62/step → ~70 steps to 1e-12) and
    // dt·K·λ_max ≈ 60 (CG converges in ~30 iterations).
    cfg.dt = 0.05;
    cfg.max_steps = 2000;
    cfg.steady = true;
    cfg.steady_tolerance = 1e-12;
    cfg.initial_pressure = 0.0;
    cfg.boundary_faces = {
        {BoundaryId::XNeg, true, 1e5},   // 1 bar
        {BoundaryId::XPos, true, 0.0},
    };

    const auto r = solve_porous(cfg);
    REQUIRE(r.ok);
    CHECK(r.max_change < 1e-8);           // CG relative residual (direct steady solve)
    double err = 0.0;
    for (int i = 0; i < 41; ++i)
    {
        const double exact = 1e5 * (1.0 - static_cast<double>(i) / 40.0);
        const double p = r.pressure.values[static_cast<size_t>(i)];
        err = std::max(err, std::fabs(p - exact));
    }
    CHECK(err < 1e-6);
}

TEST_CASE("Porous: pressure bump spreads with exact variance 2·K·t")
{
    PorousConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {2, 0.2, 0.2}, {201, 5, 5});
    cfg.porosity = 0.2;
    cfg.compressibility = 5e-7;           // K = k/(φ·μ·c_t) = 0.01 m²/s
    cfg.dt = 0.001;
    cfg.max_steps = 500;                  // t = 0.5 s
    cfg.initial_pressure = 0.0;
    const double K = hydraulic_diffusivity(cfg);
    REQUIRE(K == doctest::Approx(0.01).epsilon(1e-9));  // 1e-12/(0.2·1e-3·5e-7)

    // per-node IC: plane bump at x = 1.0 (uniform in y/z)
    const size_t nx = static_cast<size_t>(cfg.grid.dims[0]);
    const size_t ny = static_cast<size_t>(cfg.grid.dims[1]);
    const size_t nz = static_cast<size_t>(cfg.grid.dims[2]);
    std::vector<double> bump(cfg.grid.node_count(), 0.0);
    for (size_t k = 0; k < nz; ++k)
        for (size_t j = 0; j < ny; ++j)
            bump[nx / 2 + nx * (j + ny * k)] = 1.0;
    cfg.initial_pressure_field = bump;

    const auto r = solve_porous(cfg);
    REQUIRE(r.ok);
    CHECK(r.time == doctest::Approx(0.5).epsilon(1e-9));
    double m = 0.0, mx = 0.0, mxx = 0.0;
    for (size_t k = 0; k < nz; ++k)
        for (size_t j = 0; j < ny; ++j)
            for (size_t i = 0; i < nx; ++i)
            {
                const double p = r.pressure.values[i + nx * (j + ny * k)];
                const double x = static_cast<double>(i) / (nx - 1) * 2.0;
                m += p; mx += p * x; mxx += p * x * x;
            }
    const double mean_x = mx / m;
    const double var = mxx / m - mean_x * mean_x;
    CHECK(std::fabs(mean_x - 1.0) < 1e-6);
    // exact discrete variance growth: 2·K·t
    CHECK(std::fabs(var - 2.0 * K * 0.5) < 0.05 * (2.0 * K * 0.5));
}

TEST_CASE("Porous: source-driven reservoir mass grows linearly")
{
    PorousConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.5, 0.1, 0.1}, {11, 3, 3});
    cfg.porosity = 0.2;
    cfg.compressibility = 1e-9;
    std::vector<double> src(cfg.grid.node_count(), 0.0);
    src[0] = 100.0;                       // injector at the origin node (1/s)
    cfg.source_rate = src;
    cfg.dt = 0.01;
    cfg.max_steps = 200;

    const auto r = solve_porous(cfg);
    REQUIRE(r.ok);
    // total mass M = φ·c_t·∫p dV grows at the injection rate: ∫q dV = 100·cell_vol
    // per second → ΔM = q·φ·c_t·... (mass added per second = q_inj·(φ c_t)·? )
    // The source adds pressure at rate q (1/s): Δp = q·dt per node → the mass
    // added per step = 100·(φ·c_t·cell_vol)·dt.  With p ≈ q·t at the source
    // node and diffusion spreading, the TOTAL mass rate = ∫q·(φc_t)·dV = 100·(φc_t·cell_vol) kg/s.
    const double cell_vol = cfg.grid.cell_volume();
    const double mdot = 100.0 * cfg.porosity * cfg.compressibility * cell_vol;
    const double expect = mdot * r.time;
    CHECK(r.total_mass == doctest::Approx(expect).epsilon(1e-3));
}
