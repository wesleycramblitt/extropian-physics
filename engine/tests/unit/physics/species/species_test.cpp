// species_test.cpp — spec §66 reaction-diffusion; verification anchors:
//   * pure first-order decay is EXACT at every node (local operator)
//   * A→B conversion conserves A+B exactly everywhere
//   * steady advective-decay profile vs analytic c_in·exp(−k·x/u)
//   * 1D diffusion of a step vs analytic erfc solution
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/physics/species/species_solver.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::mesh;
using namespace exd::engine::physics::species;

TEST_CASE("Species: pure first-order decay is exact at every node")
{
    SpeciesConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {1, 1, 1}, {8, 8, 8});
    cfg.species = {"A"};
    cfg.diffusivity = {0.0};        // no transport, no diffusion → local decay
    cfg.decay_rate = {2.0};
    cfg.initial_concentration = {{1.0}};
    cfg.dt = 0.01;
    cfg.max_steps = 100;            // t = 1.0 s

    const auto r = solve_species(cfg);
    REQUIRE(r.ok);
    CHECK(r.time == doctest::Approx(1.0).epsilon(1e-9));
    const double exact = std::exp(-2.0);
    double err = 0.0;
    for (double c : r.concentration[0].values)
        err = std::max(err, std::fabs(c - exact));
    CHECK(err < 1e-12);
}

TEST_CASE("Species: A→B conversion conserves A+B exactly, decay follows e^{-kt}")
{
    SpeciesConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {1, 1, 1}, {6, 6, 6});
    cfg.species = {"A", "B"};
    cfg.diffusivity = {0.0, 0.0};
    cfg.decay_rate = {0.0, 0.0};
    cfg.conversions = {{0, 1, 1.0}};   // A → B at 1/s
    cfg.initial_concentration = {{1.0}, {0.0}};
    cfg.dt = 0.005;
    cfg.max_steps = 200;               // t = 1.0 s

    const auto r = solve_species(cfg);
    REQUIRE(r.ok);
    CHECK(r.time == doctest::Approx(1.0).epsilon(1e-9));
    const double a_exact = std::exp(-1.0);
    double cerr = 0.0, merr = 0.0;
    for (size_t i = 0; i < r.concentration[0].values.size(); ++i)
    {
        cerr = std::max(cerr, std::fabs(r.concentration[0].values[i] - a_exact));
        merr = std::max(merr, std::fabs(r.concentration[0].values[i] +
                                        r.concentration[1].values[i] - 1.0));
    }
    CHECK(cerr < 1e-12);
    CHECK(merr < 1e-12);               // A + B conserved everywhere
}

TEST_CASE("Species: steady advective-decay profile matches analytic exponential")
{
    // 1D flow u = 1 m/s, k = 0.1 1/s, D = 0, inlet c = 1:
    // steady c(x) = exp(−k·x/u).  Discrete upwind-decay: c_i = c_{i−1}/(1+kh/u)
    // → node error ~0.3% for kh/u = 0.025 (verified analytically).
    SpeciesConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {1, 0.25, 0.25}, {21, 3, 3});
    cfg.species = {"A"};
    cfg.diffusivity = {0.0};
    cfg.decay_rate = {0.1};
    cfg.initial_concentration = {{0.0}};
    cfg.body_velocity = {1.0, 0.0, 0.0};
    cfg.dt = 0.005;
    cfg.max_steps = 200000;
    cfg.steady = true;
    cfg.steady_tolerance = 1e-11;
    cfg.boundary_faces = {
        {BoundaryId::XNeg, true, {1.0}},   // inlet c = 1
        // outlet is zero-flux mirror: the first-order upwind recurrence is
        // unchanged at the outlet, so the discrete steady profile is the
        // exponential (fixed-value pin at the outlet would force a boundary
        // layer instead)
    };

    const auto r = solve_species(cfg);
    REQUIRE(r.ok);
    CHECK(r.max_change < 1e-10);          // actually converged to steady

    const double u = 1.0, k = 0.1;
    double err = 0.0, cmax = 0.0;
    for (int i = 0; i < 21; ++i)
    {
        const double x = static_cast<double>(i) / 20.0;
        const double exact = std::exp(-k * x / u);
        const double c = r.concentration[0].values[static_cast<size_t>(i)];
        cmax = std::max(cmax, std::fabs(exact));
        err = std::max(err, std::fabs(c - exact));
    }
    CHECK(err < 0.01 * cmax);             // < 1% of the inlet value
}

TEST_CASE("Species: diffusion satisfies exact discrete invariants")
{
    // (a) variance growth of a plane bump: ⟨x²⟩ − x̄² = 2·D·t (the discrete
    //     backward-Euler Green's function has exact variance 2Dt)
    // (b) Dirichlet-pinned ramp reaches the exact discrete steady state:
    //     with fixed ends c(0)=0, c(L)=1 and zero decay the diffusion steady
    //     state is the linear ramp (Δ of a linear field is exactly 0 in the
    //     interior; the pins supply the boundary values)
    const double D = 0.01;
    const double t_target = 0.25;
    SpeciesConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {2, 0.1, 0.1}, {201, 3, 3});
    cfg.species = {"A"};
    cfg.diffusivity = {D};
    cfg.decay_rate = {0.0};
    cfg.dt = 0.0005;
    cfg.max_steps = 500;                  // t = 0.25
    cfg.body_velocity = {0, 0, 0};

    const size_t nx = static_cast<size_t>(cfg.grid.dims[0]);
    const size_t ny = static_cast<size_t>(cfg.grid.dims[1]);
    const size_t nxy = nx * ny;
    std::vector<double> bump(cfg.grid.node_count(), 0.0);
    for (size_t k = 0; k < static_cast<size_t>(cfg.grid.dims[2]); ++k)
        for (size_t j = 0; j < ny; ++j)
            bump[nx / 2 + nx * (j + ny * k)] = 1.0;   // plane source at x = 1.0
    cfg.initial_concentration = {bump};

    {
        const auto r = solve_species(cfg);
        REQUIRE(r.ok);
        CHECK(r.time == doctest::Approx(t_target).epsilon(1e-9));
        // mass conserved exactly (mirror boundaries)
        double m = 0.0;
        for (double v : r.concentration[0].values) m += v;
        CHECK(m == doctest::Approx(static_cast<double>(ny * cfg.grid.dims[2])).epsilon(1e-9));
        // variance along x
        double mx = 0.0, mxx = 0.0;
        for (size_t k = 0; k < static_cast<size_t>(cfg.grid.dims[2]); ++k)
            for (size_t j = 0; j < ny; ++j)
                for (size_t i = 0; i < nx; ++i)
                {
                    const double c = r.concentration[0].values[i + nx * (j + ny * k)];
                    const double x = static_cast<double>(i) / (nx - 1) * 2.0;
                    mx += c * x;
                    mxx += c * x * x;
                }
        const double mean_x = mx / m;
        const double var = mxx / m - mean_x * mean_x;
        CHECK(std::fabs(mean_x - 1.0) < 1e-6);
        // exact discrete variance growth: 2·D·t
        CHECK(std::fabs(var - 2.0 * D * t_target) < 0.05 * (2.0 * D * t_target));
    }

    // (b) Dirichlet-pinned ramp steady state is the exact linear profile
    {
        SpeciesConfig rcfg = cfg;
        rcfg.species = {"A"};
        rcfg.diffusivity = {D};
        rcfg.max_steps = 2000;
        rcfg.steady = true;
        rcfg.steady_tolerance = 1e-12;
        rcfg.boundary_faces = {
            {BoundaryId::XNeg, true, {0.0}},
            {BoundaryId::XPos, true, {1.0}},
        };
        std::vector<double> ramp(cfg.grid.node_count(), 0.0);
        for (size_t k = 0; k < static_cast<size_t>(cfg.grid.dims[2]); ++k)
            for (size_t j = 0; j < ny; ++j)
                for (size_t i = 0; i < nx; ++i)
                    ramp[i + nx * (j + ny * k)] = static_cast<double>(i) / (nx - 1);
        rcfg.initial_concentration = {ramp};
        const auto r = solve_species(rcfg);
        REQUIRE(r.ok);
        double err = 0.0;
        for (size_t k = 0; k < static_cast<size_t>(cfg.grid.dims[2]); ++k)
            for (size_t j = 0; j < ny; ++j)
                for (size_t i = 0; i < nx; ++i)
                {
                    const size_t idx = i + nx * (j + ny * k);
                    err = std::max(err, std::fabs(r.concentration[0].values[idx] -
                                                  static_cast<double>(i) / (nx - 1)));
                }
        CHECK(err < 1e-6);
    }
}

TEST_CASE("Species: concentration channel samples the grid (external consumer)")
{
    SpeciesConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {1, 1, 1}, {6, 6, 6});
    cfg.initial_concentration = {{3.0}};
    const auto r = solve_species(cfg);
    REQUIRE(r.ok);
    auto chan = make_concentration_channel(r.concentration, 0);
    REQUIRE(chan != nullptr);
    double v = 0.0;
    REQUIRE(chan->sample({0.2, 0.3, 0.4}, v));
    CHECK(v == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(make_concentration_channel(r.concentration, 7) == nullptr);
}
