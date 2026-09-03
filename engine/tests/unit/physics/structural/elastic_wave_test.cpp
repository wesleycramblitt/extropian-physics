// elastic_wave_test.cpp — transient elasticity (W14): a longitudinal
// displacement pulse propagates at the P-wave speed c_p = sqrt((λ+2μ)/ρ);
// mechanical energy is conserved (velocity-Verlet on a linear system).
#include <exd/engine/physics/structural/elasticity.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::physics::structural;

namespace {
/// Long bar along x (4 m, 4×4 cross-section); gaussian u_x pulse at
/// x0 = 2.0 m.  The run is short enough that NO wave reaches a face, so the
/// probe measures the free-flight P-wave arrival and the energy conservation
/// is measured in unreflected flight.
ElasticityConfig make_wave_config(double h = 0.01, int max_steps = 300)
{
    ElasticityConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {h, h, h};
    cfg.grid.dims = {401, 4, 4};
    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.3;
    cfg.density = 7850.0;
    cfg.transient = true;
    cfg.dt = 0.0;                        // CFL-adaptive
    cfg.max_steps = max_steps;
    const size_t N = 401 * 4 * 4ull;
    cfg.initial_displacement.assign(3 * N, 0.0);
    for (int k = 0; k < 4; ++k)
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 401; ++i)
            {
                const double x = static_cast<double>(i) * h;
                const double r = (x - 2.0) / 0.02;
                const size_t idx = static_cast<size_t>(i) + 401ull * (static_cast<size_t>(j) + 4ull * k);
                cfg.initial_displacement[0 * N + idx] = 1e-5 * std::exp(-0.5 * r * r);
            }
    cfg.probe_index = 275;               // x = 2.75 m
    return cfg;
}
} // namespace

TEST_CASE("Elastic waves: pulse travels at the P-wave speed")
{
    ModelStatus st;
    ElasticityConfig cfg = make_wave_config();
    const double cp = p_wave_speed(cfg.material, cfg.density);
    // sanity: ~5856 m/s for steel
    CHECK(cp == doctest::Approx(5856.0).epsilon(0.02));

    const auto r = solve_elasticity(cfg, st);
    REQUIRE(r.ok);
    CHECK(r.steps == cfg.max_steps);
    CHECK(r.time_elapsed > 0.0);

    // arrival at the probe: time of the first max |u_x|
    REQUIRE(r.probe_history.size() == cfg.max_steps + 1);
    double arrival = 0.0;
    double peak = -1.0;
    for (size_t t = 1; t < r.probe_history.size(); ++t)
    {
        if (std::fabs(r.probe_history[t]) > peak)
        {
            peak = std::fabs(r.probe_history[t]);
            arrival = static_cast<double>(t) * r.dt_used;
        }
    }
    REQUIRE(peak > 0.0);
    const double distance = 2.75 - 2.0;
    const double expected = distance / cp;
    CHECK(std::fabs(arrival - expected) / expected < 0.05);
    CHECK(r.dt_used <= 0.3 * 0.01 / cp);   // CFL clamp was active
}

TEST_CASE("Elastic waves: mechanical energy conserved within 5%")
{
    ModelStatus st;
    const auto r = solve_elasticity(make_wave_config(), st);
    REQUIRE(r.ok);
    CHECK(r.energy_drift < 0.05);
    CHECK(r.kinetic_energy > 0.0);        // the pulse carries KE during the run
    CHECK(r.max_displacement > 0.0);
}

TEST_CASE("Elastic waves: velocity output is populated and finite")
{
    ModelStatus st;
    const auto r = solve_elasticity(make_wave_config(), st);
    REQUIRE(r.ok);
    CHECK(r.velocity.values.size() == 3 * 401 * 4 * 4ull);
    bool finite = true;
    for (double v : r.velocity.values)
        finite = finite && std::isfinite(v);
    CHECK(finite);
}
