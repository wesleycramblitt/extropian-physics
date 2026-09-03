// thermal_solver_test.cpp
// Phase I thermal domain: 1D conduction linearity, parabolic source profile,
// advection upwind sanity, insulated zero-flux, 3D linear-z sanity and the
// grid-channel adapter smoke test.

#include <exd/engine/physics/thermal/thermal_solver.hpp>
#include <exd/engine/coupling/field_channels.hpp>   // test-only: wrap the result grid

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace exd::engine;
using namespace exd::engine::physics::thermal;
using exd::engine::coupling::StructuredScalarGrid;

namespace {

// 1D-in-3D slab configured for conduction tests along x: nx interior nodes
// across 1.0 m (dx = 0.05), thin transverse axes (ny = nz = 2, insulated).
ThermalConfig make_slab_config()
{
    ThermalConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 2, 2};
    config.material.conductivity = 50.0;
    config.material.density = 7800.0;
    config.material.specific_heat = 500.0;
    config.boundary_kind = {
        ThermalBoundaryKind::FixedValue,  // +x
        ThermalBoundaryKind::FixedValue,  // -x
        ThermalBoundaryKind::Insulated,   // +y
        ThermalBoundaryKind::Insulated,   // -y
        ThermalBoundaryKind::Insulated,   // +z
        ThermalBoundaryKind::Insulated,   // -z
    };
    config.boundary_values = {300.0, 300.0, 300.0, 300.0, 300.0, 300.0};
    return config;
}

// Value at the node (i, j, k) of the result grid.
double node_value(const ThermalResult& result, int i, int j, int k)
{
    const int nx = result.temperature.dims[0];
    const int ny = result.temperature.dims[1];
    const std::size_t I = static_cast<std::size_t>(i + nx * (j + ny * k));
    return result.temperature.values[I];
}

} // anonymous namespace

TEST_CASE("thermal: pure conduction between fixed walls is linear")
{
    ThermalConfig config = make_slab_config();
    config.boundary_values = {300.0, 400.0, 300.0, 300.0, 300.0, 300.0}; // +x=300, -x=400
    config.source_density = 0.0;

    ModelStatus status;
    const ThermalResult result = solve_thermal(config, status);
    REQUIRE(result.ok);
    REQUIRE(status.ok);

    // L = 1.0 m, T(x) = 400 - 100*x/L.  Check 3 interior nodes.
    CHECK(node_value(result, 5, 0, 0) == doctest::Approx(375.0).epsilon(0.5 / 375.0));
    CHECK(node_value(result, 10, 0, 0) == doctest::Approx(350.0).epsilon(0.5 / 350.0));
    CHECK(node_value(result, 15, 0, 0) == doctest::Approx(325.0).epsilon(0.5 / 325.0));

    // Full linearity within 0.5 K along the whole line.
    double worst = 0.0;
    for (int i = 0; i < 21; ++i)
    {
        const double analytic = 400.0 - 100.0 * (0.05 * i);
        worst = std::max(worst, std::abs(node_value(result, i, 0, 0) - analytic));
    }
    CHECK(worst < 0.5);

    // No transverse gradient (uniform across y/z).
    CHECK(node_value(result, 10, 1, 1) == doctest::Approx(node_value(result, 10, 0, 0)).epsilon(1e-9));
}

TEST_CASE("thermal: uniform source + insulated sides produces the parabolic profile")
{
    ThermalConfig config = make_slab_config();
    config.source_density = 8000.0;      // W/m^3
    config.tolerance = 1e-13;            // keep the symmetric pair tight

    ModelStatus status;
    const ThermalResult result = solve_thermal(config, status);
    REQUIRE(result.ok);

    // T(x) = 300 + qdot/(2k)*x*(L-x), L = 1.0 m, k = 50, qdot = 8000.
    // Vertex (center, x = 0.5) = 320 K.
    const double center = node_value(result, 10, 0, 0);
    CHECK(center == doctest::Approx(320.0).epsilon(0.01));   // within 1%

    // Symmetry: T(0.25L) == T(0.75L).
    const double q = 8000.0, k = 50.0;
    const double expected025 = 300.0 + q / (2.0 * k) * (0.25 * 0.75);
    CHECK(node_value(result, 5, 0, 0) == doctest::Approx(expected025).epsilon(1e-6));
    CHECK(std::abs(node_value(result, 5, 0, 0) - node_value(result, 15, 0, 0)) < 1e-9);

    // End nodes stay pinned at 300.
    CHECK(node_value(result, 0, 0, 0) == doctest::Approx(300.0));
    CHECK(node_value(result, 20, 0, 0) == doctest::Approx(300.0));

    // Flux bookkeeping: qdot * (nx-1)(ny-1)(nz-1) cells * dx*dy*dz.
    // 20*1*1 cells of volume 0.05^3 = 1.25e-4 -> 8000*20*1.25e-4 = 20 W.
    CHECK(result.total_power == doctest::Approx(20.0).epsilon(1e-12));
}

TEST_CASE("thermal: advection shifts the profile (upwind sanity)")
{
    // Same slab as the parabolic test but with a small +x velocity (Pe < 1).
    ThermalConfig config = make_slab_config();
    config.source_density = 8000.0;
    config.tolerance = 1e-12;

    // Pe = rho*cp*|u|*dx/(2k) = (7800*500)*u*0.05/100 = 1950*u.
    config.body_velocity = {2e-4, 0.0, 0.0};   // Pe ~= 0.39 < 1

    ModelStatus status;
    const ThermalResult adv = solve_thermal(config, status);
    REQUIRE(adv.ok);
    REQUIRE(status.warnings.empty());          // Pe < 1: no warning

    // Pure-conduction baseline for the comparison.
    ThermalConfig baseline = config;
    baseline.body_velocity = {0.0, 0.0, 0.0};
    ModelStatus status2;
    const ThermalResult cond = solve_thermal(baseline, status2);
    REQUIRE(cond.ok);

    // The peak heats the downstream side, so the center cools relative to
    // the pure-conduction vertex, and the half-plane symmetry breaks.
    const double center_adv = node_value(adv, 10, 0, 0);
    const double center_cond = node_value(cond, 10, 0, 0);
    CHECK(center_adv < center_cond);
    CHECK(center_adv > 0.0);                   // bounded/positive
    CHECK(node_value(adv, 10, 0, 0) > 300.0);  // still heated well above the walls

    const double asym = std::abs(node_value(adv, 5, 0, 0) - node_value(adv, 15, 0, 0));
    CHECK(asym > 1e-3);                        // downstream shift broke symmetry

    // Peak moves downstream: T downstream of center > T upstream symmetric point.
    CHECK(node_value(adv, 15, 0, 0) > node_value(adv, 5, 0, 0));
    // And a small u truly advects: with u=0 the two are exactly equal.
    CHECK(node_value(cond, 15, 0, 0) == doctest::Approx(node_value(cond, 5, 0, 0)).epsilon(1e-9));

    // Pe > 1 fires the numerical-diffusion warning (validation and solve).
    ThermalConfig fast = config;
    fast.body_velocity = {1e-2, 0.0, 0.0};    // Pe ~= 19.5 > 1
    std::string error;
    std::vector<std::string> warnings;
    CHECK(validate_thermal_config(fast, error, warnings));
    bool found = false;
    for (const std::string& w : warnings)
        if (w.find("advection-dominated") != std::string::npos)
            found = true;
    CHECK(found);

    ModelStatus status3;
    const ThermalResult fast_result = solve_thermal(fast, status3);
    REQUIRE(fast_result.ok);
    found = false;
    for (const std::string& w : fast_result.status.warnings)
        if (w.find("advection-dominated") != std::string::npos)
            found = true;
    CHECK(found);
}

TEST_CASE("thermal: insulated face zero-flux gives a flat profile")
{
    ThermalConfig config = make_slab_config();
    // Left (-x) fixed at 350, right (+x) insulated (zero normal flux).
    config.boundary_kind = {
        ThermalBoundaryKind::Insulated,   // +x
        ThermalBoundaryKind::FixedValue,  // -x
        ThermalBoundaryKind::Insulated,   // +y
        ThermalBoundaryKind::Insulated,   // -y
        ThermalBoundaryKind::Insulated,   // +z
        ThermalBoundaryKind::Insulated,   // -z
    };
    config.boundary_values = {300.0, 350.0, 300.0, 300.0, 300.0, 300.0};
    config.source_density = 0.0;

    ModelStatus status;
    const ThermalResult result = solve_thermal(config, status);
    REQUIRE(result.ok);

    // No source, zero flux on the right: the only steady profile is a flat
    // 350 everywhere (max - min < 0.01 K).
    const double t_min = result.min_temperature;
    const double t_max = result.max_temperature;
    CHECK(t_max - t_min < 0.01);
    CHECK(t_min == doctest::Approx(350.0).epsilon(0.01 / 350.0));
}

TEST_CASE("thermal: 3D steady conduction linear along z")
{
    // 11^3 box.  Bottom (z=0) fixed 300, top (z=Lz) fixed 400; the side faces
    // are insulated so the field is 1D in z and the profile is exactly
    // T(z) = 300 + 100*z/Lz.  (With side faces FIXED the finite box would
    // pull the center below 350; the linear assertion requires the
    // insulated sides used here.)
    ThermalConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {11, 11, 11};
    config.boundary_kind = {
        ThermalBoundaryKind::Insulated,   // +x
        ThermalBoundaryKind::Insulated,   // -x
        ThermalBoundaryKind::Insulated,   // +y
        ThermalBoundaryKind::Insulated,   // -y
        ThermalBoundaryKind::FixedValue,  // +z
        ThermalBoundaryKind::FixedValue,  // -z
    };
    config.boundary_values = {300.0, 300.0, 300.0, 300.0, 400.0, 300.0};

    ModelStatus status;
    const ThermalResult result = solve_thermal(config, status);
    REQUIRE(result.ok);

    // Center node (i=j=k=5): T = 300 + 100*(0.25/0.5) = 350.
    CHECK(node_value(result, 5, 5, 5) == doctest::Approx(350.0).epsilon(0.5 / 350.0));

    // Bottom corners stay pinned at 300.
    CHECK(node_value(result, 0, 0, 0) == doctest::Approx(300.0));
    CHECK(node_value(result, 10, 0, 0) == doctest::Approx(300.0));
    CHECK(node_value(result, 0, 10, 0) == doctest::Approx(300.0));
    CHECK(node_value(result, 10, 10, 0) == doctest::Approx(300.0));

    // No transverse variation at mid-height.
    CHECK(node_value(result, 0, 0, 5) == doctest::Approx(node_value(result, 5, 5, 5)).epsilon(1e-9));
    CHECK(node_value(result, 10, 10, 5) == doctest::Approx(node_value(result, 5, 5, 5)).epsilon(1e-9));
}

TEST_CASE("thermal: validation reports errors, never throws")
{
    ThermalConfig config = make_slab_config();

    // dims < 2 per axis
    config.grid.dims = {1, 3, 3};
    std::string error;
    std::vector<std::string> warnings;
    CHECK_FALSE(validate_thermal_config(config, error, warnings));
    CHECK(error.find("dims") != std::string::npos);
    ModelStatus status;
    const ThermalResult bad1 = solve_thermal(config, status);
    CHECK_FALSE(bad1.ok);
    CHECK_FALSE(status.ok);
    CHECK(status.error.find("dims") != std::string::npos);
    CHECK(bad1.temperature.values.empty());

    // conductivity <= 0
    config = make_slab_config();
    config.material.conductivity = 0.0;
    CHECK_FALSE(validate_thermal_config(config, error, warnings));
    const ThermalResult bad2 = solve_thermal(config, status);
    CHECK_FALSE(bad2.ok);
    CHECK(status.error.find("conductivity") != std::string::npos);

    // dt <= 0
    config = make_slab_config();
    config.dt = 0.0;
    CHECK_FALSE(validate_thermal_config(config, error, warnings));
    const ThermalResult bad3 = solve_thermal(config, status);
    CHECK_FALSE(bad3.ok);
    CHECK(status.error.find("dt") != std::string::npos);

    // max_steps == 0
    config = make_slab_config();
    config.max_steps = 0;
    CHECK_FALSE(validate_thermal_config(config, error, warnings));

    // spacing <= 0 and density <= 0 also rejected
    config = make_slab_config();
    config.grid.spacing = {0.0, 0.05, 0.05};
    CHECK_FALSE(validate_thermal_config(config, error, warnings));
    config = make_slab_config();
    config.material.density = -1.0;
    CHECK_FALSE(validate_thermal_config(config, error, warnings));
}

TEST_CASE("thermal: channel adapter smoke test")
{
    // Linear conduction field as the channel payload.
    ThermalConfig config = make_slab_config();
    config.boundary_values = {300.0, 400.0, 300.0, 300.0, 300.0, 300.0};
    ModelStatus status;
    const ThermalResult result = solve_thermal(config, status);
    REQUIRE(result.ok);

    auto channel = exd::engine::coupling::make_scalar_grid_field(result.temperature);
    REQUIRE(channel != nullptr);

    double value = 0.0;
    // Sample exactly on a node: analytic value 400 - 100*(0.5/1.0) = 350.
    CHECK(channel->sample({0.5, 0.0, 0.0}, value));
    CHECK(value == doctest::Approx(node_value(result, 10, 0, 0)).epsilon(1e-12));
    CHECK(value == doctest::Approx(350.0).epsilon(0.5 / 350.0));

    // Sample near the middle of a cell: trilinear value is the linear
    // interpolation of the bracketing nodes, hence monotone between them.
    // Cell from x=0.10 (node 2) to x=0.15 (node 3): sample at x=0.125, f=0.5.
    CHECK(channel->sample({0.125, 0.0, 0.0}, value));
    const double lo = node_value(result, 2, 0, 0);
    const double hi = node_value(result, 3, 0, 0);
    CHECK(value == doctest::Approx(0.5 * (lo + hi)).epsilon(1e-12));
    CHECK(value >= std::min(lo, hi) - 1e-12);
    CHECK(value <= std::max(lo, hi) + 1e-12);
    CHECK(value == doctest::Approx(400.0 - 100.0 * (0.125 / 1.0)).epsilon(0.5 / 400.0));

    // Out of bounds is handled gracefully by the channel.
    CHECK_FALSE(channel->sample({2.0, 0.0, 0.0}, value));
}
// ─────────────────────────────────────────────────────────────
// W11: transient implicit solver + velocity-channel advection
// ─────────────────────────────────────────────────────────────

TEST_CASE("Thermal: transient CN matches the Fourier series and converges to steady")
{
    ThermalConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {21, 2, 2}; // L = 1.0 m rod
    cfg.material.conductivity = 50.0;
    cfg.material.density = 1.0;      // alpha = k/(rho*cp) = 0.05 m^2/s
    cfg.material.specific_heat = 1000.0;
    for (int f = 0; f < 6; ++f)
        cfg.boundary_kind[static_cast<size_t>(f)] = ThermalBoundaryKind::Insulated;
    cfg.boundary_kind[0] = ThermalBoundaryKind::FixedValue; // +x = 400
    cfg.boundary_kind[1] = ThermalBoundaryKind::FixedValue; // -x = 300
    cfg.boundary_values[0] = 400.0;
    cfg.boundary_values[1] = 300.0;
    cfg.initial_temperature = 300.0;
    cfg.transient = true;
    cfg.dt = 0.1;
    cfg.end_time = 2.0;   // intermediate: Fourier series check
    cfg.tolerance = 1e-8;

    ModelStatus status;
    const auto res = simulate_thermal(cfg, status);
    REQUIRE(status.ok);
    REQUIRE(res.ok);

    // Analytic: T(x,t) = 300 + 100x + sum_n Bn sin(n pi x) e^(-alpha (n pi)^2 t),
    // Bn = 200 (-1)^n / (n pi) for the 300 K initial condition.
    const double alpha = 0.05, t = 2.0, xmid = 0.5;
    double analytic = 300.0 + 100.0 * xmid;
    for (int n = 1; n <= 5; n += 2)
    {
        const double bn = 200.0 * std::pow(-1.0, n) / (n * 3.14159265358979);
        analytic += bn * std::sin(n * 3.14159265358979 * xmid) *
                    std::exp(-alpha * (n * 3.14159265358979) * (n * 3.14159265358979) * t);
    }
    const size_t mid = static_cast<size_t>(10 + 21 * (0 + 2 * 0)); // node (10,0,0)
    CHECK(res.temperature.values[mid] == doctest::Approx(analytic).epsilon(2e-3));

    // Longer run converges to the steady linear profile.
    cfg.end_time = 200.0;
    const auto res_ss = simulate_thermal(cfg, status);
    REQUIRE(status.ok);
    for (int i = 0; i <= 20; i += 5)
    {
        const double x = 0.05 * i;
        const size_t I = static_cast<size_t>(i);
        CHECK(res_ss.temperature.values[I] == doctest::Approx(300.0 + 100.0 * x).epsilon(2e-3));
    }
}

TEST_CASE("Thermal: velocity channel advects like body velocity (effective-Pe profile)")
{
    struct UniformFlow : exd::engine::coupling::IVectorField3D
    {
        std::array<double, 3> v{0.01, 0.0, 0.0};
        bool sample(const std::array<double, 3>&, std::array<double, 3>& out) const override
        {
            out = v;
            return true;
        }
    } flow;

    auto make_cfg = [&](bool use_channel) {
        ThermalConfig cfg;
        cfg.grid.origin = {0, 0, 0};
        cfg.grid.spacing = {0.05, 0.05, 0.05};
        cfg.grid.dims = {21, 2, 2};
        cfg.material.conductivity = 50.0;
        cfg.material.density = 1.0;
        cfg.material.specific_heat = 1000.0; // rho*cp = 1000
        for (int f = 0; f < 6; ++f)
            cfg.boundary_kind[static_cast<size_t>(f)] = ThermalBoundaryKind::Insulated;
        cfg.boundary_kind[0] = ThermalBoundaryKind::FixedValue;
        cfg.boundary_kind[1] = ThermalBoundaryKind::FixedValue;
        cfg.boundary_values[0] = 400.0;
        cfg.boundary_values[1] = 300.0;
        if (use_channel)
            cfg.velocity_channel = &flow;
        else
            cfg.body_velocity = flow.v;
        return cfg;
    };

    ModelStatus status;
    const auto res_body = solve_thermal(make_cfg(false), status);
    REQUIRE(status.ok);
    const auto res_chan = solve_thermal(make_cfg(true), status);
    REQUIRE(status.ok);

    // The two paths agree exactly.
    CHECK(res_body.temperature.values == res_chan.temperature.values);

    // First-order upwind adds numerical diffusion k_eff = k + rho*cp*|u|*h/2:
    // the profile matches the effective-Pe exponential.
    const double rho_cp = 1000.0, u = 0.01, k = 50.0, h = 0.05, L = 1.0;
    const double k_eff = k + rho_cp * std::fabs(u) * h / 2.0;
    const double pe = rho_cp * std::fabs(u) * L / k_eff;
    for (int i = 0; i <= 20; i += 5)
    {
        const double x = 0.05 * i;
        const double expected = 300.0 + 100.0 * (std::exp(pe * x / L) - 1.0) /
                                              (std::exp(pe) - 1.0);
        const size_t I = static_cast<size_t>(i);
        CHECK(res_body.temperature.values[I] == doctest::Approx(expected).epsilon(3e-3));
    }
}

TEST_CASE("Thermal: set_temperature_point pins a node through subsequent steps")
{
    ThermalConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {21, 2, 2};
    cfg.material.conductivity = 50.0;
    cfg.material.density = 1.0;
    cfg.material.specific_heat = 1.0; // alpha = 50
    for (int f = 0; f < 6; ++f)
        cfg.boundary_kind[static_cast<size_t>(f)] = ThermalBoundaryKind::Insulated;
    cfg.boundary_kind[0] = ThermalBoundaryKind::FixedValue;
    cfg.boundary_kind[1] = ThermalBoundaryKind::FixedValue;
    cfg.boundary_values[0] = 400.0;
    cfg.boundary_values[1] = 300.0;
    cfg.initial_temperature = 350.0;

    ThermalState state;
    ModelStatus status;
    REQUIRE(init_thermal_state(state, cfg, status));

    const std::array<double, 3> p = {0.5, 0.0, 0.0}; // node (10,0,0)
    REQUIRE(set_temperature_point(state, p, 200.0, status)); // cold spot
    REQUIRE(advance_thermal(state, 1e-3, cfg, status));
    REQUIRE(advance_thermal(state, 1e-3, cfg, status));

    const size_t I = static_cast<size_t>(10);
    CHECK(state.temperature.values[I] == doctest::Approx(200.0).epsilon(1e-12));
    // Its neighbors moved toward the spot (heat flows out of the pin).
    CHECK(state.temperature.values[I - 1] < 350.0);
    CHECK(state.temperature.values[I + 1] < 350.0);
    // Out-of-bounds write fails cleanly.
    const std::array<double, 3> bad = {5.0, 0.0, 0.0};
    CHECK(!set_temperature_point(state, bad, 200.0, status));
}
