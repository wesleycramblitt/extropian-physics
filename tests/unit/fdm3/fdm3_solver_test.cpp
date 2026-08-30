#include <doctest/doctest.h>
#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include "../../../../src/fluid/fdm3/fdm3_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using namespace exd::physics::fluid::fdm3;
using namespace exd::physics;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Runs one SIMPLE step directly on an internal grid (mirrors the solver's
// step sequence).  Used by tests that need a custom initial field.
void run_simple_step(FDM3Grid& grid, const FDM3Config& config, double dt) {
    std::copy(grid.u.begin(), grid.u.end(), grid.u_old.begin());
    std::copy(grid.v.begin(), grid.v.end(), grid.v_old.begin());
    std::copy(grid.w.begin(), grid.w.end(), grid.w_old.begin());
    std::copy(grid.p.begin(), grid.p.end(), grid.p_old.begin());

    apply_time_integration(grid, config, dt);
    compute_pressure_rhs(grid, config, dt);
    std::fill(grid.p.begin(), grid.p.end(), 0.0);
    solve_pressure_poisson(grid, config);
    std::copy(grid.p.begin(), grid.p.end(), grid.p_prime.begin());
    std::copy(grid.p_old.begin(), grid.p_old.end(), grid.p.begin());
    correct_velocity(grid, config, dt);
    apply_boundary_conditions(grid, config);
}

// Volume-averaged kinetic energy of an interior field.
double kinetic_energy(const FDM3Grid& g) {
    double e = 0.0;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const size_t id = g.idx(i, j, k);
                e += 0.5 * (g.u[id] * g.u[id] + g.v[id] * g.v[id] + g.w[id] * g.w[id]);
            }
    return e / static_cast<double>(g.nx * g.ny * g.nz);
}

// ── fixtures ────────────────────────────────────────────────────────

FDM3Config uniform_flow_config() {
    FDM3Config c;
    c.nx = 12; c.ny = 12; c.nz = 12;
    c.lx = 1.0; c.ly = 1.0; c.lz = 1.0;
    c.rho = 1.225;
    c.mu = 1e-3;
    c.dt = 0.5 * (1.0 / 12.0) * 0.5;   // 0.5 * dx * cfl_target
    c.max_steps = 50;
    c.time_integration = TimeIntegration::ForwardEuler;
    c.advection_scheme = AdvectionScheme::Upwind;
    c.initial_u = 1.0;
    c.pressure_max_iterations = 400;
    c.pressure_tolerance = 1e-8;
    c.sor_omega = 1.5;
    c.velocity_under_relaxation = 0.7;
    c.pressure_under_relaxation = 0.3;
    c.convergence_window = 1000;       // never triggers early
    c.convergence_tolerance = 1e-12;
    c.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::Inlet, 1.0, 0.0, 0.0},
        {BoundaryFace::XMax, FDMBoundaryType::Outlet},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMax, FDMBoundaryType::Symmetry},
    };
    return c;
}

FDM3Config poiseuille_config() {
    FDM3Config c;
    c.nx = 16; c.ny = 8; c.nz = 16;
    c.lx = 1.0; c.ly = 0.5; c.lz = 1.0;
    c.rho = 1.0;
    c.mu = 0.1;
    c.dt = 0.005;
    c.max_steps = 2000;
    c.time_integration = TimeIntegration::Heun;
    c.advection_scheme = AdvectionScheme::Central;
    c.pressure_max_iterations = 400;
    c.pressure_tolerance = 1e-8;
    c.sor_omega = 1.5;
    c.velocity_under_relaxation = 0.7;
    c.pressure_under_relaxation = 0.3;
    c.convergence_window = 1000;       // never triggers early
    c.convergence_tolerance = 1e-12;
    c.initial_u = 0.0;
    c.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::FixedPressure, 0.0, 0.0, 0.0, 0.01, BoundaryFace::XMax},
        {BoundaryFace::XMax, FDMBoundaryType::FixedPressure, 0.0, 0.0, 0.0, 0.00, BoundaryFace::XMin},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMin, FDMBoundaryType::Wall},
        {BoundaryFace::ZMax, FDMBoundaryType::Wall},
    };
    return c;
}

} // anonymous namespace

// ── Uniform-flow preservation (smoke test for 3D indexing/ghosts) ───

TEST_CASE("fdm3: uniform flow is preserved exactly through the grid") {
    auto config = uniform_flow_config();
    FDM3Solver solver;
    ModelStatus status;
    REQUIRE(solver.initialize(config, status));
    REQUIRE(status.ok);

    for (int n = 0; n < 50; ++n) {
        REQUIRE(solver.step(config.dt, status));
        REQUIRE(status.ok);
    }

    const auto& field = solver.field();
    REQUIRE(field.u.size() == static_cast<size_t>(config.nx) * config.ny * config.nz);

    double max_err_u = 0.0, max_err_v = 0.0, max_err_w = 0.0;
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i) {
                const size_t id = field.index(i, j, k);
                max_err_u = std::max(max_err_u, std::abs(field.u[id] - 1.0));
                max_err_v = std::max(max_err_v, std::abs(field.v[id]));
                max_err_w = std::max(max_err_w, std::abs(field.w[id]));
            }
    CHECK(max_err_u < 0.05);
    CHECK(max_err_v < 0.05);
    CHECK(max_err_w < 0.05);
    CHECK(solver.last_step().divergence < 1e-4);
}

// ── 3D plane Poiseuille ─────────────────────────────────────────────

TEST_CASE("fdm3: pressure-driven channel develops a parabolic profile") {
    auto config = poiseuille_config();
    auto result = solve_fdm3(config);
    REQUIRE(result.valid);

    const double dp = 0.01, mu = 0.1, h = 1.0, L = 1.0;
    const double u_max_analytic = dp * h * h / (8.0 * mu * L);  // = 0.0125

    // Gather the streamwise velocity profile at the x/y mid-plane.
    const int i_mid = config.nx / 2;
    const int j_mid = config.ny / 2;
    std::vector<double> u_profile(config.nz);
    for (int k = 0; k < config.nz; ++k)
        u_profile[k] = result.field.u[result.field.index(i_mid, j_mid, k)];

    const double u_center = *std::max_element(u_profile.begin(), u_profile.end());
    CHECK(u_center == doctest::Approx(u_max_analytic).epsilon(0.20));

    // Parabolic shape: symmetric about the mid-plane and matching the
    // analytical profile at the quarter points (mirror pairs across z=0.5).
    const int k_low = 4;                    // z ≈ 0.28
    const int k_high = config.nz - 1 - k_low;   // z ≈ 0.72 (mirror of k_low)
    const double z_low = result.field.z[k_low];
    const double z_high = result.field.z[k_high];
    const double para_low = u_max_analytic * 4.0 * z_low * (1.0 - z_low);
    const double para_high = u_max_analytic * 4.0 * z_high * (1.0 - z_high);

    CHECK(std::abs(u_profile[k_low] - para_low) < 0.25 * u_max_analytic);
    CHECK(std::abs(u_profile[k_high] - para_high) < 0.25 * u_max_analytic);
    CHECK(std::abs(u_profile[k_low] - u_profile[k_high]) < 0.1 * u_center);
}

// ── Taylor-Green vortex decay ───────────────────────────────────────

TEST_CASE("fdm3: Taylor-Green vortex energy decays like exp(-2*nu*k^2*t)") {
    const int n = 16;
    const double u0 = 0.5;
    const double nu = 0.01;
    const double dt = 0.005;
    const int steps = 110;   // t = 0.55 s

    FDM3Config config;
    config.nx = n; config.ny = n; config.nz = n;
    config.lx = 1.0; config.ly = 1.0; config.lz = 1.0;
    config.rho = 1.0;
    config.mu = nu;
    config.dt = dt;
    config.time_integration = TimeIntegration::RK4;
    config.advection_scheme = AdvectionScheme::Central;
    config.pressure_max_iterations = 600;
    config.pressure_tolerance = 1e-9;
    config.sor_omega = 1.5;
    // The energy-decay measurement demands a full projection each step:
    // under-relaxed velocity corrections (alpha_u < 1) leave a residual
    // divergence that systematically distorts the decay rate of a free
    // decaying mode.  Other tests use the config default (0.7).
    config.velocity_under_relaxation = 1.0;
    config.pressure_under_relaxation = 0.3;
    config.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::XMax},
        {BoundaryFace::XMax, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::XMin},
        {BoundaryFace::YMin, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::YMax},
        {BoundaryFace::YMax, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::YMin},
        {BoundaryFace::ZMin, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::ZMax},
        {BoundaryFace::ZMax, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::ZMin},
    };

    // Divergence-free single-mode field (u,v,w mode k in every direction).
    FDM3Grid grid;
    initialize_grid(grid, config);
    const double k = 2.0 * kPi;
    for (int kk = 1; kk <= n; ++kk)
        for (int j = 1; j <= n; ++j)
            for (int i = 1; i <= n; ++i) {
                const double x = (i - 0.5) * grid.dx;
                const double y = (j - 0.5) * grid.dy;
                const double z = (kk - 0.5) * grid.dz;
                const size_t id = grid.idx(i, j, kk);
                grid.u[id] = u0 * std::cos(k * x) * std::sin(k * y) * std::sin(k * z);
                grid.v[id] = u0 * std::sin(k * x) * std::cos(k * y) * std::sin(k * z);
                grid.w[id] = -2.0 * u0 * std::sin(k * x) * std::sin(k * y) * std::cos(k * z);
            }
    apply_boundary_conditions(grid, config);

    const double e0 = kinetic_energy(grid);
    REQUIRE(e0 > 0.0);

    double t = 0.0;
    for (int s = 0; s < steps; ++s) {
        run_simple_step(grid, config, dt);
        t += dt;
    }

    const double et = kinetic_energy(grid);
    const double k2 = 3.0 * k * k;
    const double expected = e0 * std::exp(-2.0 * nu * k2 * t);
    const double measured = et / expected;
    INFO("E0=", e0, " Et=", et, " expected=", expected, " measured/expected=", measured);
    CHECK(measured > 0.85);
    CHECK(measured < 1.15);
}

// ── fdm3 field adapter (coupling seam) ──────────────────────────────

TEST_CASE("fdm3 field adapter trilinearly samples cell-centered data") {
    auto config = uniform_flow_config();
    FDM3Solver solver;
    ModelStatus status;
    REQUIRE(solver.initialize(config, status));
    REQUIRE(solver.step(config.dt, status));

    auto field = exd::physics::coupling::make_fdm3_field_adapter(solver);
    REQUIRE(field != nullptr);

    const auto& f = solver.field();
    std::array<double, 3> v{0.0, 0.0, 0.0};
    double p = 0.0;

    // Sampling exactly on cell centers returns the stored values.
    for (int i = 0; i < config.nx; i += 3)
        for (int j = 0; j < config.ny; j += 3)
            for (int k = 0; k < config.nz; k += 3) {
                REQUIRE(field->sample({f.x[i], f.y[j], f.z[k]}, v, p));
                CHECK(v[0] == doctest::Approx(f.u[f.index(i, j, k)]));
                CHECK(v[1] == doctest::Approx(f.v[f.index(i, j, k)]));
                CHECK(v[2] == doctest::Approx(f.w[f.index(i, j, k)]));
            }

    // Interior, off-node sample is smooth (uniform field => 1.0).
    REQUIRE(field->sample({0.42, 0.33, 0.55}, v, p));
    CHECK(v[0] == doctest::Approx(1.0).epsilon(1e-6));

    // Out-of-bounds in any axis reports false.
    CHECK_FALSE(field->sample({-0.1, 0.5, 0.5}, v, p));
    CHECK_FALSE(field->sample({1.1, 0.5, 0.5}, v, p));
    CHECK_FALSE(field->sample({0.5, -0.1, 0.5}, v, p));
    CHECK_FALSE(field->sample({0.5, 0.5, 1.1}, v, p));

    CHECK(field->density() == doctest::Approx(config.rho));
    CHECK(field->viscosity() == doctest::Approx(config.mu));
}