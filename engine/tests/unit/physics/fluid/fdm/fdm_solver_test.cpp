#include <doctest/doctest.h>
#include <exd/engine/physics/fluid/fdm/fdm_config.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_result.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_solver.hpp>
#include <cmath>

using namespace exd::engine::physics::fluid::fdm;

namespace {
FDMConfig channel_flow() {
    FDMConfig c;
    c.nx = 16;
    c.ny = 16;
    c.lx = 1.0;
    c.ly = 1.0;
    c.rho = 1.0;
    c.mu = 0.01;
    c.dt = 0.001;
    c.max_steps = 50;
    c.time_integration = TimeIntegration::Heun;
    c.advection_scheme = AdvectionScheme::Upwind;
    c.convergence_tolerance = 1e-4;
    c.convergence_window = 50;
    c.pressure_max_iterations = 200;
    c.pressure_tolerance = 1e-6;
    c.sor_omega = 1.2;
    c.velocity_under_relaxation = 0.5;
    c.pressure_under_relaxation = 0.1;
    c.initial_u = 0.0;
    c.initial_v = 0.0;
    c.initial_p = 0.0;
    c.boundary_conditions = {
        {BoundaryEdge::Left,   FDMBoundaryType::Inlet, 1.0, 0.0},
        {BoundaryEdge::Right,  FDMBoundaryType::Outlet},
        {BoundaryEdge::Bottom, FDMBoundaryType::Wall},
        {BoundaryEdge::Top,    FDMBoundaryType::Wall},
    };
    return c;
}
} // anonymous namespace

TEST_CASE("FDM solver: channel flow produces valid result") {
    auto config = channel_flow();
    auto result = solve_fdm(config);

    REQUIRE(result.valid);
    CHECK(result.total_steps > 0);
    CHECK(std::isfinite(result.max_velocity));
    CHECK(result.max_velocity >= 0.0);

    // Field data should be populated
    CHECK(result.field.nx == config.nx);
    CHECK(result.field.ny == config.ny);
    CHECK(result.field.u.size() == static_cast<size_t>(config.nx * config.ny));
}

TEST_CASE("FDM solver: invalid config returns error") {
    FDMConfig config;
    config.nx = 0;
    auto result = solve_fdm(config);
    CHECK_FALSE(result.valid);
    CHECK_FALSE(result.error.empty());
}

TEST_CASE("FDM solver: all advection schemes produce finite results") {
    AdvectionScheme schemes[] = {
        AdvectionScheme::Central,
        AdvectionScheme::Upwind,
        AdvectionScheme::Hybrid
    };
    for (auto scheme : schemes) {
        auto config = channel_flow();
        config.advection_scheme = scheme;
        config.max_steps = 10;
        auto result = solve_fdm(config);
        REQUIRE(result.valid);
        CHECK(std::isfinite(result.max_velocity));
    }
}

TEST_CASE("FDM solver: field coordinates are correct") {
    FDMConfig config;
    config.nx = 4;
    config.ny = 4;
    config.lx = 2.0;
    config.ly = 1.0;
    config.max_steps = 1;
    config.boundary_conditions = {
        {BoundaryEdge::Left,   FDMBoundaryType::Wall},
        {BoundaryEdge::Right,  FDMBoundaryType::Wall},
        {BoundaryEdge::Bottom, FDMBoundaryType::Wall},
        {BoundaryEdge::Top,    FDMBoundaryType::Wall},
    };
    auto result = solve_fdm(config);
    REQUIRE(result.valid);

    double dx = 2.0 / 4.0;
    CHECK(result.field.x[0] == doctest::Approx(0.5 * dx));
    CHECK(result.field.x[3] == doctest::Approx(3.5 * dx));

    double dy = 1.0 / 4.0;
    CHECK(result.field.y[0] == doctest::Approx(0.5 * dy));
    CHECK(result.field.y[3] == doctest::Approx(3.5 * dy));
}

TEST_CASE("FDM solver: lid-driven cavity") {
    FDMConfig config;
    config.nx = 16;
    config.ny = 16;
    config.lx = 1.0;
    config.ly = 1.0;
    config.rho = 1.0;
    config.mu = 0.01;
    config.dt = 0.001;
    config.max_steps = 50;
    config.time_integration = TimeIntegration::Heun;
    config.advection_scheme = AdvectionScheme::Upwind;
    config.pressure_max_iterations = 200;
    config.pressure_tolerance = 1e-6;
    config.sor_omega = 1.2;
    config.velocity_under_relaxation = 0.5;
    config.pressure_under_relaxation = 0.1;
    config.boundary_conditions = {
        {BoundaryEdge::Left,   FDMBoundaryType::Wall},
        {BoundaryEdge::Right,  FDMBoundaryType::Wall},
        {BoundaryEdge::Bottom, FDMBoundaryType::Wall},
        {BoundaryEdge::Top,    FDMBoundaryType::Inlet, 1.0, 0.0},
    };
    auto result = solve_fdm(config);
    REQUIRE(result.valid);
    CHECK(std::isfinite(result.max_velocity));
    CHECK(result.max_velocity >= 0.0);
    CHECK(result.reynolds_number >= 0.0);
}
