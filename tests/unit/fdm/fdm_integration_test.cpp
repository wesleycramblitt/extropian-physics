#include <doctest/doctest.h>
#include <exd/physics/fluid/fdm/fdm_config.hpp>
#include <exd/physics/fluid/fdm/fdm_result.hpp>
#include <exd/physics/fluid/fdm/fdm_solver.hpp>
#include "../../../../src/fluid/fdm/fdm_internal.hpp"
#include <cmath>
#include <vector>

using namespace exd::physics::fluid::fdm;

namespace {
FDMConfig base_config() {
    FDMConfig c;
    c.nx = 16;
    c.ny = 16;
    c.lx = 1.0;
    c.ly = 1.0;
    c.rho = 1.0;
    c.mu = 0.01;
    c.dt = 0.001;
    c.max_steps = 10;
    c.convergence_tolerance = 1e-10;
    c.pressure_max_iterations = 200;
    c.pressure_tolerance = 1e-6;
    c.sor_omega = 1.2;
    c.velocity_under_relaxation = 0.5;
    c.pressure_under_relaxation = 0.1;
    c.simple_max_iterations = 20;
    c.simple_tolerance = 1e-4;
    c.boundary_conditions = {
        {BoundaryEdge::Left,   FDMBoundaryType::Inlet,  1.0, 0.0},
        {BoundaryEdge::Right,  FDMBoundaryType::Outlet},
        {BoundaryEdge::Bottom, FDMBoundaryType::Wall},
        {BoundaryEdge::Top,    FDMBoundaryType::Wall},
    };
    return c;
}
} // namespace

TEST_CASE("Forward Euler: velocity remains finite") {
    auto config = base_config();
    config.time_integration = TimeIntegration::ForwardEuler;
    config.advection_scheme = AdvectionScheme::Upwind;
    config.max_steps = 20;
    auto result = solve_fdm(config);
    REQUIRE(result.valid);
    for (const auto& h : result.history) {
        CHECK(std::isfinite(h.max_velocity));
    }
}

TEST_CASE("Heun: velocity remains finite") {
    auto config = base_config();
    config.time_integration = TimeIntegration::Heun;
    config.advection_scheme = AdvectionScheme::Upwind;
    config.max_steps = 20;
    auto result = solve_fdm(config);
    REQUIRE(result.valid);
    for (const auto& h : result.history) {
        CHECK(std::isfinite(h.max_velocity));
    }
}

TEST_CASE("RK4: velocity remains finite") {
    auto config = base_config();
    config.time_integration = TimeIntegration::RK4;
    config.advection_scheme = AdvectionScheme::Upwind;
    config.max_steps = 20;
    auto result = solve_fdm(config);
    REQUIRE(result.valid);
    for (const auto& h : result.history) {
        CHECK(std::isfinite(h.max_velocity));
    }
}

TEST_CASE("Crank-Nicolson: velocity remains finite") {
    auto config = base_config();
    config.time_integration = TimeIntegration::CrankNicolson;
    config.advection_scheme = AdvectionScheme::Upwind;
    config.max_steps = 20;
    auto result = solve_fdm(config);
    REQUIRE(result.valid);
    for (const auto& h : result.history) {
        CHECK(std::isfinite(h.max_velocity));
    }
}

TEST_CASE("All integration methods produce finite results") {
    TimeIntegration methods[] = {
        TimeIntegration::ForwardEuler,
        TimeIntegration::Heun,
        TimeIntegration::RK4,
        TimeIntegration::CrankNicolson
    };
    for (auto method : methods) {
        auto config = base_config();
        config.time_integration = method;
        config.max_steps = 10;
        auto result = solve_fdm(config);
        REQUIRE(result.valid);
        CHECK(std::isfinite(result.max_velocity));
    }
}
