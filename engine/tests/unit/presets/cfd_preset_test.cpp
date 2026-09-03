// cfd_preset_test.cpp — spec §14/§52: the incompressible CFD preset builds
// a valid default channel-flow configuration and solves it.
#include <exd/engine/presets/cfd/incompressible.hpp>
#include <doctest/doctest.h>

using namespace exd::engine;
using namespace exd::engine::presets::cfd;

TEST_CASE("CFD preset: default channel flow config is valid and solves")
{
    ModelStatus st;
    IncompressibleCfdConfig cfg;
    auto flow = make_channel_flow(cfg, st);
    REQUIRE(st.ok);
    CHECK(flow.boundary_conditions.size() == 6);   // inlet + outlet + 4 walls

    const auto result = exd::engine::physics::fluid::fdm3::solve_fdm3(flow);
    REQUIRE(result.valid);
    CHECK(result.history.size() > 0);
    // settled duct: residual below the convergence tolerance
    CHECK(result.history.back().residual_u < 1e-3);
}

TEST_CASE("CFD preset: caller overrides are honored (presets are defaults)")
{
    ModelStatus st;
    IncompressibleCfdConfig cfg;
    cfg.nx = 10; cfg.ny = 4; cfg.nz = 4;
    auto flow = make_channel_flow(cfg, st);
    REQUIRE(st.ok);
    CHECK(flow.nx == 10);
}
