// cht_preset_test.cpp — spec §68: "Conjugate Heat Transfer" preset through
// the full configuration pipeline (validate → coupling graph → execution
// graph → allocate state → run).  Expected steady state: exact joined linear
// profile 400 → 350 → 300 across the interface.
#include <exd/engine/coupling/pipeline.hpp>
#include <exd/engine/presets/multiphysics/conjugate_heat_transfer.hpp>
#include <doctest/doctest.h>

using namespace exd::engine;
using namespace exd::engine::coupling;
using namespace exd::engine::presets::multiphysics;

TEST_CASE("CHT preset: pipeline-driven two-slab coupling reaches 350 K interface")
{
    ModelStatus st;
    ConjugateHeatTransferConfig cfg;
    ConjugateHeatTransfer cht;
    Simulation sim;

    REQUIRE(cht.configure(cfg, sim, st));
    CHECK(sim.configured());
    REQUIRE(sim.run(st));
    CHECK(st.ok);

    // node temperature along x in slab A: node 0 = 400 K fixed, node 20 = interface
    const double t_iface = cht.state_a.temperature.values[20];
    CHECK(t_iface == doctest::Approx(350.0).epsilon(0.02));
    CHECK(cht.state_b.temperature.values[0] == doctest::Approx(350.0).epsilon(0.02));
    CHECK(cht.state_a.temperature.values[10] == doctest::Approx(375.0).epsilon(0.02));
    CHECK(cht.state_b.temperature.values[10] == doctest::Approx(325.0).epsilon(0.02));
    CHECK(sim.report().total_exchanges > 0);
}

TEST_CASE("CHT preset: state allocation + execution graph exist")
{
    ModelStatus st;
    ConjugateHeatTransfer cht;
    Simulation sim;
    REQUIRE(cht.configure(ConjugateHeatTransferConfig{}, sim, st));
    // the pipeline allocated a state (module declarations recorded)
    CHECK(sim.state().field_count() >= 2);
    // and built an execution graph with one compute node per module
    CHECK(sim.execution_graph().node_count() >= 2);
}

TEST_CASE("CHT preset: invalid config rejected before execution")
{
    ModelStatus st;
    ConjugateHeatTransfer cht;
    Simulation sim;
    ConjugateHeatTransferConfig bad;
    bad.nodes_per_slab = 1;   // dims must be >= 2
    CHECK(!cht.configure(bad, sim, st));
    CHECK(!st.ok);
}
