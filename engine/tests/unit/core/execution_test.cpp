// execution_test.cpp — spec §38: execution graph DAG, topological order, cycles.
#include <exd/engine/core/execution.hpp>
#include <doctest/doctest.h>

using namespace exd::engine::core;

TEST_CASE("ExecutionGraph: topological order + execute")
{
    ModelStatus st;
    ExecutionGraph g;
    std::string log;

    GraphNode velocity;
    velocity.id = "compute_velocity";
    velocity.kind = NodeKind::Compute;
    velocity.step = [&](double, double dt, ModelStatus&) { log += "v"; return true; };

    GraphNode flux;
    flux.id = "compute_flux";
    flux.kind = NodeKind::Compute;
    flux.depends_on = {"compute_velocity"};
    flux.step = [&](double, double dt, ModelStatus&) { log += "f"; return true; };

    GraphNode residual;
    residual.id = "compute_residual";
    residual.kind = NodeKind::Compute;
    residual.depends_on = {"compute_flux"};
    residual.step = [&](double, double dt, ModelStatus&) { log += "r"; return true; };

    GraphNode out;
    out.id = "output";
    out.kind = NodeKind::Output;
    out.depends_on = {"compute_residual"};
    out.step = [&](double, double dt, ModelStatus&) { log += "o"; return true; };

    REQUIRE(g.add_node(velocity, st));
    REQUIRE(g.add_node(flux, st));
    REQUIRE(g.add_node(residual, st));
    REQUIRE(g.add_node(out, st));
    REQUIRE(g.finalize(st));
    REQUIRE(g.execute(0.0, 0.01, st));
    CHECK(log == "vfro");
}

TEST_CASE("ExecutionGraph: duplicate node rejected; cycle detected")
{
    ModelStatus st;
    ExecutionGraph g;
    GraphNode a;
    a.id = "a";
    a.step = [](double, double, ModelStatus&) { return true; };
    GraphNode b;
    b.id = "a";                      // duplicate
    b.step = [](double, double, ModelStatus&) { return true; };
    REQUIRE(g.add_node(a, st));
    CHECK(!g.add_node(b, st));
    CHECK(!st.ok);
    st = ModelStatus{};

    ExecutionGraph cyc;
    GraphNode c1, c2;
    c1.id = "c1"; c1.depends_on = {"c2"};
    c2.id = "c2"; c2.depends_on = {"c1"};
    c1.step = c2.step = [](double, double, ModelStatus&) { return true; };
    REQUIRE(cyc.add_node(c1, st));
    REQUIRE(cyc.add_node(c2, st));
    CHECK(!cyc.finalize(st));
    CHECK(!st.ok);

    // unknown dependency rejected
    ModelStatus st2;
    ExecutionGraph g3;
    GraphNode d;
    d.id = "d";
    d.depends_on = {"missing"};
    d.step = [](double, double, ModelStatus&) { return true; };
    REQUIRE(g3.add_node(d, st2));
    CHECK(!g3.finalize(st2));
    CHECK(!st2.ok);
}
