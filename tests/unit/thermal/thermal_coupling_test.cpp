// thermal_coupling_test.cpp
// W11 acceptance: two transient thermal slabs (real domains) coupled
// through the CoupledSimulation framework with relaxed, read-back-aware
// links.  Steady state of the operator-split problem is the exact joined
// linear profile 400 -> 350 -> 300 across the interface.

#include <exd/physics/coupling/coupling_manager.hpp>
#include <exd/physics/thermal/thermal_solver.hpp>

#include <doctest/doctest.h>

#include <array>
#include <functional>
#include <memory>

using namespace exd::physics;
using namespace exd::physics::thermal;
using namespace exd::physics::coupling;

namespace {

// 1D rod config: x face pair fixed/insulated, y/z faces insulated (2-node
// axes).
ThermalConfig make_slab_config(const std::array<double, 3>& origin,
                               double left_value, double right_value,
                               bool left_fixed, bool right_fixed,
                               double initial_temperature)
{
    ThermalConfig cfg;
    cfg.grid.origin = origin;
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {21, 2, 2}; // L = 1.0 m rod
    cfg.material.conductivity = 50.0;
    cfg.material.density = 1.0;      // fast thermal diffusivity alpha = 50 m^2/s
    cfg.material.specific_heat = 1.0;
    cfg.initial_temperature = initial_temperature;
    cfg.tolerance = 1e-7;
    cfg.transient = true;
    for (int f = 0; f < 6; ++f)
        cfg.boundary_kind[static_cast<size_t>(f)] = ThermalBoundaryKind::Insulated;
    cfg.boundary_kind[0] = right_fixed ? ThermalBoundaryKind::FixedValue
                                       : ThermalBoundaryKind::Insulated; // +x
    cfg.boundary_kind[1] = left_fixed ? ThermalBoundaryKind::FixedValue
                                      : ThermalBoundaryKind::Insulated; // -x
    cfg.boundary_values[0] = right_value;
    cfg.boundary_values[1] = left_value;
    return cfg;
}

double node_temperature(const ThermalState& s, int i)
{
    const int nx = s.temperature.dims[0];
    return s.temperature.values[static_cast<size_t>(i)];
}

} // anonymous namespace

TEST_CASE("Thermal coupling: two slabs exchange interface temperature on CoupledSimulation")
{
    // Slab A: x in [0, 1], left face fixed at 400, right face = interface.
    ThermalState state_a;
    ThermalState state_b;
    ModelStatus st;

    const ThermalConfig cfg_a = make_slab_config({0, 0, 0}, 400.0, 0.0,
                                                 true, false, 400.0);
    const ThermalConfig cfg_b = make_slab_config({1, 0, 0}, 0.0, 300.0,
                                                 false, true, 300.0);
    REQUIRE(init_thermal_state(state_a, cfg_a, st));
    REQUIRE(init_thermal_state(state_b, cfg_b, st));

    // Channels + sinks + read-back for the two domains.
    auto chan_a = make_temperature_channel(state_a);
    auto chan_b = make_temperature_channel(state_b);

    CoupledDomainSpec dom_a;
    dom_a.name = "slab_a";
    dom_a.dt = 2.0e-4;
    dom_a.step = [&](double dt) { return advance_thermal(state_a, dt, cfg_a, st); };
    dom_a.handle.name = "slab_a";
    dom_a.handle.scalar_channel = [&](std::string_view) -> const IScalarField3D* {
        return chan_a.get();
    };
    dom_a.handle.scalar_write = [&](std::string_view, const std::array<double, 3>& p,
                                    double v) { return set_temperature_point(state_a, p, v, st); };
    dom_a.handle.scalar_read = [&](std::string_view, const std::array<double, 3>& p,
                                   double& v) {
        bool ok = false;
        const IScalarField3D* c = chan_a.get();
        if (c->sample(p, v)) ok = true;
        return ok;
    };

    CoupledDomainSpec dom_b;
    dom_b.name = "slab_b";
    dom_b.dt = 2.0e-4;
    dom_b.step = [&](double dt) { return advance_thermal(state_b, dt, cfg_b, st); };
    dom_b.handle.name = "slab_b";
    dom_b.handle.scalar_channel = [&](std::string_view) -> const IScalarField3D* {
        return chan_b.get();
    };
    dom_b.handle.scalar_write = [&](std::string_view, const std::array<double, 3>& p,
                                    double v) { return set_temperature_point(state_b, p, v, st); };
    dom_b.handle.scalar_read = [&](std::string_view, const std::array<double, 3>& p,
                                   double& v) {
        bool ok = false;
        const IScalarField3D* c = chan_b.get();
        if (c->sample(p, v)) ok = true;
        return ok;
    };

    // Interface coordinate shared by both slabs.
    const std::array<double, 3> iface = {1.0, 0.0, 0.0};

    CouplingLink link_a_to_b;
    link_a_to_b.id = "a_to_b";
    link_a_to_b.source_domain = "slab_a";
    link_a_to_b.source_channel = "temperature";
    link_a_to_b.target_domain = "slab_b";
    link_a_to_b.target_channel = "interface";
    link_a_to_b.relaxation = 0.7;
    link_a_to_b.sub_iterations = 20;
    link_a_to_b.probe_points = {iface};

    CouplingLink link_b_to_a;
    link_b_to_a.id = "b_to_a";
    link_b_to_a.source_domain = "slab_b";
    link_b_to_a.source_channel = "temperature";
    link_b_to_a.target_domain = "slab_a";
    link_b_to_a.target_channel = "interface";
    link_b_to_a.relaxation = 0.7;
    link_b_to_a.sub_iterations = 20;
    link_b_to_a.probe_points = {iface};

    // EXPLICIT (staggered) exchange: one relaxed, read-back-aware pass per
    // macro-step.  The per-link implicit sub-iteration would converge each
    // interface to its SOURCE value (400) instead of the coupled fixed
    // point (350); the alternating relaxed sweep across macro-steps is the
    // operator-split semantics that preserves the interface average.
    CoupledSimulation sim(CoupledSimulation::Config{
        /*implicit=*/false, /*relaxation=*/0.7, /*tolerance=*/1e-6,
        /*max_sub_iterations=*/40});
    REQUIRE(sim.add_domain(dom_a, st));
    REQUIRE(sim.add_domain(dom_b, st));
    REQUIRE(sim.add_link(link_a_to_b, st));
    REQUIRE(sim.add_link(link_b_to_a, st));

    // Diffusion time scale tau = L^2/alpha = 1^2/50 = 0.02 s; run 10 tau.
    const auto report = sim.run(0.2, st);
    REQUIRE(st.ok);
    REQUIRE(report.ok);

    // Steady operator-split state: joined linear profile 400 -> 350 -> 300.
    const double t_iface = node_temperature(state_a, 20); // x = 1.0
    CHECK(t_iface == doctest::Approx(350.0).epsilon(0.02));
    CHECK(node_temperature(state_b, 0) == doctest::Approx(350.0).epsilon(0.02));
    CHECK(node_temperature(state_a, 10) == doctest::Approx(375.0).epsilon(0.02));
    CHECK(node_temperature(state_b, 10) == doctest::Approx(325.0).epsilon(0.02));
    // The interface pins agree to the exchange resolution.
    CHECK(std::fabs(node_temperature(state_a, 20) - node_temperature(state_b, 0)) < 1e-6);
    // The coupled run actually exchanged data.
    CHECK(report.total_exchanges > 0);
}
