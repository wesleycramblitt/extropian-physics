#pragma once

// ─────────────────────────────────────────────────────
// Configuration pipeline (implementation_spec §54, §71).
//
//   Parse → Construct → Validate → Resolve Defaults →
//   Build Coupling Graph → Build Execution Graph →
//   Allocate State → Initialize Backend → Execute
//
// `Simulation` is the ultimate abstraction: modules +
// couplings + fidelity + backend resolve into one state,
// one coupling graph (CouplingManager links) and one
// execution graph (module steps + exchanges).  Both the
// expert path (build pieces by hand) and the preset path
// (select a preset, then run) resolve into this same
// runtime (spec §70).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/execution.hpp>
#include <exd/engine/core/state.hpp>
#include <exd/engine/coupling/coupling_manager.hpp>
#include <exd/engine/coupling/rules.hpp>
#include <exd/engine/fidelity/profiles.hpp>
#include <exd/engine/mesh/boundary.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace exd::engine::coupling {

/// One physical module in a simulation (spec §11, §12).
struct SimulationModule
{
    std::string name;                                  // unique across the simulation
    std::string physics;                               // e.g. "thermal"
    std::vector<std::string> supported_discretizations = {"FDM"};
    std::string preferred_discretization = "FDM";
    std::vector<core::FieldMetadata> state_requirements;  // fields to allocate
    std::vector<mesh::BoundaryCondition> boundary_conditions;
    double dt = 0.0;                                   // per-domain timestep (multi-rate)
    std::function<bool(double dt, core::ModelStatus&)> step;  // advances the module
    CouplingManager::DomainHandle handle;              // channels/sinks for coupling
};

struct SimulationConfig
{
    std::vector<SimulationModule> modules;
    std::vector<CouplingContract> contracts;
    std::string backend = "cpu";
    fidelity::FidelityProfile fidelity = fidelity::profile(fidelity::FidelityLevel::Balanced);
    double time_horizon = 0.0;
    bool implicit_coupling = false;                    // iterative/staggered (§19)
};

/// The §71 ultimate abstraction, v1: validates, allocates state from module
/// requirements, builds the coupling graph (links from contracts), builds
/// the execution graph (module steps + exchanges), and runs.
class Simulation
{
public:
    // ── pipeline ──
    bool configure(const SimulationConfig& config, core::ModelStatus& status);
    bool run(core::ModelStatus& status);
    bool step(double dt, core::ModelStatus& status);   // single macro-step (expert mode)

    // ── results/views (spec §44 external consumer contract) ──
    core::State& state() { return state_; }
    const core::State& state() const { return state_; }
    const CoupledSimulation::RunReport& report() const { return report_; }
    const core::ExecutionGraph& execution_graph() const { return graph_; }
    const RuleRegistry& rules() const { return rules_; }
    bool configured() const { return configured_; }

private:
    bool allocate_state(core::ModelStatus& status);
    bool build_links(core::ModelStatus& status);
    bool build_graph(core::ModelStatus& status);

    core::State state_;
    RuleRegistry rules_;
    SimulationConfig config_;
    std::unique_ptr<CoupledSimulation> driver_;
    core::ExecutionGraph graph_;
    CoupledSimulation::RunReport report_;
    bool configured_ = false;
    double t_ = 0.0;
};

inline bool Simulation::configure(const SimulationConfig& config, core::ModelStatus& status)
{
    config_ = config;

    // ── Validate (§55): modules unique, contracts valid, rules pass ──
    for (size_t i = 0; i < config_.modules.size(); ++i)
    {
        for (size_t j = i + 1; j < config_.modules.size(); ++j)
        {
            if (config_.modules[i].name == config_.modules[j].name)
            {
                status.ok = false;
                status.error = "simulation: duplicate module '" + config_.modules[i].name + "'";
                return false;
            }
        }
        if (config_.modules[i].name.empty() || !config_.modules[i].step)
        {
            status.ok = false;
            status.error = "simulation: module must have a name and a step function";
            return false;
        }
        // per-module discretization validation (§13)
        ValidationContext ctx;
        ctx.module_name = config_.modules[i].name;
        ctx.discretization = config_.modules[i].preferred_discretization;
        ctx.mesh_family = "structured";
        ctx.backend = config_.backend;
        DiscretizationRule dr(config_.modules[i].supported_discretizations);
        if (!dr.check(ctx, status)) return false;
    }
    for (auto& contract : config_.contracts)
    {
        ValidationContext ctx;
        ctx.contract = &contract;
        ctx.module_name = contract.source_domain + "→" + contract.destination_domain;
        ctx.backend = config_.backend;
        if (!rules_.validate(ctx, status))
            return false;
    }

    // ── Build coupling graph (§54: contracts → links) ──
    if (!build_links(status)) return false;

    // ── Allocate state from module requirements ──
    if (!allocate_state(status)) return false;

    // ── Build execution graph ──
    if (!build_graph(status)) return false;

    // ── Ready ──
    configured_ = true;
    return true;
}

inline bool Simulation::build_links(core::ModelStatus& status)
{
    driver_ = std::make_unique<CoupledSimulation>();
    CoupledSimulation::Config dcfg;
    dcfg.implicit = config_.implicit_coupling;
    dcfg.tolerance = config_.fidelity.coupling_tolerance;
    dcfg.max_sub_iterations =
        config_.fidelity.coupling_iterations > 0 ? config_.fidelity.coupling_iterations : 1;
    dcfg.relaxation = 1.0;
    driver_->set_config(dcfg);

    for (auto& m : config_.modules)
    {
        CoupledDomainSpec spec;
        spec.name = m.name;
        spec.dt = m.dt;
        spec.handle = m.handle;
        auto step_fn = m.step;
        spec.step = [step_fn](double dt) {
            core::ModelStatus st;
            return step_fn(dt, st);
        };
        if (!driver_->add_domain(spec, status)) return false;
    }
    for (auto& c : config_.contracts)
    {
        CouplingLink link;
        link.id = c.id;
        link.source_domain = c.source_domain;
        link.source_channel = c.source_quantity;
        link.target_domain = c.destination_domain;
        link.target_channel = c.destination_quantity;
        link.vector_field = (c.rank == core::FieldRank::Vector);
        link.interpolation = InterpolationMode::Nearest;
        link.interval = c.execution_interval;
        link.relaxation = c.coupling_strength;
        link.sub_iterations = std::max(1, config_.fidelity.coupling_iterations);
        link.probe_points = c.probe_points;
        if (link.probe_points.empty())
        {
            status.ok = false;
            status.error = "simulation: contract '" + c.id +
                           "' has no probe_points (spec §17: spatial mapping required)";
            return false;
        }
        if (!driver_->add_link(link, status)) return false;
    }
    return true;
}

inline bool Simulation::allocate_state(core::ModelStatus& status)
{
    state_ = core::State("simulation");
    for (auto& m : config_.modules)
    {
        state_.add_entity_set(core::EntitySet(m.name + "_entities",
                                              core::EntityKind::Cells, 0), status);
        for (auto& meta : m.state_requirements)
        {
            // size 0 → module owns its storage; state records the declaration
            core::FieldMetadata declared = meta;
            (void)declared;
            state_.add_field(meta, 0, status);
            if (!status.ok) return false;
        }
    }
    return true;
}

inline bool Simulation::build_graph(core::ModelStatus& status)
{
    graph_ = core::ExecutionGraph();
    for (auto& m : config_.modules)
    {
        core::GraphNode node;
        node.id = m.name;
        node.kind = core::NodeKind::Compute;
        auto step_fn = m.step;
        node.step = [step_fn](double, double dt, core::ModelStatus& st) { return step_fn(dt, st); };
        if (!graph_.add_node(std::move(node), status)) return false;
    }
    return graph_.finalize(status);
}

inline bool Simulation::run(core::ModelStatus& status)
{
    if (!configured_)
    {
        status.ok = false;
        status.error = "simulation: configure() must be called before run()";
        return false;
    }
    if (config_.time_horizon <= 0.0)
    {
        status.ok = false;
        status.error = "simulation: time_horizon must be > 0";
        return false;
    }
    report_ = driver_->run(config_.time_horizon, status);
    return report_.ok;
}

inline bool Simulation::step(double dt, core::ModelStatus& status)
{
    if (!configured_)
    {
        status.ok = false;
        status.error = "simulation: configure() must be called before step()";
        return false;
    }
    return graph_.execute(t_, dt, status);
}

} // namespace exd::engine::coupling
