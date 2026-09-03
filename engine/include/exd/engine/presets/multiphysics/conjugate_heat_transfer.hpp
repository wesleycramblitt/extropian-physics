#pragma once

// ─────────────────────────────────────────────────────
// Conjugate Heat Transfer preset (implementation_spec
// §68: "Conjugate Heat Transfer" preset example).
//
// ASSEMBLY, not a solver: two transient thermal modules
// (FDM) coupled through interface temperature contracts,
// resolved by the Simulation pipeline (validate → coupling
// graph → execution graph → allocate → run).  The preset
// only declares modules, couplings, and defaults — all
// numerical work stays in the thermal module.
// ─────────────────────────────────────────────────────

#include <exd/engine/coupling/pipeline.hpp>
#include <exd/engine/physics/thermal/thermal_solver.hpp>

#include <array>
#include <functional>
#include <memory>

namespace exd::engine::presets::multiphysics {

struct ConjugateHeatTransferConfig
{
    double rod_length = 1.0;             // total length; split at the middle
    int nodes_per_slab = 21;             // node count along x (>= 2)
    double left_temperature = 400.0;     // fixed T at the left end (K)
    double right_temperature = 300.0;    // fixed T at the right end (K)
    double conductivity = 50.0;          // W/(m·K)
    double density = 1.0;                // kg/m³ (fast diffusivity for the test)
    double specific_heat = 1.0;          // J/(kg·K)
    double dt = 2.0e-4;                  // per-domain timestep (multi-rate allowed)
    double time_horizon = 0.2;           // run duration (s)
    double interface_relaxation = 0.7;   // coupling strength (relaxation)
    int sub_iterations = 20;
    bool implicit_coupling = false;      // staggered (operator-split) by default
    double solver_tolerance = 1e-7;
};

struct ConjugateHeatTransferResult
{
    bool ok = false;
    double interface_temperature = 0.0;  // converged interface temperature (K)
    double left_midpoint_temperature = 0.0;
    double right_midpoint_temperature = 0.0;
    size_t exchanges = 0;
    core::ModelStatus status;
};

/// Configure a `Simulation` for the two-slab conjugate heat transfer case.
/// The caller owns the module state (via the returned state handles) and can
/// run the simulation with Simulation::run.  This function is the preset's
/// "Parse → Construct → Validate → Build Coupling Graph → Allocate State"
/// pipeline entry.
struct ConjugateHeatTransfer
{
    physics::thermal::ThermalState state_a;
    physics::thermal::ThermalState state_b;
    physics::thermal::ThermalConfig cfg_a, cfg_b;
    std::unique_ptr<coupling::IScalarField3D> channel_a, channel_b;

    bool configure(const ConjugateHeatTransferConfig& cfg,
                   coupling::Simulation& sim,
                   core::ModelStatus& status);
};

inline bool slab_config(const ConjugateHeatTransferConfig& cfg,
                        const std::array<double, 3>& origin,
                        double left_value, double right_value,
                        bool left_fixed, bool right_fixed,
                        double initial, physics::thermal::ThermalConfig& out,
                        core::ModelStatus& status)
{
    out.grid.origin = origin;
    out.grid.spacing = {cfg.rod_length * 0.5 / (cfg.nodes_per_slab - 1),
                        cfg.rod_length * 0.5 / (cfg.nodes_per_slab - 1),
                        cfg.rod_length * 0.5 / (cfg.nodes_per_slab - 1)};
    out.grid.dims = {cfg.nodes_per_slab, 2, 2};
    out.material.conductivity = cfg.conductivity;
    out.material.density = cfg.density;
    out.material.specific_heat = cfg.specific_heat;
    out.initial_temperature = initial;
    out.tolerance = cfg.solver_tolerance;
    out.transient = true;
    for (int f = 0; f < 6; ++f)
        out.boundary_kind[static_cast<size_t>(f)] = physics::thermal::ThermalBoundaryKind::Insulated;
    out.boundary_kind[0] = right_fixed ? physics::thermal::ThermalBoundaryKind::FixedValue
                                       : physics::thermal::ThermalBoundaryKind::Insulated;  // +x
    out.boundary_kind[1] = left_fixed ? physics::thermal::ThermalBoundaryKind::FixedValue
                                      : physics::thermal::ThermalBoundaryKind::Insulated;  // -x
    out.boundary_values[0] = right_value;
    out.boundary_values[1] = left_value;
    for (int a = 0; a < 3; ++a)
    {
        if (out.grid.dims[a] < 2 || !(out.grid.spacing[a] > 0.0))
        {
            status.ok = false;
            status.error = "cht preset: invalid thermal grid (dims/spacing)";
            return false;
        }
    }
    return true;
}

inline bool ConjugateHeatTransfer::configure(const ConjugateHeatTransferConfig& cfg,
                                             coupling::Simulation& sim,
                                             core::ModelStatus& status)
{
    const double half = cfg.rod_length * 0.5;
    const std::array<double, 3> origin_a = {0, 0, 0};
    const std::array<double, 3> origin_b = {half, 0, 0};

    if (!slab_config(cfg, origin_a, cfg.left_temperature, 0.0,
                     true, false, cfg.left_temperature, cfg_a, status)) return false;
    if (!slab_config(cfg, origin_b, 0.0, cfg.right_temperature,
                     false, true, cfg.right_temperature, cfg_b, status)) return false;

    if (!physics::thermal::init_thermal_state(state_a, cfg_a, status)) return false;
    if (!physics::thermal::init_thermal_state(state_b, cfg_b, status)) return false;
    channel_a = physics::thermal::make_temperature_channel(state_a);
    channel_b = physics::thermal::make_temperature_channel(state_b);
    if (!channel_a || !channel_b) { status.ok = false; status.error = "cht preset: channels"; return false; }

    coupling::SimulationConfig simcfg;
    simcfg.time_horizon = cfg.time_horizon;
    simcfg.implicit_coupling = cfg.implicit_coupling;
    simcfg.fidelity.coupling_tolerance = 1e-6;
    simcfg.fidelity.coupling_iterations = cfg.sub_iterations;

    const std::array<double, 3> iface = {half, 0.0, 0.0};

    // ── module A: left slab, right face = interface ──
    coupling::SimulationModule mod_a;
    mod_a.name = "slab_a";
    mod_a.physics = "thermal";
    mod_a.supported_discretizations = {"FDM"};
    mod_a.preferred_discretization = "FDM";
    mod_a.dt = cfg.dt;
    // declared state requirement (§11: modules declare state; the pipeline
    // allocates it — §54 "Allocate State")
    mod_a.state_requirements.push_back(core::FieldMetadata{
        .name = "slab_a.temperature", .rank = core::FieldRank::Scalar,
        .components = 1, .units = core::units::kelvin,
        .location = core::FieldLocation::Node, .domain = "thermal"});
    mod_a.step = [this](double dt, core::ModelStatus& st) {
        return physics::thermal::advance_thermal(state_a, dt, cfg_a, st);
    };
    mod_a.handle.name = "slab_a";
    mod_a.handle.scalar_channel = [this](std::string_view) -> const coupling::IScalarField3D* {
        return channel_a.get();
    };
    mod_a.handle.scalar_write = [this](std::string_view, const std::array<double, 3>& p,
                                       double v) {
        core::ModelStatus st_local;
        return physics::thermal::set_temperature_point(state_a, p, v, st_local);
    };
    mod_a.handle.scalar_read = [this](std::string_view, const std::array<double, 3>& p,
                                      double& v) {
        return channel_a->sample(p, v);
    };

    // ── module B: right slab, left face = interface ──
    coupling::SimulationModule mod_b;
    mod_b.name = "slab_b";
    mod_b.physics = "thermal";
    mod_b.dt = cfg.dt;
    mod_b.state_requirements.push_back(core::FieldMetadata{
        .name = "slab_b.temperature", .rank = core::FieldRank::Scalar,
        .components = 1, .units = core::units::kelvin,
        .location = core::FieldLocation::Node, .domain = "thermal"});
    mod_b.step = [this](double dt, core::ModelStatus& st) {
        return physics::thermal::advance_thermal(state_b, dt, cfg_b, st);
    };
    mod_b.handle.name = "slab_b";
    mod_b.handle.scalar_channel = [this](std::string_view) -> const coupling::IScalarField3D* {
        return channel_b.get();
    };
    mod_b.handle.scalar_write = [this](std::string_view, const std::array<double, 3>& p,
                                       double v) {
        core::ModelStatus st_local;
        return physics::thermal::set_temperature_point(state_b, p, v, st_local);
    };
    mod_b.handle.scalar_read = [this](std::string_view, const std::array<double, 3>& p,
                                      double& v) {
        return channel_b->sample(p, v);
    };

    // ── contracts (spec §17: explicit coupling contracts) ──
    coupling::CouplingContract c_a2b;
    c_a2b.id = "slab_a_to_slab_b";
    c_a2b.source_domain = "slab_a";
    c_a2b.source_quantity = "temperature";
    c_a2b.destination_domain = "slab_b";
    c_a2b.destination_quantity = "interface";
    c_a2b.units = core::units::kelvin;
    c_a2b.rank = core::FieldRank::Scalar;
    c_a2b.source_association = core::FieldLocation::Node;
    c_a2b.destination_association = core::FieldLocation::Node;
    c_a2b.mapping = coupling::MappingKind::SameMesh;   // both slabs share the interface lattice
    c_a2b.interpolation = coupling::InterpolationKind::Nearest;
    c_a2b.temporal = coupling::TemporalBehavior::Explicit;
    c_a2b.coupling_strength = cfg.interface_relaxation;
    c_a2b.probe_points = {iface};

    coupling::CouplingContract c_b2a = c_a2b;
    c_b2a.id = "slab_b_to_slab_a";
    c_b2a.source_domain = "slab_b";
    c_b2a.destination_domain = "slab_a";

    simcfg.modules = {mod_a, mod_b};
    simcfg.contracts = {c_a2b, c_b2a};

    return sim.configure(simcfg, status);
}

} // namespace exd::engine::presets::multiphysics
