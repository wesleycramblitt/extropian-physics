#pragma once

// Coupled turbine-in-grid driver (Wave 4).
//
// Wires the actuator-disk-style turbine↔CFD coupling:
//
//   fdm3 solver field → disk-plane sampling (fdm3 adapter)
//     → local blade-element forces (no induction solve)
//     → per-station ring torque/thrust
//     → Gaussian-smeared body force applied to the fluid
//     → rotor rotational dynamics against the generator load
//
// Conventions (see default_grid_config below):
//   * inflow travels −Z (Inlet on ZMax with w_value = −V_inf)
//   * rotation axis = +Z (grid +Z), torque positive drives +ω
//   * rotor plane sits at grid z = rotor_origin[2]
//   * the fluid receives the NEGATED per-blade force: the disk decelerates
//     the inflow (+e_z body force) and counter-swirls it (−e_t), which spins
//     the rotor up to +ω.

#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/engine/output/field_writer.hpp>
#include <exd/engine/output/output_policy.hpp>
#include <exd/engine/physics/rigid_body/moment_model.hpp>
#include <exd/engine/core/model_status.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::presets::turbine {

// Aliases so the driver body reads naturally (fdm3:: / fdm::).
namespace fdm3 = exd::engine::physics::fluid::fdm3;
namespace fdm = exd::engine::physics::fluid::fdm;

// ─────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────

struct CoupledTurbineConfig
{
    fdm3::FDM3Config grid;                  // caller-configured BCs (see default_grid_config)
    exd::geometry::TurbineDefinition turbine;
    int element_count = 16;                 // ≥ 4 blade elements
    double rotor_inertia = 1000.0;          // kg·m² about the axis (> 0)
    exd::engine::physics::rigid_body::CurveMomentConfig generator; // opposing T(ω); empty = free spin
    std::array<double, 3> rotor_origin = {0, 0, 0}; // axis point in GRID coords
                                                        // (rotor plane z)
    uint32_t fluid_steps_per_exchange = 15;     // fluid steps per coupling exchange
    double force_relaxation = 0.4;              // β in (0,1]; under-relaxation of the body force
    double ramp_time_s = 0.5;                   // forcing ramp; must be ≥ 10·window (validated)
    double smear_cells = 2.5;                   // Gaussian sigma in cells, ≥ 1 (axial smear)
    double dt = 0.0;                            // 0 → grid.dt
    uint64_t max_steps = 10000;                 // ≥ 1 fluid steps
    bool record_history = true;
    uint64_t history_interval = 1;              // ≥ 1
    exd::engine::output::IFieldWriter* field_writer = nullptr;         // non-owning; stamped per policy
    exd::engine::output::OutputScheduler* output_scheduler = nullptr; // non-owning; null → stamps
                                                        // every field_stamp_interval
    uint32_t field_stamp_interval = 100;            // used when output_scheduler is null

    std::string csv_path;                            // rotor machine-state CSV;
                                                     // empty = no CSV output
};

// ─────────────────────────────────────────────────────
// Step / result records
// ─────────────────────────────────────────────────────

struct CoupledTurbineStep
{
    double t = 0.0;        // fluid time (s)
    double omega = 0.0;    // rotor speed (rad/s)
    double angle_rad = 0.0;
    double torque = 0.0;   // ring torque of the last exchange (N·m)
    double axial_force = 0.0; // ring axial force of the last exchange (N)
    double power = 0.0;    // torque · ω (W)
    uint64_t exchange = 0; // number of completed exchanges
};

struct CoupledTurbineResult
{
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    fdm3::FDM3Result fluid;                 // final fluid summary
    std::vector<CoupledTurbineStep> history;
    double final_omega = 0.0;
    double final_cp = 0.0;
    double final_tsr = 0.0;
    uint64_t exchanges = 0;
    double aero_work = 0.0;      // Σ aero_power·dt          (J)
    double load_work = 0.0;      // Σ load_power·dt          (J)
    double rotor_ke_change = 0.0;// 0.5·J·(ω_end² − ω_start²) (J)
    double fluid_ke_change = 0.0;// 0.5·ρ·(KE_end − KE_start) (J)
};

// ─────────────────────────────────────────────────────
// Default grid layout for this driver (documented convention):
// inflow travels −Z; Inlet at face ZMax with w_value = −v_inf, Outlet at
// ZMin, Symmetry on the 4 lateral faces, uniform init w = −v_inf.
// Box sized from a NOMINAL unit-radius disk via the multipliers — callers
// with a real turbine should override lx/ly/lz (and rotor_origin) for their
// R_tip; the driver reads whatever box the config specifies.
// ─────────────────────────────────────────────────────
fdm3::FDM3Config default_grid_config(double v_inf, int n_per_axis = 24,
                                     double domain_radius_mult = 6.0,
                                     double domain_length_mult = 8.0);

/// Runs the coupled loop: fluid steps interlaced with disk exchanges.
/// Validate in this order: grid.validate, turbine geometry via
/// make_blade_geometry, coupling-specific checks (ramp/relaxation/stability/
/// containment).  Never throws; failures are reported through `status`.
CoupledTurbineResult run_coupled_turbine(const CoupledTurbineConfig& config,
                                         ModelStatus& status);

} // namespace exd::engine::presets::turbine