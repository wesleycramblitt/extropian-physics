#pragma once

// -----------------------------------------------------
// Single axial turbomachinery stage, mean-line Euler-work
// closure (W9).
//
// Total-state bookkeeping through a stage:
//   1. Density-velocity fixed point at the rotor leading
//      edge (compressible):  c_a = mdot/(rho*A) iterated to
//      consistency with the isentropic static inlet state.
//   2. Euler work:  delta_h0 = u * (c_w2 - c_w1) with the
//      absolute exit swirl c_w2 = u + c_a*tan(beta_2).
//      The sign of delta_h0 is geometry-emergent: positive
//      for work-in (compression), negative for work-out
//      (expansion). No mode flags anywhere.
//   3. Polytropic total-pressure bookkeeping through
//      exd::engine::physics::thermo::polytropic with a single stage efficiency.
//   4. Static exit from the total state, and relative Mach
//      at the rotor LE for the choked/envelope warning.
//
// See the mapping doc (docs/) for the closure derivation.
// -----------------------------------------------------

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/physics/thermo/eos.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::turbomachinery
{

struct StageGeometryConfig
{
    double r_hub = 0.18;       // hub radius, m (> 0)
    double r_tip = 0.22;       // tip radius, m (> r_hub)
    double alpha_1_rad = 0.35; // absolute inlet swirl angle at rotor LE, rad
    double beta_2_rad = 0.35;  // relative exit angle at rotor TE, rad

    double r_mean() const { return 0.5 * (r_hub + r_tip); }
    double flow_area() const { return M_PI * (r_tip * r_tip - r_hub * r_hub); }
    double hub_tip_ratio() const { return r_tip > 0.0 ? r_hub / r_tip : 0.0; }
};

struct StageLossConfig
{
    double polytropic_efficiency = 0.85; // in (0, 1]
};

struct StageConfig
{
    StageGeometryConfig geometry;
    StageLossConfig loss;
};

/// Validate a stage config. Fatal problems return false with `error`;
/// non-fatal observations are pushed to `warnings` (envelope guard).
bool validate_stage_config(const StageConfig& config,
                           std::string& error,
                           std::vector<std::string>& warnings);

/// Inlet state to a stage. p0/T0 are TOTAL state; c_theta is the
/// absolute swirl entering the rotor (m/s; carried between stages).
struct StageInlet
{
    double p0 = 101325.0; // total pressure, Pa
    double T0 = 288.15;   // total temperature, K
    double c_theta = 0.0; // absolute swirl, m/s
};

struct StageResult
{
    bool ok = false;
    exd::engine::core::ModelStatus status;

    double p0_out = 0.0;     // total pressure at stage exit, Pa
    double T0_out = 0.0;     // total temperature at stage exit, K
    double c_theta_out = 0.0;// absolute swirl at stage exit, m/s
    double static_p = 0.0;   // static pressure at exit, Pa
    double static_T = 0.0;   // static temperature at exit, K
    double static_rho = 0.0; // static density at exit, kg/m^3
    double delta_h0 = 0.0;   // specific total-enthalpy change, J/kg
    double torque = 0.0;     // shaft torque, N*m
    double power = 0.0;      // shaft power, W
    double pi = 0.0;         // total-pressure ratio p0_out/p0_in
    double tau = 0.0;        // total-temperature ratio T0_out/T0_in
    double mach_rel_le = 0.0;// relative Mach at rotor leading edge
    double work_coefficient = 0.0;  // delta_h0/u^2
    double flow_coefficient = 0.0;  // c_a/u
    bool choked = false;     // mach_rel_le >= 1 at rotor LE
};

/// Solve one stage. `omega` >= 0 (0 allowed), `mdot` > 0. `eos` is a
/// non-owning reference used for density and specific heats. Outcome and
/// warnings are mirrored into `status` (also stored on the result).
StageResult solve_stage(const StageConfig& config,
                        const StageInlet& inlet,
                        double omega,
                        double mdot,
                        const exd::engine::physics::thermo::IEos& eos,
                        exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::turbomachinery