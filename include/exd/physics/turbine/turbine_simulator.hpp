#pragma once

#include <exd/geometry/turbine.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/fluid/forces/force_evaluator.hpp>
#include <exd/physics/fluid/forces/flow_types.hpp>
#include <exd/physics/mechanics/dynamics.hpp>
#include <exd/physics/mechanics/moment_model.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>
#include <exd/physics/mechanics/status.hpp>

#include <string>
#include <vector>

namespace exd::physics::turbine {

// ─────────────────────────────────────────────────────
// Turbine application assembly.
//
// This is the ONLY turbine-specific layer: it maps an
// exd-geometry TurbineDefinition onto the generic
// fluid::forces + mechanics modules and wires them into
// a coupled loop:
//
//   field → sample blade surface → force evaluator
//        → integrate moments → generator moment
//        → rotational dynamics → new state
//
// The same generic modules drive propellers, pumps, fans…
// without this file.
// ─────────────────────────────────────────────────────

struct TurbineConfig
{
    // Blade discretization
    int element_count = 32;           // ≥ 4
    std::string default_airfoil = "naca0012";

    // Force evaluation. Default MomentumBalance = reduced-order standalone;
    // for CFD-coupled runs set force.type = ForceEvaluatorType::PressureIntegration
    // (with a field carrying the real pressure field).
    fluid::forces::ForceEvaluatorParams force;

    // Generator / load: opposing torque vs ω curve. Empty = no load.
    mechanics::CurveMomentConfig generator;

    // Drivetrain
    double inertia = 1000.0;          // kg·m² about the axis (> 0)
    mechanics::RotationalIntegration integration = mechanics::RotationalIntegration::Heun;

    // Time stepping (simulate_turbine)
    double dt = 0.01;                 // s (> 0)
    int max_steps = 2000;             // > 0
    double initial_omega = 0.0;       // rad/s at t = 0

    // Output
    bool record_history = true;
    int history_interval = 1;         // ≥ 1
};

struct TurbineStepResult
{
    bool ok = false;
    mechanics::ModelStatus status;
    mechanics::RotationalState state; // updated state
    mechanics::MomentResult aero;     // torque about axis + thrust
    double external_moment = 0.0;     // opposing (N·m)
    double power = 0.0;               // aero torque · ω (W)
    std::vector<mechanics::ElementForce3D> per_element;
};

struct TurbineSimResult
{
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    TurbineStepResult final_step;
    double cp = 0.0;   // power / (0.5·ρ·v³·π·R_tip²)
    double ct = 0.0;   // thrust / (0.5·ρ·v²·π·R_tip²)
    double tsr = 0.0;  // ω·R_tip / v_axial
    double total_energy = 0.0; // Σ aero_power·dt (J)
    std::vector<TurbineStepResult> history;
};

/// Map an exd-geometry TurbineDefinition onto the generic blade geometry.
/// 3D convention follows exd-geometry: rotation axis = +Z, radial = (r),
/// rotor plane at z = z_r (mid LE/TE axial coordinate).
fluid::forces::BladeGeometry make_blade_geometry(const exd::geometry::TurbineDefinition& turbine,
                                                 int element_count,
                                                 const std::string& default_airfoil,
                                                 std::vector<std::string>& warnings);

/// Generator load curve: T(ω) = P_rated / (η·ω), clamped flat below min_omega.
/// Returns a CurveMomentConfig with `points` samples on [min_omega, ...].
mechanics::CurveMomentConfig make_generator_curve(double rated_power,
                                                  double efficiency,
                                                  double min_omega,
                                                  int points = 32);

/// One coupled turbine step against a flow field (CFD or analytic/uniform).
/// Samples the blade surface at the current azimuth, evaluates forces,
/// integrates moments, applies the generator load, advances the rotor state.
/// `state` is updated in place; the return value carries the full step data.
TurbineStepResult step_turbine(mechanics::RotationalState& state,
                               const exd::geometry::TurbineDefinition& turbine,
                               const coupling::IFlowField3D& flow,
                               const TurbineConfig& config);

/// Standalone simulation against a uniform freestream field.
TurbineSimResult simulate_turbine(const exd::geometry::TurbineDefinition& turbine,
                                  const fluid::forces::Freestream& freestream,
                                  const TurbineConfig& config);

} // namespace exd::physics::turbine
