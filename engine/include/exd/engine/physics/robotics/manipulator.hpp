#pragma once

// ─────────────────────────────────────────────────────
// Serial manipulator dynamics (implementation_spec §30
// rigid bodies + §11 control; robotic-arm use case).
//
// THE GENERAL CLASS — an N-link planar serial arm with
// revolute joints in minimal coordinates q ∈ ℝⁿ:
//
//   M(q)·q̈ + C(q,q̇)·q̇ + G(q) = τ
//
// configured entirely by the user's data (CAD-derived
// link lengths/masses, joint stops, actuator limits,
// controller gains).  Two links, six links, a hand on
// the end — the same solver runs.
//
// Dynamics: point masses at the distal end of each link,
// mass matrix from the Jacobians M = Σ mᵢ·Jᵢᵀ·Jᵢ,
// Coriolis vector from the Christoffel symbols
// (finite-difference ∂M/∂q), gravity from the Jacobians.
// Per-joint actuator torque limits and optional soft
// joint stops.  Controls: PD + anti-windup per joint.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <vector>

namespace exd::engine::physics::robotics {

/// One link of the serial chain: the rigid body between joint i and i+1.
struct LinkSpec
{
    double length = 0.5;   // (m)
    double mass = 1.0;     // concentrated at the distal end (kg)
};

struct SerialManipulatorConfig
{
    std::vector<LinkSpec> links;         // n links → n revolute joints
    double gravity = 0.0;                // m/s² (0 = horizontal plane)

    // per-joint limits (rad); empty entry = no stop for that joint
    std::vector<double> q_min;
    std::vector<double> q_max;
    // per-joint actuator torque limits (N·m); empty entry = unlimited
    std::vector<double> torque_max;
    double stop_stiffness = 500.0;       // N·m/rad
    double stop_damping = 20.0;          // N·m·s/rad
};

struct ManipulatorState
{
    std::vector<double> q;               // joint angles (rad)
    std::vector<double> dq;              // joint rates (rad/s)
};

/// PD + gravity-feedforward joint controller (per joint).
struct PidGains
{
    double kp = 100.0;
    double kd = 10.0;
    double ki = 0.0;
    double integral_clamp = 10.0;        // anti-windup
};

/// One PD-controlled joint torque; `integral` is the accumulated anti-windup
/// state (per joint, mutated in place).
double pd_torque(double q, double dq, double q_ref, double dq_ref,
                 const PidGains& gains, double dt, double& integral);

/// Mass matrix M(q) (n×n, symmetric, flattened row-major).
void mass_matrix(const SerialManipulatorConfig& c, const std::vector<double>& q,
                 std::vector<double>& m);

/// Coriolis/centrifugal vector C(q,q̇)·q̇.
void coriolis(const SerialManipulatorConfig& c, const std::vector<double>& q,
              const std::vector<double>& dq, std::vector<double>& cc);

/// Gravity torque vector G(q).
void gravity(const SerialManipulatorConfig& c, const std::vector<double>& q,
             std::vector<double>& g);

/// End-effector position (m).
void forward_kinematics(const SerialManipulatorConfig& c, const std::vector<double>& q,
                        double& x, double& y);

/// Advance one step with the given joint torques (RK4; torque clamped to
/// ±torque_max; joint stops add a penalty beyond the limits).  Deterministic.
bool step_manipulator(const SerialManipulatorConfig& config, ManipulatorState& state,
                      const std::vector<double>& tau, double dt,
                      exd::engine::core::ModelStatus& status);

/// Kinetic energy ½·q̇ᵀ·M·q̇ (J) — the horizontal-plane free-motion invariant.
double kinetic_energy(const SerialManipulatorConfig& c, const ManipulatorState& s);

/// Validate the config (link count, vector lengths).
bool validate_manipulator(const SerialManipulatorConfig& config,
                          exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::robotics
