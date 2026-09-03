#pragma once

// ─────────────────────────────────────────────────────
// Articulated planar arm (implementation_spec §30 rigid
// bodies + §11 control; robotic-arm use case).
//
// A two-link arm with REVOLUTE joints in minimal
// coordinates q = (θ₁, θ₂), point masses at the elbow and
// tip.  Dynamics (standard 2R form):
//
//   M(q)·q̈ + C(q,q̇)·q̇ + G(q) = τ
//
// with per-joint actuator torque limits and OPTIONAL
// joint stops (soft, stiff-spring penalty at q_min/q_max).
//
// Controls: PD + gravity feedforward with integral
// anti-windup — a first-class "special joints and
// controls" primitive.  Verified: horizontal-plane free
// spin conserves kinetic energy; PD setpoint tracking
// converges with near-zero steady error; the joint stops
// bound the motion.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <cstdint>

namespace exd::engine::physics::robotics {

struct TwoLinkArmConfig
{
    double l1 = 0.5;                 // link-1 length (m)
    double l2 = 0.4;                 // link-2 length (m)
    double m1 = 1.0;                 // point mass at the elbow (kg)
    double m2 = 1.0;                 // point mass at the tip (kg)
    double gravity = 0.0;            // m/s² (0 = horizontal plane)
    double q1_min = 0.0, q1_max = 0.0;   // joint-1 stops (rad); equal → no stop
    double q2_min = 0.0, q2_max = 0.0;   // joint-2 stops (rad)
    double t1_max = 1e6;             // joint-1 actuator torque limit (N·m)
    double t2_max = 1e6;             // joint-2 actuator torque limit (N·m)
    double stop_stiffness = 500.0;   // joint-stop penalty (N·m/rad)
    double stop_damping = 20.0;      // joint-stop damping (N·m·s/rad)
};

struct ArmState
{
    double q1 = 0.0;                 // joint angles (rad)
    double q2 = 0.0;
    double dq1 = 0.0;                // joint rates (rad/s)
    double dq2 = 0.0;
};

/// PD + gravity-feedforward joint controller (per joint).
struct PidGains
{
    double kp = 100.0;               // proportional (N·m/rad)
    double kd = 10.0;                // derivative (N·m·s/rad)
    double ki = 0.0;                 // integral (N·m/(rad·s)); anti-windup clamped
    double integral_clamp = 10.0;    // anti-windup limit for the integral state
};

/// One PD-controlled joint: torque = kp·e + kd·de/dt + ki·∫e.  `integral` is
/// the controller's accumulated integral state (per joint), mutated in place
/// with anti-windup clamping.
double pd_torque(double q, double dq, double q_ref, double dq_ref,
                 const PidGains& gains, double dt, double& integral);

/// Mass matrix M(q) (2×2, symmetric).
void arm_mass(const TwoLinkArmConfig& c, double q2,
              double& m11, double& m12, double& m22);

/// Coriolis/centrifugal torque vector C(q,q̇)·q̇.
void arm_coriolis(const TwoLinkArmConfig& c, double q2, double dq1, double dq2,
                  double& c1, double& c2);

/// Gravity torque vector G(q).
void arm_gravity(const TwoLinkArmConfig& c, double q1, double q2,
                 double& g1, double& g2);

/// End-effector position (m) for forward kinematics.
void arm_forward_kinematics(const TwoLinkArmConfig& c, double q1, double q2,
                            double& x, double& y);

/// Advance one step with the given joint torques (clamped to ±t_max; joint
/// stops add a penalty when a limit is crossed).  RK4 on (q, q̇).  The
/// kinetic energy of a free horizontal-plane spin is printed in the result
/// via the caller's energy bookkeeping — this function returns the state.
/// Deterministic, no exceptions.
bool step_arm(const TwoLinkArmConfig& config, ArmState& state,
              double tau1, double tau2, double dt, core::ModelStatus& status);

/// Kinetic energy ½·q̇ᵀ·M(q)·q̇ (J) — the free-motion invariant in the
/// horizontal plane.
double arm_kinetic_energy(const TwoLinkArmConfig& c, const ArmState& s);

} // namespace exd::engine::physics::robotics
