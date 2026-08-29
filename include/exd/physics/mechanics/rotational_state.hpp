#pragma once

#include <array>
#include <span>
#include <vector>

namespace exd::physics::mechanics {

// ─────────────────────────────────────────────────────
// 3D rotational kinematics for a body rotating about a
// fixed axis. Domain-agnostic: turbines, propellers,
// pumps, fans, flywheels, drivetrains.
// ─────────────────────────────────────────────────────

/// A rotation axis in 3D space: a point on the axis plus a direction.
struct RotationAxis
{
    std::array<double, 3> origin = {0.0, 0.0, 0.0};    // point on the axis (m)
    std::array<double, 3> direction = {0.0, 0.0, 1.0}; // axis direction (normalized on use)
};

/// 1-DOF rotational state about a fixed axis — kinematics only.
struct RotationalState
{
    double omega = 0.0;     // rad/s, signed by right-hand rule about axis direction
    double angle_rad = 0.0; // total rotation angle (rad)
};

/// Generic 3D force and moment applied to one element of a body.
struct ElementForce3D
{
    double r = 0.0;                            // spanwise radius (diagnostic, m)
    std::array<double, 3> ref = {0, 0, 0};     // application point (m)
    std::array<double, 3> force = {0, 0, 0};   // total force (N)
    std::array<double, 3> force_pressure = {0, 0, 0}; // pressure part (N, diagnostic)
    std::array<double, 3> force_shear = {0, 0, 0};    // shear part (N, diagnostic)
    std::array<double, 3> moment = {0, 0, 0};         // moment about ref (N·m)
};

/// Resultant of a set of element forces expressed about an axis.
struct MomentResult
{
    bool valid = false;                        // false when axis direction is (near) zero
    double torque = 0.0;                       // moment about the axis (N·m)
    double axial_force = 0.0;                  // force along the axis direction (N)
    std::array<double, 3> total_force = {0, 0, 0};  // full 3D resultant (N)
    std::array<double, 3> total_moment = {0, 0, 0}; // full 3D resultant moment about origin (N·m)
};

/// Normalize `v` in place. Returns false when |v| is (near) zero.
bool normalize(std::array<double, 3>& v);

/// Integrate 3D element forces into moments about an axis (pure geometry).
///
///   torque      = Σ ((ref_i − origin) × force_i) · direction
///   axial_force = Σ force_i · direction
///   total_moment about origin includes element moments (rigid-body transfer).
///
/// `MomentResult::valid` is false when the axis direction is (near) zero.
MomentResult integrate_moment(std::span<const ElementForce3D> forces,
                              const RotationAxis& axis);

} // namespace exd::physics::mechanics
