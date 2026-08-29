#pragma once

#include "quaternion.hpp"
#include "status.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string_view>

namespace exd::physics::mechanics {

// ─────────────────────────────────────────────────────
// General 6-DOF rigid body dynamics (no external math lib).
// Translation in the world frame; rotation in the body
// frame with quaternion orientation. The 1-DOF machine-axis
// stack (RotatingState/IRotationalDynamics) stays the fast
// path for turbines/engines/compressors; this module covers
// free bodies (pistons, loads, drones, contacts in Phase J).
// ─────────────────────────────────────────────────────

struct RigidBodyState
{
    std::array<double, 3> position = {0, 0, 0};        // CoM, world (m)
    std::array<double, 4> orientation = {1, 0, 0, 0};  // unit quaternion (w,x,y,z)
    std::array<double, 3> linear_velocity = {0, 0, 0}; // world (m/s)
    std::array<double, 3> angular_velocity = {0, 0, 0};// body frame (rad/s)
};

struct RigidBodyForces
{
    std::array<double, 3> force = {0, 0, 0};  // world frame (N)
    std::array<double, 3> torque = {0, 0, 0}; // body frame (N·m)
};

struct RigidBodyConfig
{
    double mass = 0.0;                                // kg (> 0)
    std::array<double, 3> inertia_principal = {0, 0, 0}; // principal moments
                                                      // about the CoM (kg·m², each > 0)
};

enum class RigidBodyIntegration : uint8_t
{
    SymplecticEuler, // semi-implicit; bounded energy drift on Hamiltonian systems
    RK4,             // fourth-order explicit; renormalizes orientation each step
};

class IRigidBodyDynamics
{
public:
    virtual ~IRigidBodyDynamics() = default;
    virtual std::string_view name() const = 0;

    /// Advance `state` by `dt` under constant `loads`.
    virtual RigidBodyState advance(double dt, const RigidBodyForces& loads,
                                   const RigidBodyState& state,
                                   ModelStatus& status) const = 0;
};

std::unique_ptr<IRigidBodyDynamics> make_rigid_body_dynamics(
    RigidBodyIntegration method, const RigidBodyConfig& config);

// ── Convenience load helpers (world-frame force/torque sources) ──

/// Gravity: F = m · g (g usually (0, 0, -9.81) m/s²).
std::array<double, 3> force_gravity(double mass, const std::array<double, 3>& g);

/// Linear spring anchored at `anchor`: F = -k · (position - anchor).
std::array<double, 3> force_linear_spring(const std::array<double, 3>& position,
                                          const std::array<double, 3>& anchor,
                                          double stiffness);

} // namespace exd::physics::mechanics
