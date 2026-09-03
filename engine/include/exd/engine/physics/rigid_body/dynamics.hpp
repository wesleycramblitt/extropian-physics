#pragma once

#include "rotational_state.hpp"
#include "status.hpp"

#include <cstdint>
#include <memory>
#include <string_view>

namespace exd::engine::physics::rigid_body {

/// Time integration for rotational dynamics.
enum class RotationalIntegration : uint8_t
{
    Euler, // angle += ω·dt            (first order)
    Heun,  // angle += (ω + ω_new)·dt/2 (second order)
};

/// Rigid rotor about a fixed axis: J·α = M_net.
struct RigidRotorConfig
{
    double inertia = 0.0;   // kg·m² about the axis (validated > 0)
    RotationalIntegration integration = RotationalIntegration::Heun;
};

class IRotationalDynamics
{
public:
    virtual ~IRotationalDynamics() = default;
    virtual std::string_view name() const = 0;

    /// Advance `state` by `dt` under a constant `net_moment`.
    /// Positive net moment accelerates in the +axis direction.
    virtual RotationalState advance(double dt, double net_moment,
                                    const RotationalState& state,
                                    ModelStatus& status) const = 0;
};

std::unique_ptr<IRotationalDynamics> make_rigid_rotor_dynamics(const RigidRotorConfig& config);

} // namespace exd::engine::physics::rigid_body
