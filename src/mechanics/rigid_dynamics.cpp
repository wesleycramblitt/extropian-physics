// Rigid-rotor rotational dynamics about a fixed axis: J*alpha = M_net.
// Euler and Heun time integration, validated inertia and step size.

#include <exd/physics/mechanics/dynamics.hpp>

#include <memory>

namespace exd::physics::mechanics
{

namespace
{

class RigidRotorDynamics final : public IRotationalDynamics
{
public:
    explicit RigidRotorDynamics(RigidRotorConfig config) : config_(config) {}

    std::string_view name() const override { return "rigid_rotor"; }

    RotationalState advance(double dt, double net_moment,
                            const RotationalState& state,
                            ModelStatus& status) const override
    {
        if (config_.inertia <= 0.0)
        {
            status.ok = false;
            status.error = "rigid rotor: inertia must be positive";
            return state;
        }

        if (dt <= 0.0)
        {
            status.ok = false;
            status.error = "dt must be positive";
            return state;
        }

        const double alpha = net_moment / config_.inertia;

        RotationalState next = state;
        next.omega += alpha * dt;

        if (config_.integration == RotationalIntegration::Euler)
        {
            next.angle_rad += state.omega * dt;
        }
        else // Heun (trapezoidal average of omega over the step)
        {
            next.angle_rad += 0.5 * (state.omega + next.omega) * dt;
        }

        return next;
    }

private:
    RigidRotorConfig config_;
};

} // anonymous namespace

// ── Factory functions ─────────────────────────────────────────────

std::unique_ptr<IRotationalDynamics> make_rigid_rotor_dynamics(const RigidRotorConfig& config)
{
    return std::make_unique<RigidRotorDynamics>(config);
}

} // namespace exd::physics::mechanics