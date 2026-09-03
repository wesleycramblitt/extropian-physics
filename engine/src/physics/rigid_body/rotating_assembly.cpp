// Generic rotating-machine assembly: combines an injected aero force
// provider, an external moment model and rotational dynamics into a
// single state-advancing step.

#include <exd/engine/physics/rigid_body/rotating_assembly.hpp>

#include <utility>

namespace exd::engine::physics::rigid_body
{

RotatingAssembly::RotatingAssembly(ForceProvider provider,
                                   std::unique_ptr<IMomentModel> external,
                                   std::unique_ptr<IRotationalDynamics> dynamics)
    : provider_(std::move(provider)),
      external_(std::move(external)),
      dynamics_(std::move(dynamics))
{
}

AssemblyStepResult RotatingAssembly::step(double dt, const RotationalState& state)
{
    AssemblyStepResult result;
    result.state = state;

    // 1. Aero loads for the current state.
    const AeroResult aero = provider_(state, result.status);
    if (!result.status.ok)
    {
        result.ok = false;
        return result;
    }

    // Invalid aero moments (e.g. degenerate axis) mean no usable torque:
    // warn, treat the torque as zero and still run the step.
    const double aero_torque = aero.moments.valid ? aero.moments.torque : 0.0;
    if (!aero.moments.valid)
        result.status.warnings.push_back("aero moments invalid (axis?)");

    // 2. External (opposing) moment.
    const double external_moment = external_->moment(state, result.status);
    if (!result.status.ok)
    {
        result.ok = false;
        return result;
    }

    // 3. Net driving moment and integration.
    const double net_moment = aero_torque - external_moment;
    result.state = dynamics_->advance(dt, net_moment, state, result.status);
    if (!result.status.ok)
    {
        result.ok = false;
        return result;
    }

    // 4. Report derived quantities on the pre-step state.
    const double omega_before = state.omega;

    result.ok = true;
    result.aero = aero.moments;
    result.external_moment = external_moment;
    result.net_moment = net_moment;
    result.aero_power = aero_torque * omega_before;
    result.mechanical_power = net_moment * omega_before;
    result.per_element = aero.per_element;

    return result;
}

} // namespace exd::engine::physics::rigid_body