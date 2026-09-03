#pragma once

#include "dynamics.hpp"
#include "moment_model.hpp"
#include "rotational_state.hpp"
#include "status.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace exd::engine::physics::rigid_body {

// ─────────────────────────────────────────────────────
// Generic rotating machine assembly.
//
// One step:  provider(ω) → aero loads
//            + external moment model  → net moment
//            → rotational dynamics    → new state
//
// The force provider is injected as a callback so this
// module has zero knowledge of fluids, blades, or CFD —
// it is the reusable "lego" for ANY rotating machine.
// ─────────────────────────────────────────────────────

/// Aero-side loads produced by the injected provider.
struct AeroResult
{
    MomentResult moments;                    // torque about axis + axial force
    std::vector<ElementForce3D> per_element; // diagnostics (may be empty)
};

struct AssemblyStepResult
{
    bool ok = false;
    ModelStatus status;
    RotationalState state;          // updated state
    MomentResult aero;              // aero moments from the provider
    double external_moment = 0.0;   // opposing (positive), N·m
    double net_moment = 0.0;        // aero.torque − external_moment, N·m
    double aero_power = 0.0;        // aero.torque · ω_before (W)
    double mechanical_power = 0.0;  // net_moment · ω_before (W)
    std::vector<ElementForce3D> per_element;
};

class RotatingAssembly
{
public:
    /// Produces aero moments for the current rotational state.
    using ForceProvider = std::function<AeroResult(const RotationalState&, ModelStatus&)>;

    RotatingAssembly(ForceProvider provider,
                     std::unique_ptr<IMomentModel> external,
                     std::unique_ptr<IRotationalDynamics> dynamics);

    /// Advance one timestep from `state`.
    AssemblyStepResult step(double dt, const RotationalState& state);

private:
    ForceProvider provider_;
    std::unique_ptr<IMomentModel> external_;
    std::unique_ptr<IRotationalDynamics> dynamics_;
};

} // namespace exd::engine::physics::rigid_body
