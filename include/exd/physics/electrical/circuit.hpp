#pragma once

#include <exd/physics/mechanics/moment_model.hpp>
#include <exd/physics/model_status.hpp>

#include <memory>
#include <string>

namespace exd::physics::electrical {

// ─────────────────────────────────────────────────────
// Lumped electrical machines — the coupling payoff of the
// electrical domain: motors/generators plug into the
// generic mechanics stack.
//
//   DcMotorModel        — dynamic (armature current integrated
//                         with the shared integrator module)
//   make_dc_motor_moment — quasi-steady T(ω) as IMomentModel
//                         (steady-state armature current)
// ─────────────────────────────────────────────────────

struct DcMotorConfig
{
    double kt = 0.0;      // torque constant (N·m/A)
    double ke = 0.0;      // back-emf constant (V·s/rad)
    double R = 0.0;       // armature resistance (Ω, > 0)
    double L = 0.0;       // armature inductance (H, >= 0; 0 = quasi-steady)
    double v_supply = 0.0;// applied armature voltage (V)
};

/// Dynamic DC motor: integrates the armature current
///   L·di/dt + R·i = v_supply − ke·ω
/// where ω comes from the caller (explicit coupling per step). When L == 0
/// the current is evaluated quasi-steadily: i = (v_supply − ke·ω)/R.
class DcMotorModel
{
public:
    explicit DcMotorModel(const DcMotorConfig& config);

    /// Advance the armature current by `dt` at the given shaft speed and
    /// return the electromagnetic torque T = kt·i (N·m; sign follows i:
    /// positive when motoring/assisting, negative when generating/braking).
    double step(double dt, double omega, exd::physics::ModelStatus& status);

    /// Torque from the last step (N·m).
    double torque() const;
    /// Armature current from the last step (A).
    double current() const;
    /// Config sanity (R > 0, L >= 0, kt/ke consistent).
    bool valid() const;
    /// Zero the internal current.
    void reset();

private:
    DcMotorConfig config_;
    double current_ = 0.0;
    double torque_ = 0.0;
    std::string error_;
};

/// Quasi-steady motor as a generic mechanics load (IMomentModel convention:
/// positive = opposing rotation): T(ω) = kt·(ke·ω − v_supply)/R.
/// Positive when generating/braking (ke·ω > v_supply), negative when
/// motoring/assisting (ke·ω < v_supply). Null when R <= 0.
std::unique_ptr<mechanics::IMomentModel> make_dc_motor_moment(const DcMotorConfig& config);

} // namespace exd::physics::electrical
