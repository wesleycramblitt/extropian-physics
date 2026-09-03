// Lumped electrical machines (Phase D part 3): a dynamic DC motor model
// whose armature current follows L·di/dt + R·i = v_supply − ke·ω, plus the
// quasi-steady equivalent exposed as a generic mechanics load moment
// (exd::engine::physics::rigid_body::IMomentModel). Concrete model classes are hidden in an
// anonymous namespace and reachable only through DcMotorModel /
// make_dc_motor_moment.

#include <exd/engine/physics/electromagnetics/circuit.hpp>

#include <cmath>
#include <memory>
#include <string_view>

namespace exd::engine::physics::electromagnetics
{

// ── Dynamic DC motor ───────────────────────────────────────────────────
// Integrates the armature current over [t, t+dt] with the exact solution of
//   L·di/dt + R·i = v_supply − ke·ω
// at constant shaft speed: steady current i_inf = (v_supply − ke·ω)/R and
//   i(t+dt) = i(t) + (i_inf − i(t))·(1 − exp(−R·dt/L)).
// L == 0 selects the quasi-steady path i = i_inf (no energy stored in the
// armature). The stored `error_` member (declared in the public contract)
// mirrors the status.error message on failure for diagnostics.

DcMotorModel::DcMotorModel(const DcMotorConfig& config)
    : config_(config),
      current_(0.0),
      torque_(0.0)
{
}

bool DcMotorModel::valid() const
{
    // ke may be zero: a pure torque actuator with no back-emf.
    return config_.R > 0.0 && config_.L >= 0.0 && config_.kt > 0.0 && config_.ke >= 0.0;
}

double DcMotorModel::step(double dt, double omega, exd::engine::core::ModelStatus& status)
{
    if (dt <= 0.0)
    {
        status.ok = false;
        status.error = "DC motor: dt must be positive";
        error_ = status.error;
        return torque_;
    }

    if (!valid())
    {
        status.ok = false;
        status.error = "DC motor: invalid config (R must be positive, L >= 0, kt > 0)";
        error_ = status.error;
        return 0.0;
    }

    status.ok = true;
    status.error.clear();

    const double i_inf = (config_.v_supply - config_.ke * omega) / config_.R;
    if (config_.L == 0.0)
    {
        current_ = i_inf;
    }
    else
    {
        const double decay = std::exp(-config_.R * dt / config_.L);
        current_ = current_ * decay + i_inf * (1.0 - decay);
    }

    torque_ = config_.kt * current_;
    return torque_;
}

double DcMotorModel::torque() const
{
    return torque_;
}

double DcMotorModel::current() const
{
    return current_;
}

void DcMotorModel::reset()
{
    current_ = 0.0;
    torque_ = 0.0;
}

// ── Quasi-steady motor as a mechanics load ─────────────────────────────
// T(ω) = kt·(ke·ω − v_supply)/R per the IMomentModel convention (positive =
// opposes rotation): positive when generating/braking (ke·ω > v_supply),
// negative when motoring/assisting. In a RotatingAssembly the net moment is
// aero − external, so an assisting (negative) motor torque adds to the drive.

namespace
{

class DcMotorMomentModel final : public exd::engine::physics::rigid_body::IMomentModel
{
public:
    explicit DcMotorMomentModel(DcMotorConfig config) : config_(config) {}

    std::string_view name() const override { return "dc_motor"; }

    double moment(const exd::engine::physics::rigid_body::RotationalState& state,
                  exd::engine::physics::rigid_body::ModelStatus& status) const override
    {
        if (config_.R <= 0.0)
        {
            status.ok = false;
            status.error = "DC motor: R must be positive";
            return 0.0;
        }
        return config_.kt * (config_.ke * state.omega - config_.v_supply) / config_.R;
    }

private:
    DcMotorConfig config_;
};

} // anonymous namespace

// ── Factory functions ──────────────────────────────────────────────────

std::unique_ptr<exd::engine::physics::rigid_body::IMomentModel> make_dc_motor_moment(const DcMotorConfig& config)
{
    if (config.R <= 0.0)
        return nullptr;
    return std::make_unique<DcMotorMomentModel>(config);
}

} // namespace exd::engine::physics::electromagnetics