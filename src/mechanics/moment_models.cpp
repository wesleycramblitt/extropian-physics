// External (load) moment models: constant, linear and piecewise-linear
// T(omega) laws. Concrete classes are hidden in an anonymous namespace and
// reachable only through the make_* factory functions.

#include <exd/physics/mechanics/moment_model.hpp>

#include <memory>
#include <vector>

namespace exd::physics::mechanics
{

namespace
{

// ── ConstantMomentModel ───────────────────────────────────────────
// Constant opposing torque, independent of rotational state.

class ConstantMomentModel final : public IMomentModel
{
public:
    explicit ConstantMomentModel(ConstantMomentConfig config) : config_(config) {}

    std::string_view name() const override { return "constant"; }

    double moment(const RotationalState&, ModelStatus&) const override
    {
        return config_.torque;
    }

private:
    ConstantMomentConfig config_;
};

// ── LinearMomentModel ─────────────────────────────────────────────
// Opposing torque linear in omega: T = k*omega + offset.

class LinearMomentModel final : public IMomentModel
{
public:
    explicit LinearMomentModel(LinearMomentConfig config) : config_(config) {}

    std::string_view name() const override { return "linear"; }

    double moment(const RotationalState& state, ModelStatus&) const override
    {
        return config_.k * state.omega + config_.offset;
    }

private:
    LinearMomentConfig config_;
};

// ── CurveMomentModel ──────────────────────────────────────────────
// Piecewise-linear T(omega) with flat clamps beyond the table ends.
// The table is validated on every evaluation (no-exception doctrine):
// mismatched sizes, a non-increasing omega column and a missing omega
// column are config errors. An empty table is a documented "no load"
// case and returns zero torque without error.

class CurveMomentModel final : public IMomentModel
{
public:
    explicit CurveMomentModel(CurveMomentConfig config) : config_(std::move(config)) {}

    std::string_view name() const override { return "curve"; }

    double moment(const RotationalState& state, ModelStatus& status) const override
    {
        const std::vector<double>& omega = config_.omega_pts;
        const std::vector<double>& torque = config_.torque_pts;

        if (omega.size() != torque.size())
        {
            status.ok = false;
            status.error = "curve moment model: omega_pts and torque_pts must have the same size";
            return 0.0;
        }

        if (omega.empty())
            return 0.0; // documented "no load" model

        for (size_t i = 1; i < omega.size(); ++i)
        {
            if (!(omega[i] > omega[i - 1]))
            {
                status.ok = false;
                status.error = "curve moment model: omega_pts must be strictly increasing";
                return 0.0;
            }
        }

        const double w = state.omega;
        if (w <= omega.front())
            return torque.front();
        if (w >= omega.back())
            return torque.back();

        for (size_t i = 1; i < omega.size(); ++i)
        {
            if (w <= omega[i])
            {
                const double t = (w - omega[i - 1]) / (omega[i] - omega[i - 1]);
                return torque[i - 1] + t * (torque[i] - torque[i - 1]);
            }
        }

        return torque.back(); // unreachable; omega is bracketed above
    }

private:
    CurveMomentConfig config_;
};

} // anonymous namespace

// ── Factory functions ─────────────────────────────────────────────

std::unique_ptr<IMomentModel> make_constant_moment(const ConstantMomentConfig& config)
{
    return std::make_unique<ConstantMomentModel>(config);
}

std::unique_ptr<IMomentModel> make_linear_moment(const LinearMomentConfig& config)
{
    return std::make_unique<LinearMomentModel>(config);
}

std::unique_ptr<IMomentModel> make_curve_moment(const CurveMomentConfig& config)
{
    return std::make_unique<CurveMomentModel>(config);
}

} // namespace exd::physics::mechanics