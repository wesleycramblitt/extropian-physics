// PI controller with output clamping and optional anti-windup. The concrete
// model is hidden in an anonymous namespace and reachable only through the
// make_pi_controller factory function.

#include <exd/engine/physics/control/controller.hpp>

#include <memory>

namespace exd::engine::physics::control
{

namespace
{

// ── PiController ──────────────────────────────────────────────
// out = kp·err + integral, clamped to [clamp_min, clamp_max]. With
// anti-windup enabled the saturating integration step is not stored, so the
// integral never accumulates while the output is saturated and the controller
// keeps full authority once the error reverses.

class PiController final : public IController
{
public:
    explicit PiController(PiControllerConfig config)
        : kp_(config.kp), ki_(config.ki),
          clamp_min_(config.clamp_min), clamp_max_(config.clamp_max),
          anti_windup_(config.anti_windup), integral_(0.0)
    {
    }

    std::string_view name() const override { return "pi"; }

    double update(double setpoint, double measurement, double dt,
                  ModelStatus& status) override
    {
        if (dt <= 0.0)
        {
            status.ok = false;
            status.error = "PI controller: dt must be positive";
            return 0.0;
        }

        const double err = setpoint - measurement;
        const double p = kp_ * err;
        const double integral_step = ki_ * err * dt;
        integral_ += integral_step;
        double out = p + integral_;

        const bool clamped = (out < clamp_min_) || (out > clamp_max_);
        if (out < clamp_min_)
            out = clamp_min_;
        if (out > clamp_max_)
            out = clamp_max_;

        if (clamped && anti_windup_)
            integral_ -= integral_step; // do not accumulate while saturated

        return out;
    }

    void reset() override { integral_ = 0.0; }

private:
    double kp_;
    double ki_;
    double clamp_min_;
    double clamp_max_;
    bool anti_windup_;
    double integral_;
};

} // anonymous namespace

// ── Factory ───────────────────────────────────────────────────

std::unique_ptr<IController> make_pi_controller(const PiControllerConfig& config)
{
    // clamp_min <= clamp_max is intentionally NOT validated here: the caller
    // is responsible for supplying a sane clamping range. The factory is
    // deliberately unconditional so the clamping policy stays composable.
    return std::make_unique<PiController>(config);
}

} // namespace exd::engine::physics::control