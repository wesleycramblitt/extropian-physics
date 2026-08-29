#pragma once

#include <exd/physics/model_status.hpp>

#include <memory>
#include <string_view>

namespace exd::physics::control {

// ─────────────────────────────────────────────────────
// Controllers: pitch/torque control for turbines, speed
// governors for engines, load regulation for compressors.
// First variant: PI with anti-windup.
// ─────────────────────────────────────────────────────

class IController
{
public:
    virtual ~IController() = default;
    virtual std::string_view name() const = 0;

    /// Compute the control effort from setpoint and measurement.
    /// dt must be positive; returns 0.0 and sets status.error otherwise.
    virtual double update(double setpoint, double measurement, double dt,
                          exd::physics::ModelStatus& status) = 0;

    /// Clear internal state (integral terms).
    virtual void reset() = 0;
};

struct PiControllerConfig
{
    double kp = 0.0;        // proportional gain
    double ki = 0.0;        // integral gain
    double clamp_min = -1e30; // output lower limit
    double clamp_max = 1e30;  // output upper limit
    bool anti_windup = true;  // stop integrating while the output is clamped
};

std::unique_ptr<IController> make_pi_controller(const PiControllerConfig& config);

} // namespace exd::physics::control
