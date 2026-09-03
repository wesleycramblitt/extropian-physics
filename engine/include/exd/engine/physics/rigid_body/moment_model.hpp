#pragma once

#include "rotational_state.hpp"
#include "status.hpp"

#include <memory>
#include <string_view>
#include <vector>

namespace exd::engine::physics::rigid_body {

// ─────────────────────────────────────────────────────
// External (load) moments acting on a rotating body.
// Application data expressed as generic moment laws:
// generator torque, brake, pump head curve, bearing
// friction — all are just T(ω) models.
// ─────────────────────────────────────────────────────

class IMomentModel
{
public:
    virtual ~IMomentModel() = default;
    virtual std::string_view name() const = 0;

    /// Load moment magnitude opposing rotation (positive opposes +omega), N·m.
    virtual double moment(const RotationalState& state, ModelStatus& status) const = 0;
};

/// Constant opposing torque (N·m).
struct ConstantMomentConfig
{
    double torque = 0.0;
};

/// Linear opposing torque: T = k·ω + offset (k in N·m·s/rad, offset in N·m).
struct LinearMomentConfig
{
    double k = 0.0;
    double offset = 0.0;
};

/// Piecewise-linear opposing torque vs ω. `omega_pts` strictly increasing.
struct CurveMomentConfig
{
    std::vector<double> omega_pts;  // rad/s
    std::vector<double> torque_pts; // N·m
};

std::unique_ptr<IMomentModel> make_constant_moment(const ConstantMomentConfig& config);
std::unique_ptr<IMomentModel> make_linear_moment(const LinearMomentConfig& config);
std::unique_ptr<IMomentModel> make_curve_moment(const CurveMomentConfig& config);

} // namespace exd::engine::physics::rigid_body
