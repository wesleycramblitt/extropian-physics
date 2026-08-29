#pragma once

#include "flow_types.hpp"

#include <exd/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>
#include <exd/physics/mechanics/status.hpp>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace exd::physics::fluid::forces {

/// Airfoil polar database used by coefficient-based evaluators
/// (shared with the reduced-order BEM solver).
using PolarDatabase = exd::physics::fluid::reduced_order::bem::PolarDatabase;

// ─────────────────────────────────────────────────────
// Force evaluation: flow field + blade → 3D forces.
// The coupling point between CFD and mechanics.
// ─────────────────────────────────────────────────────

class IForceEvaluator
{
public:
    virtual ~IForceEvaluator() = default;
    virtual std::string_view name() const = 0;

    /// Compute per-element 3D forces on `blade` from the sampled `flow`.
    ///
    /// @param blade   blade geometry (stations)
    /// @param flow    flow sampled at the blade surface (points grouped by
    ///                element_index); coefficient models use only the first
    ///                point's velocity as the uniform inflow
    /// @param omega   current rotational speed (rad/s) — for relative velocity
    /// @param axis    rotation axis — for relative velocity and position
    /// @param per_element  cleared and filled with one entry per used station
    /// @param status  warnings/errors (never throws)
    virtual void compute(const BladeGeometry& blade,
                         const SurfaceFlow& flow,
                         double omega,
                         const mechanics::RotationAxis& axis,
                         std::vector<mechanics::ElementForce3D>& per_element,
                         mechanics::ModelStatus& status) const = 0;
};

enum class ForceEvaluatorType : uint8_t
{
    PressureIntegration, // integrate p·n + τ over the sampled surface (CFD-coupled)
    MomentumBalance,     // blade-element momentum balance, reduced-order (standalone)
    TableLookup,         // coefficient-based blade-element without momentum coupling
};

struct PressureIntegrationConfig
{
    bool include_shear = true; // add shear_traction contribution
};

struct MomentumBalanceConfig
{
    double under_relaxation = 0.25; // (0,1]
    double tolerance = 1e-5;        // induction convergence tolerance
    int max_iterations = 100;       // induction iteration limit
};

struct TableLookupConfig
{
    int max_iterations = 20; // local-flow fixed-point iterations
    double tolerance = 1e-5;
};

/// Bundle of all force-evaluator options, dispatched by `type`.
struct ForceEvaluatorParams
{
    ForceEvaluatorType type = ForceEvaluatorType::MomentumBalance;
    PressureIntegrationConfig pressure;
    MomentumBalanceConfig momentum;
    TableLookupConfig table;
    const PolarDatabase* polars = nullptr; // for momentum/table models;
                                                // nullptr → built-in polars
};

/// Dispatcher factory (repo convention). Never returns nullptr; failures are
/// reported through `status`.
std::unique_ptr<IForceEvaluator> make_force_evaluator(const ForceEvaluatorParams& params,
                                                      mechanics::ModelStatus& status);

// ─────────────────────────────────────────────────────
// Rotor frame and blade-surface construction (3D).
// ─────────────────────────────────────────────────────

/// Orthonormal frame at azimuth θ: e_z = axis direction; e_r, e_t span the
/// rotation plane (positive θ advances e_r → e_t by the right-hand rule).
struct RotorFrame
{
    std::array<double, 3> e_z = {0, 0, 1};
    std::array<double, 3> e_r = {1, 0, 0};
    std::array<double, 3> e_t = {0, 1, 0};
};

/// Build the frame at azimuth `angle_rad` about `axis`.
/// Sets `status.ok = false` when the axis direction is (near) zero.
RotorFrame make_rotor_frame(double angle_rad, const mechanics::RotationAxis& axis,
                            mechanics::ModelStatus& status);

/// Flat-plate lifting-surface samples of all blades at one azimuth.
struct BladeSurface
{
    std::vector<std::array<double, 3>> points;          // N×3 (m)
    std::vector<std::array<double, 3>> normals;         // N×3 outward unit normals
    std::vector<double> areas;                          // N (m²)
    std::vector<int32_t> element_index;                 // N → station index
};

/// Sample the blade surface: per station, one point per side (upper/lower)
/// at mid-chord, separated by thickness_ratio·chord along the face normal.
/// area = dr·chord/2 per side. All blades are placed at their azimuth:
/// θ_b = angle_rad + b·2π/blade_count.
BladeSurface build_blade_surface(const BladeGeometry& blade, double angle_rad,
                                 const mechanics::RotationAxis& axis);

} // namespace exd::physics::fluid::forces
