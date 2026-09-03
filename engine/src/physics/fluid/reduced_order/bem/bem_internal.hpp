#pragma once

#include <exd/geometry/spline.hpp>
#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_result.hpp>

#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::reduced_order::bem
{

struct BladeElementInput
{
    double r = 0.0;
    double dr = 0.0;
    double chord = 0.0;
    double beta_deg = 0.0;
    double span = 0.0;
    double blade_count = 0.0;
    double r_hub = 0.0;
    double r_tip = 0.0;
    std::string airfoil;
};

struct BladeGeometry
{
    double z_r = 0.0;
    double r_hub = 0.0;
    double r_tip = 0.0;
    double blade_count = 0.0;
    double rpm = 0.0;
    double omega = 0.0;
    std::vector<BladeElementInput> elements;
    exd::geometry::MonotoneCubicSpline shroud_spline;
    exd::geometry::MonotoneCubicSpline hub_spline;
};

struct GeometryResult
{
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    BladeGeometry geometry;
};

GeometryResult build_blade_geometry(const exd::geometry::TurbineDefinition& turbine,
                                    const BEMSolverConfig& config);

struct StationState
{
    double a = 0.0;
    double a_prime = 0.0;
    double phi_rad = 0.0;
    double alpha_deg = 0.0;
    double W = 0.0;
    double Va = 0.0;
    double cl = 0.0;
    double cd = 0.0;
    double Re = 0.0;
    double F = 1.0;
    bool converged = false;
    uint32_t iterations = 0;
};

/// Tip/hub loss correction model interface.
/// Implementations multiply the blade element force by F in [0,1] to account
/// for finite blade count and tip/hub vortex effects.
class LossModel
{
public:
    virtual ~LossModel() = default;
    virtual double loss_factor(double sin_phi, double blade_count,
                               double radius, double r_tip, double r_hub,
                               double axial_induction = 0.0,
                               double chord = 0.0,
                               double sigma = 0.0) const = 0;
};

/// Axial induction solver interface.
/// Implementations iterate the coupled momentum/blade-element equations to
/// find the axial (a) and tangential (a') induction factors at a radial station.
class InductionModel
{
public:
    virtual ~InductionModel() = default;
    virtual StationState solve(const BladeElementInput& elem,
                               double v_rotor, double omega,
                               const OperatingConditions& conditions,
                               const PolarDatabase& polars,
                               const BEMSolverConfig& config,
                               const LossModel& loss_model,
                               std::vector<std::string>& warnings) const = 0;
};

// ── Factory functions ──────────────────────────────────────────────
// Loss model factories.
std::unique_ptr<LossModel> make_prandtl_loss_model();
std::unique_ptr<LossModel> make_du_selig_loss_model();
std::unique_ptr<LossModel> make_chaviaropoulos_loss_model();

std::unique_ptr<LossModel> make_loss_model(LossCorrection type);

// Induction model factories.
std::unique_ptr<InductionModel> make_standard_induction_model();
std::unique_ptr<InductionModel> make_glauert_iterative_induction_model();
std::unique_ptr<InductionModel> make_snel_induction_model();

std::unique_ptr<InductionModel> make_induction_model(InductionCorrection type);

// ── Per-subsystem interfaces (used by bem_solver.cpp) ─────────────

// Duct acceleration model (§6.1).
struct DuctInput
{
    double z_u = 0.0;
    double z_r = 0.0;
    double k_duct = 0.0;
    double v_inf = 0.0;
};

struct DuctOutput
{
    bool ok = false;
    std::string error;
    double area_ratio = 1.0;
    double m_duct = 1.0;
    double v_rotor = 0.0;
};

DuctOutput compute_duct_state(const exd::geometry::TurbineDefinition& turbine,
                              const BladeGeometry& geo,
                              const BEMSolverConfig& config,
                              const OperatingConditions& conditions);

// Maximum shroud radius over an axial domain [z_lo, z_hi].
double max_shroud_radius_over_domain(const exd::geometry::MonotoneCubicSpline& spline,
                                     double z_lo, double z_hi);

// Hull drag model (§8).
struct HullInput
{
    double r_shroud_max = 0.0;
    double hull_cd = 0.0;
    double v_inf = 0.0;
    double rho = 0.0;
    double mu = 0.0;
};

HullForces compute_hull_forces(const BladeGeometry& geo,
                               const HullInput& input,
                               std::vector<std::string>& warnings);

// Engineering-approximate axial velocity / pressure field (§9).
FlowFieldGrid compute_flow_field(const BladeGeometry& geo,
                                 const DuctOutput& duct,
                                 const std::vector<RadialStation>& radial,
                                 const OperatingConditions& conditions,
                                 const BEMSolverConfig& config);

} // namespace exd::engine::physics::fluid::reduced_order::bem
