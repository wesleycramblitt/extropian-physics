#include "bem_internal.hpp"

#include <exd/geometry/turbine.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/engine/physics/fluid/reduced_order/bem/bem_solver.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::reduced_order::bem
{

namespace
{

inline bool is_invalid_config(const BEMSolverConfig& cfg, std::string& reason)
{
    if (cfg.element_count < 4)
    {
        reason = "element_count < 4";
        return true;
    }
    if (cfg.k_duct < 0.0 || cfg.k_duct > 1.0)
    {
        reason = "k_duct outside [0,1]";
        return true;
    }
    if (cfg.under_relaxation <= 0.0 || cfg.under_relaxation > 1.0)
    {
        reason = "under_relaxation outside (0,1]";
        return true;
    }
    if (cfg.glauert_threshold <= 0.0 || cfg.glauert_threshold >= 1.0)
    {
        reason = "glauert_threshold outside (0,1)";
        return true;
    }
    if (cfg.wake_expansion < 0.0 || std::isnan(cfg.wake_expansion))
    {
        reason = "wake_expansion negative or NaN";
        return true;
    }
    if (cfg.upstream_extent < 0.0 || cfg.wake_length < 0.0 ||
        cfg.wake_decay_length < 0.0 || cfg.wake_radius_initial < 0.0)
    {
        reason = "negative length sentinel";
        return true;
    }
    return false;
}

} // namespace

TurbineResult solve_turbine(const exd::geometry::TurbineDefinition& turbine,
                            const OperatingConditions& conditions,
                            const PolarDatabase& polars,
                            const BEMSolverConfig& config)
{
    TurbineResult result;
    auto& warnings = result.warnings;

    std::string reason;
    if (is_invalid_config(config, reason))
    {
        result.error = "invalid BEMSolverConfig: " + reason;
        return result;
    }

    if (conditions.v_inf <= 0.0)
    {
        result.error = "v_inf <= 0";
        return result;
    }

    if (config.row_index >= turbine.blade_rows.size())
    {
        result.error = "row_index out of range";
        return result;
    }

    if (polars.empty())
    {
        result.error = "empty polar database";
        return result;
    }

    auto geo_res = build_blade_geometry(turbine, config);
    warnings.insert(warnings.end(), geo_res.warnings.begin(), geo_res.warnings.end());
    if (!geo_res.ok)
    {
        result.error = geo_res.error;
        return result;
    }

    const auto& geo = geo_res.geometry;
    const double rpm = conditions.rpm_override.value_or(geo.rpm);
    if (rpm <= 0.0)
    {
        result.error = "rpm <= 0";
        return result;
    }
    const double omega = rpm * 2.0 * M_PI / 60.0;

    // Duct acceleration.
    auto duct = compute_duct_state(turbine, geo, config, conditions);
    if (!duct.ok)
    {
        result.error = duct.error;
        return result;
    }
    const double V_rotor = duct.v_rotor;

    result.duct = {config.k_duct, duct.area_ratio, duct.m_duct, V_rotor};

    // Create correction models from config.
    auto induction = make_induction_model(config.induction_correction);
    auto loss = make_loss_model(config.loss_correction);

    // Element loop.
    double torque = 0.0;
    double thrust = 0.0;
    bool all_converged = true;

    result.radial.reserve(geo.elements.size());

    for (const auto& elem : geo.elements)
    {
        auto state = induction->solve(elem, V_rotor, omega, conditions, polars,
                                      config, *loss, warnings);
        all_converged = all_converged && state.converged;

        const double sin_phi = std::sin(state.phi_rad);
        const double cos_phi = std::cos(state.phi_rad);
        const double dynamic_pressure = 0.5 * conditions.rho * state.W * state.W;
        const double dL = dynamic_pressure * elem.chord * state.cl * elem.dr;
        const double dD = dynamic_pressure * elem.chord * state.cd * elem.dr;

        const double dQ = elem.blade_count * elem.r * (dL * sin_phi - dD * cos_phi);
        const double dT = elem.blade_count * (dL * cos_phi + dD * sin_phi);
        const double dp = dT / (2.0 * M_PI * elem.r * elem.dr);

        torque += dQ;
        thrust += dT;

        RadialStation rs;
        rs.radius_m = elem.r;
        rs.span = elem.span;
        rs.axial_velocity = state.Va;
        rs.relative_velocity = state.W;
        rs.induction_axial = state.a;
        rs.induction_tangential = state.a_prime;
        rs.inflow_angle_deg = state.phi_rad * 180.0 / M_PI;
        rs.angle_of_attack_deg = state.alpha_deg;
        rs.reynolds = state.Re;
        rs.cl = state.cl;
        rs.cd = state.cd;
        rs.lift_per_m = dL / elem.dr;
        rs.drag_per_m = dD / elem.dr;
        rs.torque_per_m = dQ / elem.dr;
        rs.thrust_per_m = dT / elem.dr;
        rs.pressure_jump = dp;
        rs.converged = state.converged;
        rs.iterations = state.iterations;
        result.radial.push_back(rs);
    }

    const double power = torque * omega;
    const double lambda_tsr = omega * geo.r_tip / conditions.v_inf;

    double A_ref = M_PI * geo.r_tip * geo.r_tip;
    if (config.reference_area == ReferenceArea::Annulus)
        A_ref = M_PI * (geo.r_tip * geo.r_tip - geo.r_hub * geo.r_hub);

    const double denom_power = 0.5 * conditions.rho * A_ref * conditions.v_inf *
                                conditions.v_inf * conditions.v_inf;
    const double denom_force = 0.5 * conditions.rho * A_ref * conditions.v_inf * conditions.v_inf;
    const double cp = (denom_power > 0.0) ? power / denom_power : 0.0;
    const double ct = (denom_force > 0.0) ? thrust / denom_force : 0.0;

    result.rotor = {torque, thrust, power, cp, ct, cp / (16.0 / 27.0)};

    // Hull drag.
    const double r_shroud_max = max_shroud_radius_over_domain(geo.shroud_spline,
        geo.z_r - (config.wake_length > 0.0 ? config.wake_length : 8.0 * geo.r_tip),
        geo.z_r + (config.upstream_extent > 0.0 ? config.upstream_extent : 3.0 * geo.r_tip));

    HullInput hull_input;
    hull_input.r_shroud_max = r_shroud_max;
    hull_input.hull_cd = config.hull_cd;
    hull_input.v_inf = conditions.v_inf;
    hull_input.rho = conditions.rho;
    hull_input.mu = conditions.mu;
    result.hull = compute_hull_forces(geo, hull_input, warnings);

    const double T_net = thrust - result.hull.drag;
    result.system = {T_net, power, result.rotor.efficiency * ((thrust > 0.0) ? (T_net / thrust) : 0.0)};

    if (config.include_flow_field)
    {
        result.flow_field = compute_flow_field(geo, duct, result.radial, conditions, config);
    }
    else
    {
        // Empty grid vectors when the field is disabled.
        result.flow_field = FlowFieldGrid{};
    }

    result.valid = true;
    result.converged = all_converged;
    return result;
}

} // namespace exd::engine::physics::fluid::reduced_order::bem
