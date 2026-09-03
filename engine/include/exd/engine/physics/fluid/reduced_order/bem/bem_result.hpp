#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::reduced_order::bem
{

struct RadialStation
{
    double radius_m = 0.0;
    double span = 0.0;
    double axial_velocity = 0.0;
    double relative_velocity = 0.0;
    double induction_axial = 0.0;
    double induction_tangential = 0.0;
    double inflow_angle_deg = 0.0;
    double angle_of_attack_deg = 0.0;
    double reynolds = 0.0;
    double cl = 0.0;
    double cd = 0.0;
    double lift_per_m = 0.0;
    double drag_per_m = 0.0;
    double torque_per_m = 0.0;
    double thrust_per_m = 0.0;
    double pressure_jump = 0.0;
    bool   converged = false;
    uint32_t iterations = 0;
};

struct RotorResults
{
    double torque = 0.0;
    double thrust = 0.0;
    double power = 0.0;
    double cp = 0.0;
    double ct = 0.0;
    double efficiency = 0.0;  // = Cp / (16/27) (Betz-relative)
};

struct HullForces
{
    double drag = 0.0;
    double pressure_drag = 0.0;
    double viscous_drag = 0.0;
    double cd = 0.0;
    double reference_area = 0.0;
};

struct DuctState
{
    double k_duct = 0.0;
    double area_ratio = 0.0;
    double m_duct = 0.0;
    double v_rotor = 0.0;
};

struct SystemResults
{
    double net_thrust = 0.0;      // T_rotor - D_hull
    double net_power = 0.0;       // = rotor power (hull drag is a support load)
    double efficiency = 0.0;      // = rotor.efficiency * (T_net/T_rotor)
                                  //   engineering FOM, documented as such
};

struct FlowFieldGrid
{
    // Engineering-approximate axial velocity / pressure field.  Not a CFD
    // solution: total pressure is not conserved downstream and the wake does
    // not recover to freestream within the grid.
    std::vector<double> z;        // axial stations, upstream (high z) -> downstream (low z)
    std::vector<double> r;        // radial stations, 0 -> R_max where
                                  //   R_max = max(R_shroud(z)) over the grid domain (>= R_tip)
    std::vector<double> velocity; // axial flow magnitude, row-major [iz*nr+ir]
    std::vector<double> pressure; // gage vs p_ref, same layout
};

struct TurbineResult
{
    bool valid = false;
    std::string error;
    bool converged = false;
    std::vector<std::string> warnings;
    RotorResults rotor;
    HullForces hull;
    DuctState duct;
    SystemResults system;
    std::vector<RadialStation> radial;
    FlowFieldGrid flow_field;
};

} // namespace exd::engine::physics::fluid::reduced_order::bem
