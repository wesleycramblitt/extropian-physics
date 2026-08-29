#include "bem_internal.hpp"

#include <exd/geometry/turbine.hpp>

#include <cmath>
#include <limits>

namespace exd::physics::fluid::reduced_order::bem
{

namespace
{

inline double shroud_area_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    const double r = static_cast<double>(spline.evaluate(static_cast<float>(z)));
    return M_PI * r * r;
}

} // namespace

DuctOutput compute_duct_state(const exd::geometry::TurbineDefinition& turbine,
                              const BladeGeometry& geo,
                              const BEMSolverConfig& config,
                              const OperatingConditions& conditions)
{
    DuctOutput out;

    // Upstream reference station: front of the shroud spline, or the inlet
    // station leading edge when the geometry provides one.
    double z_u = static_cast<double>(geo.shroud_spline.min_x());
    if (!turbine.flow_path.inlet_station.empty())
    {
        double z_front = std::numeric_limits<double>::infinity();
        for (const auto& p : turbine.flow_path.inlet_station)
            z_front = std::min(z_front, static_cast<double>(p.x));
        if (std::isfinite(z_front))
            z_u = z_front;
    }

    const double A_u = shroud_area_at(geo.shroud_spline, z_u);
    const double A_r = shroud_area_at(geo.shroud_spline, geo.z_r);

    if (A_r <= 0.0)
    {
        out.error = "rotor duct area <= 0";
        return out;
    }

    const double area_ratio = A_u / A_r;
    const double M_duct = 1.0 + config.k_duct * (area_ratio - 1.0);

    // The shroud spline is float-precision upstream; use a small epsilon margin
    // on the M_duct guard so benign round-off does not reject valid geometries.
    constexpr double M_duct_guard = 0.02;
    if (M_duct <= M_duct_guard)
    {
        out.error = "M_duct <= 0.05 (pathological duct geometry)";
        return out;
    }

    out.ok = true;
    out.area_ratio = area_ratio;
    out.m_duct = M_duct;
    out.v_rotor = M_duct * conditions.v_inf;
    return out;
}

} // namespace exd::physics::fluid::reduced_order::bem
