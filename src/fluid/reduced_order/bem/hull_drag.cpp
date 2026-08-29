#include "bem_internal.hpp"

#include <cmath>
#include <cstddef>

namespace exd::physics::fluid::reduced_order::bem
{

namespace
{

inline double shroud_radius_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    return static_cast<double>(spline.evaluate(static_cast<float>(z)));
}

inline double shroud_derivative_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    return static_cast<double>(spline.derivative(static_cast<float>(z)));
}

} // namespace

HullForces compute_hull_forces(const BladeGeometry& geo,
                               const HullInput& input,
                               std::vector<std::string>& warnings)
{
    HullForces hull{};

    // Frontal (pressure-drag) reference area: maximum shroud radius over the
    // flow-field domain.  The caller supplies R_shroud_max so it can be shared
    // with the flow-field grid sizing.
    const double A_f = M_PI * input.r_shroud_max * input.r_shroud_max;
    const double q_inf = 0.5 * input.rho * input.v_inf * input.v_inf;

    hull.cd = input.hull_cd;
    hull.reference_area = A_f;

    if (input.hull_cd <= 0.0)
    {
        // Drag components already zero.
        return hull;
    }

    const double D_hull = q_inf * input.hull_cd * A_f;

    // Axial span of the shroud spline.
    const double L = static_cast<double>(geo.shroud_spline.max_x() - geo.shroud_spline.min_x());
    const double Re_L = (L > 0.0 && input.mu > 0.0)
                            ? input.rho * input.v_inf * L / input.mu
                            : 0.0;
    const double C_f = (Re_L > 0.0) ? 0.074 * std::pow(Re_L, -0.2) : 0.0;

    // Wetted area: A_wet = 2*pi*integral r*sqrt(1 + r'^2) dz over the shroud.
    // Use a fixed 200-sample trapezoid with the spline derivative().  Guard
    // against derivative spikes by clamping |r'| to a generous maximum; the
    // geometric arc-length factor then remains bounded.
    constexpr std::size_t n_samples = 200;
    constexpr double max_slope = 100.0;
    const double z_s = static_cast<double>(geo.shroud_spline.min_x());
    const double z_e = static_cast<double>(geo.shroud_spline.max_x());

    double A_wet = 0.0;
    if (z_e > z_s)
    {
        const double dz = (z_e - z_s) / static_cast<double>(n_samples);
        for (std::size_t i = 0; i < n_samples; ++i)
        {
            const double z0 = z_s + static_cast<double>(i) * dz;
            const double z1 = z0 + dz;
            const double r0 = shroud_radius_at(geo.shroud_spline, z0);
            const double r1 = shroud_radius_at(geo.shroud_spline, z1);
            const double rp0 = shroud_derivative_at(geo.shroud_spline, z0);
            const double rp1 = shroud_derivative_at(geo.shroud_spline, z1);

            const double slope0 = std::max(-max_slope, std::min(max_slope, rp0));
            const double slope1 = std::max(-max_slope, std::min(max_slope, rp1));

            const double f0 = r0 * std::sqrt(1.0 + slope0 * slope0);
            const double f1 = r1 * std::sqrt(1.0 + slope1 * slope1);
            A_wet += 2.0 * M_PI * 0.5 * (f0 + f1) * dz;
        }
    }

    const double D_viscous = q_inf * C_f * A_wet;
    double D_pressure = D_hull - D_viscous;
    if (D_pressure < 0.0)
    {
        warnings.push_back("pressure drag clamped to 0 (viscous estimate exceeded total hull drag)");
        D_pressure = 0.0;
    }

    hull.drag = D_hull;
    hull.pressure_drag = D_pressure;
    hull.viscous_drag = D_viscous;
    return hull;
}

} // namespace exd::physics::fluid::reduced_order::bem
