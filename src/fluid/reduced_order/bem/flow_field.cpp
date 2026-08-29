#include "bem_internal.hpp"

#include <exd/physics/fluid/reduced_order/bem/bem_result.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace exd::physics::fluid::reduced_order::bem
{

namespace
{

// This file implements an analytical / engineering estimate of the axial
// velocity and pressure field.  It is NOT a CFD pressure solution: total
// pressure is not conserved downstream and the wake does not recover to
// freestream within the grid.

inline double shroud_radius_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    return static_cast<double>(spline.evaluate(static_cast<float>(z)));
}

inline double hub_radius_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    return static_cast<double>(spline.evaluate(static_cast<float>(z)));
}

inline double shroud_area_at(const exd::geometry::MonotoneCubicSpline& spline, double z)
{
    const double r = shroud_radius_at(spline, z);
    return M_PI * r * r;
}



// Piecewise-linear interpolation of a radial quantity defined at element
// midpoints.  Returns zero outside [R_hub, R_tip].  The first/last segments
// connect the hub/tip boundaries to zero so the profile is continuous.
double interpolate_radial(double r,
                          double r_hub, double r_tip,
                          const std::vector<double>& radii,
                          const std::vector<double>& values)
{
    if (r <= r_hub || r >= r_tip || radii.empty())
        return 0.0;

    // Below the first midpoint: linear from (r_hub, 0) to (radii.front(), values.front()).
    if (r <= radii.front())
    {
        const double t = (r - r_hub) / (radii.front() - r_hub);
        return t * values.front();
    }

    // Above the last midpoint: linear from (radii.back(), values.back()) to (r_tip, 0).
    if (r >= radii.back())
    {
        const double t = (r_tip - r) / (r_tip - radii.back());
        return t * values.back();
    }

    // Between two element midpoints.
    auto it = std::lower_bound(radii.begin(), radii.end(), r);
    const std::size_t i = static_cast<std::size_t>(it - radii.begin());
    const std::size_t im1 = (i == 0) ? 0 : i - 1;
    const double r0 = radii[im1];
    const double r1 = radii[i];
    const double f = (r1 - r0 > 1e-12) ? (r - r0) / (r1 - r0) : 0.0;
    return values[im1] + f * (values[i] - values[im1]);
}

} // namespace

// Maximum shroud radius over the grid axial domain [z_lo, z_hi].
// The spline evaluation clamps outside [min_x, max_x], so sampling the
// endpoints and any interior control points is sufficient.
double max_shroud_radius_over_domain(const exd::geometry::MonotoneCubicSpline& spline,
                                     double z_lo, double z_hi)
{
    double r_max = 0.0;
    r_max = std::max(r_max, shroud_radius_at(spline, z_lo));
    r_max = std::max(r_max, shroud_radius_at(spline, z_hi));

    const std::size_t n = spline.point_count();
    for (std::size_t i = 0; i < n; ++i)
    {
        const float t = static_cast<float>(i) / static_cast<float>(std::max<std::size_t>(1, n - 1));
        const float z = spline.min_x() + t * (spline.max_x() - spline.min_x());
        if (z >= static_cast<float>(z_lo) && z <= static_cast<float>(z_hi))
            r_max = std::max(r_max, static_cast<double>(spline.evaluate(z)));
    }
    return r_max;
}

FlowFieldGrid compute_flow_field(const BladeGeometry& geo,
                                 const DuctOutput& duct,
                                 const std::vector<RadialStation>& radial,
                                 const OperatingConditions& conditions,
                                 const BEMSolverConfig& config)
{
    FlowFieldGrid grid;

    const double upstream = config.upstream_extent > 0.0 ? config.upstream_extent : 3.0 * geo.r_tip;
    const double wake_len = config.wake_length > 0.0 ? config.wake_length : 8.0 * geo.r_tip;
    const double wake_decay = config.wake_decay_length > 0.0 ? config.wake_decay_length : 4.0 * geo.r_tip;
    const double wake_radius0 = config.wake_radius_initial > 0.0 ? config.wake_radius_initial : geo.r_tip;
    const double k_w = config.wake_expansion;

    const double z_hi = geo.z_r + upstream;
    const double z_lo = geo.z_r - wake_len;

    const uint32_t nz = config.field_axial_points;
    const uint32_t nr = config.field_radial_points;

    grid.z.resize(nz);
    grid.r.resize(nr);
    grid.velocity.assign(nz * nr, 0.0);
    grid.pressure.assign(nz * nr, 0.0);

    // Axial stations: index 0 is the upstream (high-z) end, index nz-1 is the
    // downstream (low-z) end.  Place the rotor plane z_r exactly on the grid
    // so the disk velocity/pressure jump are resolved without snapping.
    std::size_t i_r = 0;
    if (nz > 1)
    {
        const double frac = (z_hi - geo.z_r) / (z_hi - z_lo);
        i_r = static_cast<std::size_t>(std::round(frac * static_cast<double>(nz - 1)));
        i_r = std::clamp(i_r, std::size_t{1}, static_cast<std::size_t>(nz) - 2);
    }

    for (std::size_t i = 0; i <= i_r; ++i)
    {
        const double t = (i_r > 0) ? static_cast<double>(i) / static_cast<double>(i_r) : 0.0;
        grid.z[i] = z_hi - t * (z_hi - geo.z_r);
    }
    for (std::size_t i = i_r + 1; i < nz; ++i)
    {
        const double t = (nz - 1 > i_r)
                             ? static_cast<double>(i - i_r) / static_cast<double>(nz - 1 - i_r)
                             : 0.0;
        grid.z[i] = geo.z_r - t * (geo.z_r - z_lo);
    }

    const double r_shroud_max = max_shroud_radius_over_domain(geo.shroud_spline, z_lo, z_hi);
    const double r_max = std::max(geo.r_tip, r_shroud_max);

    for (uint32_t j = 0; j < nr; ++j)
    {
        const double t = (nr > 1) ? static_cast<double>(j) / static_cast<double>(nr - 1) : 0.0;
        grid.r[j] = t * r_max;
    }

    // Extract element radial data for interpolation.
    std::vector<double> elem_r, elem_a, elem_dp;
    elem_r.reserve(radial.size());
    elem_a.reserve(radial.size());
    elem_dp.reserve(radial.size());
    for (const auto& rs : radial)
    {
        elem_r.push_back(rs.radius_m);
        elem_a.push_back(rs.induction_axial);
        elem_dp.push_back(rs.pressure_jump);
    }

    const double v_inf = conditions.v_inf;
    const double v_rotor = duct.v_rotor;
    const double A_r = shroud_area_at(geo.shroud_spline, geo.z_r);
    const double k_duct = duct.m_duct > 0.0 ? duct.m_duct : config.k_duct;
    const double shroud_min = static_cast<double>(geo.shroud_spline.min_x());
    const double shroud_max = static_cast<double>(geo.shroud_spline.max_x());
    const double v_floor = 0.01 * v_inf;

    const bool has_hub_spline = geo.hub_spline.point_count() > 0;

    constexpr double z_eps = 1e-9;

    for (uint32_t i = 0; i < nz; ++i)
    {
        const double z = grid.z[i];
        const double r_hub_z = has_hub_spline ? hub_radius_at(geo.hub_spline, z) : 0.0;
        const bool at_disk = std::fabs(z - geo.z_r) <= z_eps;
        const bool upstream_of_disk = !at_disk && (z > geo.z_r + z_eps);
        const bool downstream_of_disk = !at_disk && (z < geo.z_r - z_eps);

        double g = 0.0;
        double R_w = wake_radius0;
        if (downstream_of_disk)
        {
            const double dz = geo.z_r - z; // positive downstream distance
            g = std::exp(-dz / wake_decay);
            R_w = wake_radius0 + k_w * dz;
            if (R_w < 1e-6) R_w = 1e-6;
        }

        for (uint32_t j = 0; j < nr; ++j)
        {
            const double r = grid.r[j];
            const std::size_t idx = static_cast<std::size_t>(i) * nr + j;

            // Solid hub body.
            if (r < r_hub_z)
            {
                grid.velocity[idx] = 0.0;
                grid.pressure[idx] = conditions.p_ref + 0.5 * conditions.rho * (v_inf * v_inf);
                continue;
            }

            const double a = interpolate_radial(r, geo.r_hub, geo.r_tip, elem_r, elem_a);
            double V = v_inf;

            if (upstream_of_disk)
            {
                const bool inside_shroud = (r <= shroud_radius_at(geo.shroud_spline, z));
                const bool within_duct = (z >= shroud_min - z_eps && z <= shroud_max + z_eps);
                if (inside_shroud && within_duct)
                {
                    const double A_z = shroud_area_at(geo.shroud_spline, z);
                    V = v_inf * (1.0 + k_duct * (A_z / A_r - 1.0));
                }
                else
                {
                    V = v_inf;
                }
            }
            else if (at_disk)
            {
                V = v_rotor * (1.0 - a);
            }
            else // downstream_of_disk
            {
                const double eta = (R_w > 1e-9) ? (r / R_w) : 0.0;
                const double gauss = std::exp(-eta * eta);
                double deficit_factor = a * (1.0 + (1.0 - g) * gauss);

                // Explicit overshoot guard: deficit factor must stay in [a, 2a].
                const double max_deficit = 2.0 * a;
                if (deficit_factor < a) deficit_factor = a;
                if (deficit_factor > max_deficit) deficit_factor = max_deficit;

                V = v_rotor * (1.0 - deficit_factor);
            }

            if (V < v_floor) V = v_floor;

            // Bernoulli pressure (gage vs p_ref).
            double p = conditions.p_ref + 0.5 * conditions.rho * (v_inf * v_inf - V * V);

            // Pressure jump across the rotor disk applied strictly downstream.
            if (downstream_of_disk)
            {
                const double dp = interpolate_radial(r, geo.r_hub, geo.r_tip, elem_r, elem_dp);
                p -= dp;
            }

            grid.velocity[idx] = V;
            grid.pressure[idx] = p;
        }
    }

    return grid;
}

} // namespace exd::physics::fluid::reduced_order::bem
