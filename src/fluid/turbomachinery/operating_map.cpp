// operating_map.cpp
// Rectangular (omega x mdot) sweep over a stage stack with surge/choke line
// extraction and bilinear sampling.

#include <exd/physics/fluid/turbomachinery/operating_map.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace exd::physics::fluid::turbomachinery
{

bool validate_map_sweep_config(const MapSweepConfig& config,
                               std::string& error,
                               std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.omega_count < 3)
    {
        error = "operating_map: omega_count must be >= 3";
        return false;
    }
    if (config.mdot_count < 3)
    {
        error = "operating_map: mdot_count must be >= 3";
        return false;
    }
    if (config.omega_min < 0.0)
    {
        error = "operating_map: omega_min must be >= 0";
        return false;
    }
    if (config.mdot_min <= 0.0)
    {
        error = "operating_map: mdot_min must be > 0";
        return false;
    }
    if (config.omega_max < config.omega_min || config.mdot_max < config.mdot_min)
    {
        error = "operating_map: sweep maxima must be >= minima";
        return false;
    }
    return true;
}

namespace
{

std::vector<double> linspace(double a, double b, int n)
{
    std::vector<double> pts;
    pts.reserve(static_cast<size_t>(n));
    const double step = (b - a) / static_cast<double>(n - 1);
    for (int i = 0; i < n; ++i)
        pts.push_back(a + step * static_cast<double>(i));
    return pts;
}

} // anonymous namespace

OperatingMap solve_operating_map(const StageStackConfig& stack_config,
                                 const StageInlet& inlet,
                                 const MapSweepConfig& sweep,
                                 const thermo::IEos& eos,
                                 exd::physics::ModelStatus& status)
{
    OperatingMap map;
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    std::string verror;
    std::vector<std::string> vwarnings;
    if (!validate_map_sweep_config(sweep, verror, vwarnings))
    {
        status.ok = false;
        status.error = verror;
        return map;
    }
    if (!validate_stage_stack_config(stack_config, verror, vwarnings))
    {
        status.ok = false;
        status.error = "operating_map: " + verror;
        return map;
    }
    status.warnings.insert(status.warnings.end(), vwarnings.begin(), vwarnings.end());

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t n = static_cast<size_t>(sweep.omega_count);
    const size_t m = static_cast<size_t>(sweep.mdot_count);

    map.omega_pts = linspace(sweep.omega_min, sweep.omega_max, sweep.omega_count);
    map.mdot_pts = linspace(sweep.mdot_min, sweep.mdot_max, sweep.mdot_count);

    map.pressure_ratio.assign(n, std::vector<double>(m, nan));
    map.temperature_ratio.assign(n, std::vector<double>(m, nan));
    map.torque.assign(n, std::vector<double>(m, nan));
    map.surge_mdot.assign(n, nan);
    map.choke_mdot.assign(n, nan);

    for (size_t i = 0; i < n; ++i)
    {
        const double omega = map.omega_pts[i];
        bool choke_found = false;

        for (size_t j = 0; j < m; ++j)
        {
            const double mdot = map.mdot_pts[j];
            ModelStatus point_status;
            StageStackResult solve =
                solve_stage_stack(stack_config, inlet, omega, mdot, eos, point_status);
            status.warnings.insert(status.warnings.end(),
                                   point_status.warnings.begin(),
                                   point_status.warnings.end());

            bool any_choked = false;
            for (const StageResult& stage : solve.per_stage)
                any_choked = any_choked || stage.choked;

            if (!solve.ok || any_choked)
            {
                // Invalidate the point (leave NaN). Record only the FIRST
                // choked point on the line; choked "clearing" at higher mass
                // flow is a mean-line artifact and is ignored.
                if (any_choked && !choke_found)
                {
                    map.choke_mdot[i] = map.mdot_pts[j];
                    choke_found = true;
                }
                continue;
            }

            map.pressure_ratio[i][j] = solve.total_pi;
            map.temperature_ratio[i][j] = solve.T0_out / inlet.T0;
            map.torque[i][j] = solve.total_torque;
        }

        // Surge line: index of the maximum pi over valid points (slope of
        // pi vs mdot flips + to - there). All-NaN lines stay NaN.
        double best_pi = -std::numeric_limits<double>::infinity();
        size_t best_j = 0;
        bool found = false;
        for (size_t j = 0; j < m; ++j)
        {
            const double pi = map.pressure_ratio[i][j];
            if (!std::isnan(pi) && pi > best_pi)
            {
                best_pi = pi;
                best_j = j;
                found = true;
            }
        }
        if (found)
            map.surge_mdot[i] = map.mdot_pts[best_j];
    }

    map.valid = true;
    return map;
}

MapSample sample_operating_map(const OperatingMap& map,
                               double omega,
                               double mdot,
                               exd::physics::ModelStatus& status)
{
    MapSample sample;
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    if (!map.valid)
    {
        status.ok = false;
        status.error = "operating_map: map is invalid";
        return sample;
    }
    const size_t n = map.omega_pts.size();
    const size_t m = map.mdot_pts.size();
    if (n == 0 || m == 0)
    {
        status.ok = false;
        status.error = "operating_map: empty grid";
        return sample;
    }

    const double om_min = map.omega_pts.front();
    const double om_max = map.omega_pts.back();
    const double md_min = map.mdot_pts.front();
    const double md_max = map.mdot_pts.back();

    // Clamp the request into grid extents.
    double o = omega;
    double md = mdot;
    bool clamped = false;
    if (o < om_min || o > om_max)
    {
        clamped = true;
        o = std::clamp(o, om_min, om_max);
    }
    if (md < md_min || md > md_max)
    {
        clamped = true;
        md = std::clamp(md, md_min, md_max);
    }
    if (clamped)
        status.warnings.push_back("operating_map: map sample clamped to grid extents");

    // Bilinear interpolation over the bracketing cell. Returns false when
    // any of the four corners is not a valid (non-NaN) point.
    auto bilinear = [&](const std::vector<std::vector<double>>& grid,
                        double& value) -> bool
    {
        if (n < 2 || m < 2)
            return false;

        size_t i0 = 0;
        while (i0 + 1 < n && map.omega_pts[i0 + 1] <= o)
            ++i0;
        if (i0 >= n - 1)
            i0 = n - 2;
        const double t = (map.omega_pts[i0 + 1] > map.omega_pts[i0])
                             ? (o - map.omega_pts[i0]) /
                                   (map.omega_pts[i0 + 1] - map.omega_pts[i0])
                             : 0.0;

        size_t j0 = 0;
        while (j0 + 1 < m && map.mdot_pts[j0 + 1] <= md)
            ++j0;
        if (j0 >= m - 1)
            j0 = m - 2;
        const double s = (map.mdot_pts[j0 + 1] > map.mdot_pts[j0])
                             ? (md - map.mdot_pts[j0]) /
                                   (map.mdot_pts[j0 + 1] - map.mdot_pts[j0])
                             : 0.0;

        const double g00 = grid[i0][j0];
        const double g10 = grid[i0 + 1][j0];
        const double g01 = grid[i0][j0 + 1];
        const double g11 = grid[i0 + 1][j0 + 1];
        if (std::isnan(g00) || std::isnan(g10) || std::isnan(g01) || std::isnan(g11))
            return false;

        value = (1.0 - t) * (1.0 - s) * g00 + t * (1.0 - s) * g10 +
                (1.0 - t) * s * g01 + t * s * g11;
        return true;
    };

    // Nearest-valid fallback over the whole grid using normalized distance.
    auto nearest_valid = [&](const std::vector<std::vector<double>>& grid,
                             double& value) -> bool
    {
        const double w_range = (om_max > om_min) ? (om_max - om_min) : 1.0;
        const double m_range = (md_max > md_min) ? (md_max - md_min) : 1.0;
        double best_dist = std::numeric_limits<double>::infinity();
        bool found = false;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                if (std::isnan(grid[i][j]))
                    continue;
                const double d_o = (map.omega_pts[i] - o) / w_range;
                const double d_m = (map.mdot_pts[j] - md) / m_range;
                const double dist = d_o * d_o + d_m * d_m;
                if (dist < best_dist)
                {
                    best_dist = dist;
                    value = grid[i][j];
                    found = true;
                }
            }
        }
        return found;
    };

    double value = 0.0;
    auto sample_grid = [&](const std::vector<std::vector<double>>& grid) -> bool
    {
        if (bilinear(grid, value))
            return true;
        sample.valid_region = false;
        if (nearest_valid(grid, value))
        {
            status.warnings.push_back(
                "operating_map: sample cell invalid; nearest-valid corner fallback used");
            return true;
        }
        return false;
    };

    if (!sample_grid(map.pressure_ratio))
    {
        status.ok = false;
        status.error = "operating_map: no valid points in map to sample";
        return sample;
    }
    sample.pressure_ratio = value;

    if (!sample_grid(map.temperature_ratio))
    {
        status.ok = false;
        status.error = "operating_map: no valid points in map to sample";
        return sample;
    }
    sample.temperature_ratio = value;

    if (!sample_grid(map.torque))
    {
        status.ok = false;
        status.error = "operating_map: no valid points in map to sample";
        return sample;
    }
    sample.torque = value;

    sample.ok = true;
    return sample;
}

} // namespace exd::physics::fluid::turbomachinery