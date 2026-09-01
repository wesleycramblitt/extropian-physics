// Phase I thermal solver: steady-state conduction + uniform-flow advection
// + source on a regular grid, SOR relaxation (see header).

#include <exd/physics/thermal/thermal_solver.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace exd::physics::thermal {

namespace {

// Flat index: i + nx*(j + ny*k).
inline std::size_t flat_index(const ThermalGridConfig& grid, int i, int j, int k)
{
    const int nx = grid.dims[0];
    const int ny = grid.dims[1];
    return static_cast<std::size_t>(i + nx * (j + ny * k));
}

// Face index constants (boundary_kind / boundary_values order).
enum : std::size_t {
    FaceXMax = 0,
    FaceXMin = 1,
    FaceYMax = 2,
    FaceYMin = 3,
    FaceZMax = 4,
    FaceZMin = 5,
};

// Core validation; shared by the public entry point and the solver so the
// reported warnings are identical either way.
bool validate_internal(const ThermalConfig& config,
                       std::string& error,
                       std::vector<std::string>& warnings)
{
    const ThermalGridConfig& g = config.grid;
    for (int a = 0; a < 3; ++a)
    {
        if (g.dims[a] < 2)
        {
            error = "thermal: grid dims must be >= 2 per axis";
            return false;
        }
        if (!(g.spacing[a] > 0.0))
        {
            error = "thermal: grid spacing must be > 0 per axis";
            return false;
        }
    }
    if (!(config.material.conductivity > 0.0))
    {
        error = "thermal: conductivity must be > 0";
        return false;
    }
    if (!(config.material.density > 0.0))
    {
        error = "thermal: material density must be > 0";
        return false;
    }
    if (!(config.material.specific_heat > 0.0))
    {
        error = "thermal: specific heat must be > 0";
        return false;
    }
    if (!(config.dt > 0.0))
    {
        error = "thermal: dt must be > 0";
        return false;
    }
    if (config.max_steps == 0)
    {
        error = "thermal: max_steps must be > 0";
        return false;
    }
    if (!(config.tolerance > 0.0))
    {
        error = "thermal: tolerance must be > 0";
        return false;
    }

    // Advection Peclet guard: Pe = rho*cp*|u|*dx_min/(2k).  When Pe crosses 1
    // the first-order upwind advection dominates and adds numerical
    // diffusion; the honest guidance is to refine the grid.
    const double speed = std::sqrt(config.body_velocity[0] * config.body_velocity[0] +
                                   config.body_velocity[1] * config.body_velocity[1] +
                                   config.body_velocity[2] * config.body_velocity[2]);
    if (speed > 0.0)
    {
        const double dx_min = std::min({g.spacing[0], g.spacing[1], g.spacing[2]});
        const double pe = config.material.density * config.material.specific_heat *
                          speed * dx_min / (2.0 * config.material.conductivity);
        if (pe > 1.0)
        {
            warnings.push_back(
                "thermal: advection-dominated with first-order upwind (numerical diffusion; refine grid)");
        }
    }

    // A fully insulated domain is a pure-Neumann problem: the field is only
    // determined up to a constant.  That is not an error, but it is worth
    // telling the caller.
    bool has_fixed = false;
    for (std::size_t f = 0; f < 6; ++f)
        if (config.boundary_kind[f] == ThermalBoundaryKind::FixedValue)
            has_fixed = true;
    if (!has_fixed)
    {
        warnings.push_back(
            "thermal: no fixed-value boundary; temperature is determined only up to a constant");
    }

    return true;
}

// Set fixed (Dirichlet) nodes and return the mean fixed value for the
// initial guess of relaxed nodes.
double mark_fixed_nodes(const ThermalConfig& config, std::vector<double>& T,
                        std::vector<bool>& fixed)
{
    const ThermalGridConfig& g = config.grid;
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    const int nz = g.dims[2];

    double mean = 0.0;
    int fixed_count = 0;

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t I = flat_index(g, i, j, k);
                bool on_fixed = false;
                double value = 0.0;

                // Priority order: +x, -x, +y, -y, +z, -z.  Any fixed face
                // pins the node (see header).
                if (i == nx - 1 && config.boundary_kind[FaceXMax] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceXMax];
                }
                else if (i == 0 && config.boundary_kind[FaceXMin] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceXMin];
                }
                else if (j == ny - 1 && config.boundary_kind[FaceYMax] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceYMax];
                }
                else if (j == 0 && config.boundary_kind[FaceYMin] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceYMin];
                }
                else if (k == nz - 1 && config.boundary_kind[FaceZMax] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceZMax];
                }
                else if (k == 0 && config.boundary_kind[FaceZMin] == ThermalBoundaryKind::FixedValue)
                {
                    on_fixed = true;
                    value = config.boundary_values[FaceZMin];
                }

                fixed[I] = on_fixed;
                if (on_fixed)
                {
                    T[I] = value;
                    mean += value;
                    ++fixed_count;
                }
            }

    return fixed_count > 0 ? mean / static_cast<double>(fixed_count) : 300.0;
}

// Neighbor value across one axis with ghost-cell mirroring for insulated
// faces.  `neighbor_low`/`neighbor_high` return the values at index-1 and
// index+1; when the index is at the domain edge the neighbor is mirrored
// from just inside (zero normal gradient), which is only applied to relaxed
// boundary nodes that are on insulated faces by construction.
inline void axis_neighbors(const ThermalGridConfig& g, const std::vector<double>& T,
                           int i, int j, int k,
                           int axis, double& low_value, double& high_value)
{
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    if (axis == 0)
    {
        low_value = (i == 0) ? T[flat_index(g, 1, j, k)] : T[flat_index(g, i - 1, j, k)];
        high_value = (i == nx - 1) ? T[flat_index(g, nx - 2, j, k)] : T[flat_index(g, i + 1, j, k)];
        return;
    }
    if (axis == 1)
    {
        low_value = (j == 0) ? T[flat_index(g, i, 1, k)] : T[flat_index(g, i, j - 1, k)];
        high_value = (j == ny - 1) ? T[flat_index(g, i, ny - 2, k)] : T[flat_index(g, i, j + 1, k)];
        return;
    }
    low_value = (k == 0) ? T[flat_index(g, i, j, 1)] : T[flat_index(g, i, j, k - 1)];
    high_value = (k == g.dims[2] - 1) ? T[flat_index(g, i, j, g.dims[2] - 2)]
                                      : T[flat_index(g, i, j, k + 1)];
}

} // anonymous namespace

bool validate_thermal_config(const ThermalConfig& config,
                             std::string& error,
                             std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();
    return validate_internal(config, error, warnings);
}

ThermalResult solve_thermal(const ThermalConfig& config, ModelStatus& status)
{
    ThermalResult result;

    std::string error;
    std::vector<std::string> warnings;
    if (!validate_internal(config, error, warnings))
    {
        result.status.ok = false;
        result.status.error = error;
        result.status.warnings = warnings;
        status = result.status;
        return result;   // result.ok stays false
    }

    const ThermalGridConfig& g = config.grid;
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    const int nz = g.dims[2];
    const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;

    const double dx = g.spacing[0];
    const double dy = g.spacing[1];
    const double dz = g.spacing[2];
    const double Dx = 1.0 / (dx * dx);
    const double Dy = 1.0 / (dy * dy);
    const double Dz = 1.0 / (dz * dz);

    const double k = config.material.conductivity;
    const double rho_cp = config.material.density * config.material.specific_heat;
    const double ux = config.body_velocity[0];
    const double uy = config.body_velocity[1];
    const double uz = config.body_velocity[2];
    const double source = config.source_density;
    constexpr double kOmega = 1.5;   // SOR relaxation factor

    std::vector<double> T(n, 0.0);
    std::vector<bool> fixed(n, false);
    const double init_mean = mark_fixed_nodes(config, T, fixed);
    for (std::size_t I = 0; I < n; ++I)
        if (!fixed[I])
            T[I] = init_mean;

    // ---- SOR relaxation -------------------------------------------------
    double max_residual = 0.0;
    uint64_t sweeps = 0;
    for (; sweeps < config.max_steps; ++sweeps)
    {
        max_residual = 0.0;
        for (int kk = 0; kk < nz; ++kk)
            for (int jj = 0; jj < ny; ++jj)
                for (int ii = 0; ii < nx; ++ii)
                {
                    const std::size_t I = flat_index(g, ii, jj, kk);
                    if (fixed[I])
                        continue;

                    double tm_x = 0.0, tp_x = 0.0;
                    double tm_y = 0.0, tp_y = 0.0;
                    double tm_z = 0.0, tp_z = 0.0;
                    axis_neighbors(g, T, ii, jj, kk, 0, tm_x, tp_x);
                    axis_neighbors(g, T, ii, jj, kk, 1, tm_y, tp_y);
                    axis_neighbors(g, T, ii, jj, kk, 2, tm_z, tp_z);

                    // Conducting part of the implicit solve numerator and the
                    // diagonal denominator (central 7-point Laplacian).
                    double numerator = k * ((tm_x + tp_x) * Dx +
                                            (tm_y + tp_y) * Dy +
                                            (tm_z + tp_z) * Dz) +
                                       source;
                    double denominator = 2.0 * k * (Dx + Dy + Dz);

                    // First-order upwind advection: -rho*cp*(u . grad T).
                    // u>0 -> dT/dx ~= (T_i - T_low)/h ; u<0 -> dT/dx ~= (T_high - T_i)/h.
                    // The +rho*cp*|u|/h diagonal keeps the relaxation stable.
                    if (ux > 0.0)
                    {
                        denominator += rho_cp * ux / dx;
                        numerator += rho_cp * ux * tm_x / dx;
                    }
                    else if (ux < 0.0)
                    {
                        denominator -= rho_cp * ux / dx;   // rho*cp*|u|/dx
                        numerator -= rho_cp * ux * tp_x / dx;
                    }
                    if (uy > 0.0)
                    {
                        denominator += rho_cp * uy / dy;
                        numerator += rho_cp * uy * tm_y / dy;
                    }
                    else if (uy < 0.0)
                    {
                        denominator -= rho_cp * uy / dy;
                        numerator -= rho_cp * uy * tp_y / dy;
                    }
                    if (uz > 0.0)
                    {
                        denominator += rho_cp * uz / dz;
                        numerator += rho_cp * uz * tm_z / dz;
                    }
                    else if (uz < 0.0)
                    {
                        denominator -= rho_cp * uz / dz;
                        numerator -= rho_cp * uz * tp_z / dz;
                    }

                    const double t_new = numerator / denominator;
                    const double residual = std::abs(t_new - T[I]);
                    max_residual = std::max(max_residual, residual);
                    T[I] = (1.0 - kOmega) * T[I] + kOmega * t_new;
                }

        if (max_residual <= config.tolerance)
            break;
    }

    // ---- Result bookkeeping ---------------------------------------------
    result.iterations = sweeps;
    result.max_residual = max_residual;

    result.max_temperature = *std::max_element(T.begin(), T.end());
    result.min_temperature = *std::min_element(T.begin(), T.end());

    // total_power = sum source*cell_volume over the (nx-1)(ny-1)(nz-1) cells.
    const double cell_volume = dx * dy * dz;
    result.total_power = source * static_cast<double>(nx - 1) *
                         static_cast<double>(ny - 1) * static_cast<double>(nz - 1) *
                         cell_volume;

    result.temperature.origin = g.origin;
    result.temperature.spacing = g.spacing;
    result.temperature.dims = g.dims;
    result.temperature.values = T;

    result.ok = true;
    result.status.ok = true;
    result.status.error.clear();
    result.status.warnings = warnings;
    status = result.status;
    return result;
}

} // namespace exd::physics::thermal