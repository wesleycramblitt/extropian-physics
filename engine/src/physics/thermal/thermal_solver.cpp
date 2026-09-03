// Phase I thermal solver: steady-state conduction + uniform-flow advection
// + source on a regular grid, SOR relaxation (see header).

#include <exd/engine/physics/thermal/thermal_solver.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace exd::engine::physics::thermal {

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

// Point of node (i,j,k).
inline std::array<double, 3> node_point(const ThermalGridConfig& g, int i, int j, int k)
{
    return {g.origin[0] + g.spacing[0] * static_cast<double>(i),
            g.origin[1] + g.spacing[1] * static_cast<double>(j),
            g.origin[2] + g.spacing[2] * static_cast<double>(k)};
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

// Per-node advecting velocity: channel sample when available (fallback to
// body_velocity on out-of-bounds), otherwise body_velocity.
inline std::array<double, 3> node_velocity(const ThermalConfig& config, int i, int j, int k,
                                           ModelStatus& status)
{
    if (config.velocity_channel == nullptr)
        return config.body_velocity;
    std::array<double, 3> v{};
    bool ok = config.velocity_channel->sample(node_point(config.grid, i, j, k), v);
    if (!ok)
    {
        if (status.warnings.empty() ||
            status.warnings.back() != "thermal: velocity channel out of bounds; using body_velocity")
            status.warnings.push_back(
                "thermal: velocity channel out of bounds; using body_velocity");
        return config.body_velocity;
    }
    return v;
}

/// One SOR relaxation sweep shared by the steady solve and the transient
/// advance.  `T_old` is the PREVIOUS time level: the implicit time term's
/// right-hand side must stay fixed at the old level while the sweep relaxes
/// `T` in place (using T itself would cancel the time term at the fixed
/// point).  `inv_time` = rho*cp/dt for transient (0 for steady).  Fixed
/// nodes hold their values; returns the max |dT| applied in the sweep.
double sor_sweep(const ThermalConfig& config, std::vector<double>& T,
                 const std::vector<double>& T_old,
                 const std::vector<bool>& fixed, double inv_time, ModelStatus& status)
{
    const ThermalGridConfig& g = config.grid;
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    const int nz = g.dims[2];
    const double dx = g.spacing[0];
    const double dy = g.spacing[1];
    const double dz = g.spacing[2];
    const double Dx = 1.0 / (dx * dx);
    const double Dy = 1.0 / (dy * dy);
    const double Dz = 1.0 / (dz * dz);

    const double k = config.material.conductivity;
    const double rho_cp = config.material.density * config.material.specific_heat;
    const bool has_source_channel = config.source_channel != nullptr;
    const double source = config.source_density;
    constexpr double kOmega = 1.5;

    double max_residual = 0.0;
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

                // Conducting part of the implicit numerator and diagonal.
                double qsource = source;
                if (has_source_channel)
                {
                    const std::array<double, 3> p{
                        g.origin[0] + ii * g.spacing[0],
                        g.origin[1] + jj * g.spacing[1],
                        g.origin[2] + kk * g.spacing[2]};
                    config.source_channel->sample(p, qsource);
                }
                double numerator = k * ((tm_x + tp_x) * Dx +
                                        (tm_y + tp_y) * Dy +
                                        (tm_z + tp_z) * Dz) +
                                   qsource;
                double denominator = 2.0 * k * (Dx + Dy + Dz);

                // First-order upwind advection with the per-node velocity.
                const std::array<double, 3> u = node_velocity(config, ii, jj, kk, status);
                if (u[0] > 0.0)
                {
                    denominator += rho_cp * u[0] / dx;
                    numerator += rho_cp * u[0] * tm_x / dx;
                }
                else if (u[0] < 0.0)
                {
                    denominator -= rho_cp * u[0] / dx;
                    numerator -= rho_cp * u[0] * tp_x / dx;
                }
                if (u[1] > 0.0)
                {
                    denominator += rho_cp * u[1] / dy;
                    numerator += rho_cp * u[1] * tm_y / dy;
                }
                else if (u[1] < 0.0)
                {
                    denominator -= rho_cp * u[1] / dy;
                    numerator -= rho_cp * u[1] * tp_y / dy;
                }
                if (u[2] > 0.0)
                {
                    denominator += rho_cp * u[2] / dz;
                    numerator += rho_cp * u[2] * tm_z / dz;
                }
                else if (u[2] < 0.0)
                {
                    denominator -= rho_cp * u[2] / dz;
                    numerator -= rho_cp * u[2] * tp_z / dz;
                }

                // Implicit time term: (rho*cp/dt)*T_new on the diagonal, the
                // PREVIOUS time level on the right-hand side.
                if (inv_time > 0.0)
                {
                    denominator += inv_time;
                    numerator += inv_time * T_old[I];
                }

                const double t_new = numerator / denominator;
                const double residual = std::abs(t_new - T[I]);
                max_residual = std::max(max_residual, residual);
                T[I] = (1.0 - kOmega) * T[I] + kOmega * t_new;
            }
    return max_residual;
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

    std::vector<double> T(n, 0.0);
    std::vector<bool> fixed(n, false);
    const double init_mean = mark_fixed_nodes(config, T, fixed);
    for (std::size_t I = 0; I < n; ++I)
        if (!fixed[I])
            T[I] = init_mean;

    // ---- SOR relaxation (shared sweep; steady: no time term) ----------
    double max_residual = 0.0;
    uint64_t sweeps = 0;
    for (; sweeps < config.max_steps; ++sweeps)
    {
        max_residual = sor_sweep(config, T, T, fixed, 0.0, status);
        if (max_residual <= config.tolerance)
            break;
    }

    // ---- Result bookkeeping ---------------------------------------------
    result.iterations = sweeps;
    result.max_residual = max_residual;

    result.max_temperature = *std::max_element(T.begin(), T.end());
    result.min_temperature = *std::min_element(T.begin(), T.end());

    // total_power = sum source*cell_volume over the (nx-1)(ny-1)(nz-1) cells.
    const bool has_source_channel = config.source_channel != nullptr;
    const double cell_volume = dx * dy * dz;
    if (has_source_channel)
    {
        // integrate the sampled source over the (nx-1)(ny-1)(nz-1) cells,
        // sampled at the cell MIN corner (the same convention the preset's
        // q-grid power integral uses)
        double power = 0.0;
        for (int kk = 0; kk < nz - 1; ++kk)
            for (int jj = 0; jj < ny - 1; ++jj)
                for (int ii = 0; ii < nx - 1; ++ii)
                {
                    const std::array<double, 3> p{
                        g.origin[0] + ii * g.spacing[0],
                        g.origin[1] + jj * g.spacing[1],
                        g.origin[2] + kk * g.spacing[2]};
                    double q = 0.0;
                    config.source_channel->sample(p, q);
                    power += q * cell_volume;
                }
        result.total_power = power;
    }
    else
        result.total_power = config.source_density * static_cast<double>(nx - 1) *
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

// ---------------------------------------------------------------------
// W11: transient stepping + coupling-writable state
// ---------------------------------------------------------------------

bool init_thermal_state(ThermalState& state, const ThermalConfig& config,
                        ModelStatus& status)
{
    std::string error;
    std::vector<std::string> warnings;
    if (!validate_internal(config, error, warnings))
    {
        status.ok = false;
        status.error = error;
        status.warnings = warnings;
        return false;
    }

    const ThermalGridConfig& g = config.grid;
    const std::size_t n = static_cast<std::size_t>(g.dims[0]) * g.dims[1] * g.dims[2];

    std::vector<double> T(n, config.initial_temperature);
    std::vector<bool> fixed(n, false);
    const double init_mean = mark_fixed_nodes(config, T, fixed);
    for (std::size_t I = 0; I < n; ++I)
        if (!fixed[I] && T[I] <= 0.0)
            T[I] = init_mean;   // seed relaxed nodes away from zero when needed

    state.time = 0.0;
    state.temperature.origin = g.origin;
    state.temperature.spacing = g.spacing;
    state.temperature.dims = g.dims;
    state.temperature.values = std::move(T);
    state.fixed = std::move(fixed);

    status.ok = true;
    status.error.clear();
    status.warnings = warnings;
    return true;
}

bool advance_thermal(ThermalState& state, double dt, const ThermalConfig& config,
                     ModelStatus& status)
{
    if (state.fixed.size() != state.temperature.values.size())
    {
        status.ok = false;
        status.error = "thermal: state/fixed size mismatch (uninitialized state?)";
        return false;
    }
    if (!(dt > 0.0))
    {
        status.ok = false;
        status.error = "thermal: dt must be > 0";
        return false;
    }

    double max_residual = 0.0;
    uint64_t sweeps = 0;
    const double inv_time = config.material.density * config.material.specific_heat / dt;
    const std::vector<double> T_old = state.temperature.values; // previous time level
    for (; sweeps < config.max_steps; ++sweeps)
    {
        max_residual = sor_sweep(config, state.temperature.values, T_old,
                                 state.fixed, inv_time, status);
        if (max_residual <= config.tolerance)
            break;
    }
    if (sweeps >= config.max_steps && max_residual > config.tolerance)
    {
        status.ok = false;
        status.error = "thermal: SOR did not converge within max_steps (residual " +
                       std::to_string(max_residual) + ")";
        return false;
    }

    state.time += dt;
    status.ok = true;
    return true;
}

bool set_temperature_point(ThermalState& state, const std::array<double, 3>& point,
                           double value, ModelStatus& status)
{
    // Adapter: ThermalGridConfig mirrors the StructuredScalarGrid layout.
    const ThermalGridConfig g{state.temperature.origin, state.temperature.spacing,
                              state.temperature.dims};
    const int nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    const int i = static_cast<int>(std::lround((point[0] - g.origin[0]) / g.spacing[0]));
    const int j = static_cast<int>(std::lround((point[1] - g.origin[1]) / g.spacing[1]));
    const int k = static_cast<int>(std::lround((point[2] - g.origin[2]) / g.spacing[2]));
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz)
    {
        status.ok = false;
        status.error = "thermal: set_temperature_point out of bounds";
        return false;
    }
    const std::size_t I = flat_index(g, i, j, k);
    state.temperature.values[I] = value;
    state.fixed[I] = true;
    status.ok = true;
    return true;
}

namespace {

/// Trilinear temperature channel over a thermal state (read-only adapter).
class TemperatureChannel final : public exd::engine::coupling::IScalarField3D
{
public:
    explicit TemperatureChannel(const ThermalState& state) : state_(state) {}

    bool sample(const std::array<double, 3>& p, double& value_out) const override
    {
        const ThermalGridConfig g{state_.temperature.origin, state_.temperature.spacing,
                                  state_.temperature.dims};
        const int nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
        const double fx = (p[0] - g.origin[0]) / g.spacing[0];
        const double fy = (p[1] - g.origin[1]) / g.spacing[1];
        const double fz = (p[2] - g.origin[2]) / g.spacing[2];
        if (fx < 0.0 || fx > static_cast<double>(nx - 1) ||
            fy < 0.0 || fy > static_cast<double>(ny - 1) ||
            fz < 0.0 || fz > static_cast<double>(nz - 1))
            return false;

        const int i0 = static_cast<int>(std::floor(fx));
        const int j0 = static_cast<int>(std::floor(fy));
        const int k0 = static_cast<int>(std::floor(fz));
        const int i1 = std::min(i0 + 1, nx - 1);
        const int j1 = std::min(j0 + 1, ny - 1);
        const int k1 = std::min(k0 + 1, nz - 1);
        const double ax = fx - static_cast<double>(i0);
        const double ay = fy - static_cast<double>(j0);
        const double az = fz - static_cast<double>(k0);
        const auto& v = state_.temperature.values;
        auto at = [&](int i, int j, int k) { return v[flat_index(g, i, j, k)]; };
        const double c00 = at(i0, j0, k0) * (1.0 - ax) + at(i1, j0, k0) * ax;
        const double c10 = at(i0, j1, k0) * (1.0 - ax) + at(i1, j1, k0) * ax;
        const double c01 = at(i0, j0, k1) * (1.0 - ax) + at(i1, j0, k1) * ax;
        const double c11 = at(i0, j1, k1) * (1.0 - ax) + at(i1, j1, k1) * ax;
        const double c0 = c00 * (1.0 - ay) + c10 * ay;
        const double c1 = c01 * (1.0 - ay) + c11 * ay;
        value_out = c0 * (1.0 - az) + c1 * az;
        return true;
    }

private:
    const ThermalState& state_;
};

} // anonymous namespace

std::unique_ptr<exd::engine::coupling::IScalarField3D> make_temperature_channel(const ThermalState& state)
{
    return std::make_unique<TemperatureChannel>(state);
}

ThermalResult simulate_thermal(const ThermalConfig& config, ModelStatus& status)
{
    if (!config.transient)
        return solve_thermal(config, status);

    ThermalState state;
    if (!init_thermal_state(state, config, status))
    {
        ThermalResult result;
        result.status = status;
        return result;
    }

    const double dt = config.dt;
    const uint64_t total = static_cast<uint64_t>(
        std::ceil(config.end_time / config.dt));
    if (total > config.max_time_steps)
    {
        status.ok = false;
        status.error = "thermal: transient would exceed max_time_steps";
        ThermalResult result;
        result.status = status;
        return result;
    }

    uint64_t steps = 0;
    for (; steps < total; ++steps)
    {
        if (!advance_thermal(state, dt, config, status))
        {
            ThermalResult result;
            result.status = status;
            return result;
        }
    }

    ThermalResult result;
    result.ok = true;
    result.status.ok = true;
    result.status.warnings = status.warnings;
    result.iterations = steps;
    result.temperature = std::move(state.temperature);
    const auto& v = result.temperature.values;
    result.max_temperature = *std::max_element(v.begin(), v.end());
    result.min_temperature = *std::min_element(v.begin(), v.end());
    const ThermalGridConfig& g = config.grid;
    result.total_power = config.source_density * g.spacing[0] * g.spacing[1] *
                         g.spacing[2] * static_cast<double>(g.dims[0] - 1) *
                         static_cast<double>(g.dims[1] - 1) *
                         static_cast<double>(g.dims[2] - 1);
    status = result.status;
    return result;
}

} // namespace exd::engine::physics::thermal