// Phase I acoustics solver: explicit leapfrog scalar wave equation on a
// regular grid with pressure-release (soft) walls and box-mode verification
// support (see header).

#include <exd/physics/acoustics/wave_solver.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace exd::physics::acoustics {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Flat index: i + nx*(j + ny*k).
inline std::size_t flat_index(const WaveGridConfig& grid, int i, int j, int k)
{
    const int nx = grid.dims[0];
    const int ny = grid.dims[1];
    return static_cast<std::size_t>(i + nx * (j + ny * k));
}

bool validate_internal(const WaveConfig& config,
                       std::string& error,
                       std::vector<std::string>& warnings)
{
    const WaveGridConfig& g = config.grid;
    for (int a = 0; a < 3; ++a)
    {
        if (g.dims[a] < 2)
        {
            error = "acoustics: grid dims must be >= 2 per axis";
            return false;
        }
        if (!(g.spacing[a] > 0.0))
        {
            error = "acoustics: grid spacing must be > 0 per axis";
            return false;
        }
    }
    if (!(config.sound_speed > 0.0))
    {
        error = "acoustics: sound speed must be > 0";
        return false;
    }
    if (config.dt < 0.0)
    {
        error = "acoustics: dt must be >= 0 (0 selects the CFL-adaptive step)";
        return false;
    }
    if (config.max_steps == 0)
    {
        error = "acoustics: max_steps must be > 0";
        return false;
    }
    const int64_t n = static_cast<int64_t>(g.dims[0]) * g.dims[1] * g.dims[2];
    if (config.probe_index < 0 || static_cast<int64_t>(config.probe_index) >= n)
    {
        error = "acoustics: probe_index out of range";
        return false;
    }

    // A requested dt above the exact von Neumann bound is stable only after
    // clamping; the solver clamps deterministically and reports the warning
    // here too so validation and the solve agree.
    if (config.dt > 0.0)
    {
        const double dx = g.spacing[0];
        const double dy = g.spacing[1];
        const double dz = g.spacing[2];
        const double cfl_limit = 1.0 /
            (config.sound_speed * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
        if (config.dt > cfl_limit)
            warnings.push_back("acoustics: CFL violation: unstable");
    }
    return true;
}

} // anonymous namespace

bool validate_wave_config(const WaveConfig& config,
                          std::string& error,
                          std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();
    return validate_internal(config, error, warnings);
}

WaveResult solve_wave(const WaveConfig& config, ModelStatus& status)
{
    WaveResult result;

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

    const WaveGridConfig& g = config.grid;
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    const int nz = g.dims[2];
    const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;

    const double dx = g.spacing[0];
    const double dy = g.spacing[1];
    const double dz = g.spacing[2];
    const double c = config.sound_speed;

    // ---- Time step ------------------------------------------------------
    const double dx_min = std::min({dx, dy, dz});
    const double adaptive_dt = 0.8 * dx_min / (c * std::sqrt(3.0));
    // Exact von Neumann bound: c*dt*sqrt(1/dx^2 + 1/dy^2 + 1/dz^2) <= 1.
    const double cfl_limit = 1.0 / (c * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));

    double dt_used = 0.0;
    if (config.dt <= 0.0)
    {
        dt_used = adaptive_dt;
    }
    else if (config.dt > cfl_limit)
    {
        dt_used = cfl_limit;   // validation already appended the warning
    }
    else
    {
        dt_used = config.dt;
    }
    result.dt_used = dt_used;

    // ---- Initial condition: box-mode seed -------------------------------
    // L_a = (dims_a - 1)*spacing_a.  Direction: dims==2 -> Neumann axis with
    // constant factor 1; otherwise Dirichlet axis with sin(mode*pi*x/L).
    const double Lx = dx * static_cast<double>(nx - 1);
    const double Ly = dy * static_cast<double>(ny - 1);
    const double Lz = dz * static_cast<double>(nz - 1);
    const double kx = (nx > 2) ? static_cast<double>(config.initial_mode[0]) * kPi / Lx : 0.0;
    const double ky = (ny > 2) ? static_cast<double>(config.initial_mode[1]) * kPi / Ly : 0.0;
    const double kz = (nz > 2) ? static_cast<double>(config.initial_mode[2]) * kPi / Lz : 0.0;

    std::vector<double> p_prev(n, 0.0);   // p at t - dt
    std::vector<double> p_cur(n, 0.0);    // p at t
    std::vector<double> p_next(n, 0.0);

    const double cell_volume = dx * dy * dz;
    const double A = config.amplitude;
    std::vector<bool> fixed(n, false);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t I = flat_index(g, i, j, k);
                const double x = dx * static_cast<double>(i);
                const double y = dy * static_cast<double>(j);
                const double z = dz * static_cast<double>(k);
                const double fx = (nx > 2) ? std::sin(kx * x) : 1.0;
                const double fy = (ny > 2) ? std::sin(ky * y) : 1.0;
                const double fz = (nz > 2) ? std::sin(kz * z) : 1.0;
                const double p0 = A * fx * fy * fz;
                p_prev[I] = p0;
                p_cur[I] = p0;   // zero initial velocity -> p(-dt) = p(0)

                // Dirichlet (pinned 0) walls on axes with 3+ nodes.
                fixed[I] = (nx > 2 && (i == 0 || i == nx - 1)) ||
                           (ny > 2 && (j == 0 || j == ny - 1)) ||
                           (nz > 2 && (k == 0 || k == nz - 1));
            }

    double energy_initial = 0.0;
    for (double v : p_cur)
        energy_initial += v * v;
    energy_initial *= cell_volume;
    result.energy_initial = energy_initial;

    double max_pressure = 0.0;
    for (double v : p_cur)
        max_pressure = std::max(max_pressure, std::abs(v));

    result.probe_history.reserve(config.max_steps + 1);
    result.probe_history.push_back(p_cur[static_cast<std::size_t>(config.probe_index)]);

    // ---- Leapfrog time integration --------------------------------------
    const double cdt2 = c * c * dt_used * dt_used;
    const double inv_dx2 = 1.0 / (dx * dx);
    const double inv_dy2 = 1.0 / (dy * dy);
    const double inv_dz2 = 1.0 / (dz * dz);

    for (uint64_t step = 0; step < config.max_steps; ++step)
    {
        std::fill(p_next.begin(), p_next.end(), 0.0);

        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const std::size_t I = flat_index(g, i, j, k);
                    if (fixed[I])
                        continue;

                    // Central second differences with Neumann mirrors on
                    // 2-node axes (neighbor at the far plane mirrors back).
                    double lap = 0.0;
                    if (nx > 2)
                    {
                        const std::size_t I_lo = flat_index(g, i - 1, j, k);
                        const std::size_t I_hi = flat_index(g, i + 1, j, k);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dx2;
                    }
                    else
                    {
                        const int im = (i == 0) ? 1 : 0;
                        const int ip = (i == 1) ? 0 : 1;
                        const std::size_t I_lo = flat_index(g, im, j, k);
                        const std::size_t I_hi = flat_index(g, ip, j, k);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dx2;
                    }
                    if (ny > 2)
                    {
                        const std::size_t I_lo = flat_index(g, i, j - 1, k);
                        const std::size_t I_hi = flat_index(g, i, j + 1, k);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dy2;
                    }
                    else
                    {
                        const int jm = (j == 0) ? 1 : 0;
                        const int jp = (j == 1) ? 0 : 1;
                        const std::size_t I_lo = flat_index(g, i, jm, k);
                        const std::size_t I_hi = flat_index(g, i, jp, k);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dy2;
                    }
                    if (nz > 2)
                    {
                        const std::size_t I_lo = flat_index(g, i, j, k - 1);
                        const std::size_t I_hi = flat_index(g, i, j, k + 1);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dz2;
                    }
                    else
                    {
                        const int km = (k == 0) ? 1 : 0;
                        const int kp = (k == 1) ? 0 : 1;
                        const std::size_t I_lo = flat_index(g, i, j, km);
                        const std::size_t I_hi = flat_index(g, i, j, kp);
                        lap += (p_cur[I_hi] - 2.0 * p_cur[I] + p_cur[I_lo]) * inv_dz2;
                    }

                    p_next[I] = 2.0 * p_cur[I] - p_prev[I] + cdt2 * lap;
                }

        p_prev.swap(p_cur);
        p_cur.swap(p_next);

        for (double v : p_cur)
            max_pressure = std::max(max_pressure, std::abs(v));
        result.probe_history.push_back(p_cur[static_cast<std::size_t>(config.probe_index)]);
    }

    // ---- Result bookkeeping ---------------------------------------------
    result.max_pressure = max_pressure;

    double energy_final = 0.0;
    for (double v : p_cur)
        energy_final += v * v;
    energy_final *= cell_volume;
    result.energy_final = energy_final;

    result.pressure_final.origin = g.origin;
    result.pressure_final.spacing = g.spacing;
    result.pressure_final.dims = g.dims;
    result.pressure_final.values = p_cur;

    result.ok = true;
    result.status.ok = true;
    result.status.error.clear();
    result.status.warnings = warnings;
    status = result.status;
    return result;
}

} // namespace exd::physics::acoustics