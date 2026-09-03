// Per-cell momentum sources for fdm3: immersed moving-solid penalty and
// Boussinesq buoyancy (W16).

#include <exd/engine/physics/fluid/fdm3/fdm3_sources.hpp>

#include <cmath>

namespace exd::engine::physics::fluid::fdm3 {

namespace {

bool valid_config(const FDM3Config& c, size_t n,
                  exd::engine::core::ModelStatus& status)
{
    if (c.nx < 2 || c.ny < 2 || c.nz < 2 || !(c.lx > 0.0) || !(c.ly > 0.0) || !(c.lz > 0.0))
    {
        status.ok = false;
        status.error = "fdm3 sources: invalid grid config";
        return false;
    }
    if (n != fdm3_cell_count(c))
    {
        status.ok = false;
        status.error = "fdm3 sources: array size does not match nx·ny·nz";
        return false;
    }
    return true;
}

/// Signed distance to the solid surface: < 0 inside, > 0 outside (m).
double solid_signed_distance(const ImmersedSolid& s, const std::array<double, 3>& c)
{
    switch (s.shape)
    {
    case ImmersedShape::Sphere:
        return std::sqrt((c[0] - s.center[0]) * (c[0] - s.center[0]) +
                         (c[1] - s.center[1]) * (c[1] - s.center[1]) +
                         (c[2] - s.center[2]) * (c[2] - s.center[2])) - s.radius;
    case ImmersedShape::Box:
    {
        const double dx = std::fabs(c[0] - s.center[0]) - s.half_extents[0];
        const double dy = std::fabs(c[1] - s.center[1]) - s.half_extents[1];
        const double dz = std::fabs(c[2] - s.center[2]) - s.half_extents[2];
        return std::max({dx, dy, dz});
    }
    }
    return 1e30;
}

/// Smooth solid fraction f ∈ [0,1]: 1 deep inside, 0 far outside, linear
/// transition over `width` centered on the surface (the smeared-IBM
/// transition that keeps the source resolvable by the pressure projection).
double solid_fraction(const ImmersedSolid& s, const std::array<double, 3>& c,
                      double width)
{
    const double d = solid_signed_distance(s, c);
    if (width <= 0.0) return d < 0.0 ? 1.0 : 0.0;
    return std::clamp(0.5 - d / width, 0.0, 1.0);
}

bool cell_in_solid(const ImmersedSolid& s, const std::array<double, 3>& c)
{
    return solid_signed_distance(s, c) < 0.0;
}

} // namespace

bool immersed_solid_mask(const FDM3Config& config,
                         const std::vector<ImmersedSolid>& solids,
                         std::vector<bool>& mask)
{
    const size_t n = fdm3_cell_count(config);
    mask.assign(n, false);
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i)
            {
                const auto c = fdm3_cell_center(config, i, j, k);
                for (const auto& s : solids)
                {
                    if (solid_signed_distance(s, c) < 0.0)
                    {
                        mask[fdm3_cell_index(config, i, j, k)] = true;
                        break;
                    }
                }
            }
    return true;
}

bool add_immersed_solid_forces(const FDM3Config& config,
                               const std::vector<ImmersedSolid>& solids,
                               std::span<const double> u,
                               std::span<const double> v,
                               std::span<const double> w,
                               std::span<double> fx,
                               std::span<double> fy,
                               std::span<double> fz,
                               exd::engine::core::ModelStatus& status)
{
    const size_t n = fdm3_cell_count(config);
    if (!valid_config(config, n, status)) return false;
    if (u.size() < n || v.size() < n || w.size() < n ||
        fx.size() < n || fy.size() < n || fz.size() < n)
    {
        status.ok = false;
        status.error = "fdm3 sources: velocity/force array size mismatch";
        return false;
    }
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i)
            {
                const size_t id = fdm3_cell_index(config, i, j, k);
                const auto c = fdm3_cell_center(config, i, j, k);
                for (const auto& s : solids)
                {
                    const double f =
                        solid_fraction(s, c, s.blend_width > 0.0
                                            ? s.blend_width
                                            : std::max({config.lx / config.nx,
                                                        config.ly / config.ny,
                                                        config.lz / config.nz}));
                    if (f <= 0.0) continue;
                    fx[id] += -s.penalty * f * (u[id] - s.velocity[0]);
                    fy[id] += -s.penalty * f * (v[id] - s.velocity[1]);
                    fz[id] += -s.penalty * f * (w[id] - s.velocity[2]);
                    break;   // union of regions: first match wins
                }
            }
    return true;
}

bool apply_kinematic_freeze(const FDM3Config& config,
                             const std::vector<ImmersedSolid>& solids,
                             double blend,
                             std::span<double> u,
                             std::span<double> v,
                             std::span<double> w,
                             exd::engine::core::ModelStatus& status)
{
    const size_t n = fdm3_cell_count(config);
    if (!valid_config(config, n, status)) return false;
    if (u.size() < n || v.size() < n || w.size() < n)
    {
        status.ok = false;
        status.error = "fdm3 sources: velocity array size mismatch";
        return false;
    }
    const double a = std::max(0.0, std::min(1.0, blend));
    const double w_default = std::max({config.lx / config.nx,
                                       config.ly / config.ny,
                                       config.lz / config.nz});
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i)
            {
                const size_t id = fdm3_cell_index(config, i, j, k);
                const auto c = fdm3_cell_center(config, i, j, k);
                for (const auto& s : solids)
                {
                    const double f = solid_fraction(s, c,
                                                    s.blend_width > 0.0
                                                        ? s.blend_width
                                                        : w_default);
                    if (f <= 0.0) continue;
                    u[id] = (1.0 - a * f) * u[id] + a * f * s.velocity[0];
                    v[id] = (1.0 - a * f) * v[id] + a * f * s.velocity[1];
                    w[id] = (1.0 - a * f) * w[id] + a * f * s.velocity[2];
                    break;
                }
            }
    return true;
}

bool add_boussinesq_forces(const FDM3Config& config,
                           std::span<const double> temperature,
                           double beta, double t_ref,
                           const std::array<double, 3>& gravity,
                           std::span<double> fx,
                           std::span<double> fy,
                           std::span<double> fz,
                           exd::engine::core::ModelStatus& status)
{
    const size_t n = fdm3_cell_count(config);
    if (!valid_config(config, n, status)) return false;
    if (temperature.size() < n || fx.size() < n || fy.size() < n || fz.size() < n)
    {
        status.ok = false;
        status.error = "fdm3 sources: temperature/force array size mismatch";
        return false;
    }
    for (size_t id = 0; id < n; ++id)
    {
        const double dT = temperature[id] - t_ref;
        fx[id] += -beta * dT * gravity[0];
        fy[id] += -beta * dT * gravity[1];
        fz[id] += -beta * dT * gravity[2];
    }
    return true;
}

bool sample_temperature_field(const FDM3Config& config,
                              const coupling::IScalarField3D& channel,
                              std::vector<double>& temperature,
                              exd::engine::core::ModelStatus& status)
{
    const size_t n = fdm3_cell_count(config);
    temperature.resize(n);
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i)
            {
                double t = 0.0;
                if (!channel.sample(fdm3_cell_center(config, i, j, k), t))
                {
                    status.ok = false;
                    status.error = "fdm3 sources: temperature channel out of "
                                   "bounds at a cell center";
                    return false;
                }
                temperature[fdm3_cell_index(config, i, j, k)] = t;
            }
    return true;
}

} // namespace exd::engine::physics::fluid::fdm3
