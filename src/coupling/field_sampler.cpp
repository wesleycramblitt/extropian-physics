// Flow-field samplers: uniform, structured-grid trilinear, and 2D-FDM
// adapted flow fields, plus the surface-sampling helpers that turn CFD-side
// fields into the SurfaceFlow data consumed by force evaluators.

#include <exd/physics/coupling/field_sampler.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace exd::physics::coupling
{

namespace
{

// ── UniformFieldImpl ──────────────────────────────────────────────
// Constant velocity / pressure field, valid at every query point.

class UniformFieldImpl final : public IFlowField3D
{
public:
    explicit UniformFieldImpl(const UniformFieldConfig& config)
        : velocity_(config.velocity),
          rho_(config.rho),
          mu_(config.mu),
          p_ref_(config.p_ref)
    {
    }

    bool sample(const std::array<double, 3>& /*p*/,
                std::array<double, 3>& velocity_out,
                double& pressure_out) const override
    {
        velocity_out = velocity_;
        pressure_out = p_ref_;
        return true;
    }

    double density() const override { return rho_; }
    double viscosity() const override { return mu_; }

private:
    std::array<double, 3> velocity_;
    double rho_;
    double mu_;
    double p_ref_;
};

// ── StructuredGridFieldImpl ───────────────────────────────────────
// Regular grid with trilinear node interpolation.  Velocity is laid out flat
// as 3·nx·ny·nz with node index (i,j,k) = i + nx·(j + ny·k) and a
// per-component start offset; pressure is nx·ny·nz with the same node index.

// Interpolate a scalar field at local cell fractions (fx, fy, fz) inside the
// cell whose lower node is (i, j, k).  `base` is the flat offset of the lower
// node and `stride` is the distance between consecutive node values (1 for
// pressure, 3 for velocity components).
inline double trilinear_scalar(const std::vector<double>& values,
                               std::size_t base, std::size_t stride,
                               int nx, int ny,
                               double fx, double fy, double fz)
{
    const std::size_t sx = stride;                                            // +i neighbor
    const std::size_t sy = stride * static_cast<std::size_t>(nx);             // +j neighbor
    const std::size_t sz = stride * static_cast<std::size_t>(nx) *
                           static_cast<std::size_t>(ny);                      // +k neighbor

    const std::size_t i000 = base;
    const std::size_t i100 = base + sx;
    const std::size_t i010 = base + sy;
    const std::size_t i110 = base + sx + sy;
    const std::size_t i001 = base + sz;
    const std::size_t i101 = base + sx + sz;
    const std::size_t i011 = base + sy + sz;
    const std::size_t i111 = base + sx + sy + sz;

    const double c00 = values[i000] * (1.0 - fx) + values[i100] * fx;
    const double c10 = values[i010] * (1.0 - fx) + values[i110] * fx;
    const double c01 = values[i001] * (1.0 - fx) + values[i101] * fx;
    const double c11 = values[i011] * (1.0 - fx) + values[i111] * fx;

    const double c0 = c00 * (1.0 - fy) + c10 * fy;
    const double c1 = c01 * (1.0 - fy) + c11 * fy;

    return c0 * (1.0 - fz) + c1 * fz;
}

class StructuredGridFieldImpl final : public IFlowField3D
{
public:
    explicit StructuredGridFieldImpl(const StructuredGridConfig& config)
        : config_(config),
          nx_(config.dims[0]),
          ny_(config.dims[1]),
          nz_(config.dims[2])
    {
    }

    bool sample(const std::array<double, 3>& p,
                std::array<double, 3>& velocity_out,
                double& pressure_out) const override
    {
        // Node-space coordinates; clamp-free bounds check — any component
        // outside [0, dims-1] is outside the field.
        double n[3];
        for (int a = 0; a < 3; ++a)
        {
            n[a] = (p[a] - config_.origin[a]) / config_.spacing[a];
            if (n[a] < 0.0 || n[a] > static_cast<double>(config_.dims[a] - 1))
                return false;
        }

        int i = static_cast<int>(n[0]);
        int j = static_cast<int>(n[1]);
        int k = static_cast<int>(n[2]);

        // Exactly on the last node: clamp to the outermost cell so the local
        // fraction is 1 and the sample stays a pure node value.
        i = std::min(i, nx_ - 2);
        j = std::min(j, ny_ - 2);
        k = std::min(k, nz_ - 2);

        const double fx = n[0] - static_cast<double>(i);
        const double fy = n[1] - static_cast<double>(j);
        const double fz = n[2] - static_cast<double>(k);

        const std::size_t node_base =
            static_cast<std::size_t>(i + nx_ * (j + ny_ * k));

        for (int c = 0; c < 3; ++c)
            velocity_out[c] = trilinear_scalar(config_.velocity,
                                               3 * node_base + static_cast<std::size_t>(c),
                                               3, nx_, ny_, fx, fy, fz);
        pressure_out = trilinear_scalar(config_.pressure, node_base, 1, nx_, ny_,
                                        fx, fy, fz);
        return true;
    }

    double density() const override { return config_.rho; }
    double viscosity() const override { return config_.mu; }

private:
    StructuredGridConfig config_;
    int nx_;
    int ny_;
    int nz_;
};

// ── FdmFieldAdapter ───────────────────────────────────────────────
// Collapses a 2D FDM field onto the z-plane: bilinear interpolation in the
// x–y plane; the z coordinate of the query point is ignored (the 2D field has
// no out-of-plane extent).  Queries outside the cell-center extent
// [x.front(), x.back()] × [y.front(), y.back()] report out of bounds.

class FdmFieldAdapter final : public IFlowField3D
{
public:
    FdmFieldAdapter(const fluid::fdm::FDMFieldData& field,
                    double rho, double mu, double p_ref)
        : field_(field),
          rho_(rho),
          mu_(mu),
          p_ref_(p_ref)
    {
    }

    bool sample(const std::array<double, 3>& p,
                std::array<double, 3>& velocity_out,
                double& pressure_out) const override
    {
        const std::size_t nx = field_.x.size();
        const std::size_t ny = field_.y.size();
        if (nx < 2 || ny < 2)
            return false;

        const double px = p[0];
        const double py = p[1];
        if (px < field_.x.front() || px > field_.x.back() ||
            py < field_.y.front() || py > field_.y.back())
            return false;

        // Bracket i so x[i] <= px <= x[i+1].  lower_bound yields the first
        // node >= px, i.e. the cell right endpoint; exactly on the last node
        // (pos == nx-1) we clamp to the last cell, giving a fraction of 1.
        auto ix = std::lower_bound(field_.x.begin(), field_.x.end(), px);
        std::size_t i = static_cast<std::size_t>(ix - field_.x.begin());
        i = std::min(i, nx - 2);
        auto iy = std::lower_bound(field_.y.begin(), field_.y.end(), py);
        std::size_t j = static_cast<std::size_t>(iy - field_.y.begin());
        j = std::min(j, ny - 2);

        const double dx = field_.x[i + 1] - field_.x[i];
        const double dy = field_.y[j + 1] - field_.y[j];
        if (dx <= 0.0 || dy <= 0.0)
            return false;

        const double fx = (px - field_.x[i]) / dx;
        const double fy = (py - field_.y[j]) / dy;

        const std::size_t id00 = field_.index(static_cast<int>(i), static_cast<int>(j));
        const std::size_t id10 = field_.index(static_cast<int>(i + 1), static_cast<int>(j));
        const std::size_t id01 = field_.index(static_cast<int>(i), static_cast<int>(j + 1));
        const std::size_t id11 = field_.index(static_cast<int>(i + 1), static_cast<int>(j + 1));

        const double w00 = (1.0 - fx) * (1.0 - fy);
        const double w10 = fx * (1.0 - fy);
        const double w01 = (1.0 - fx) * fy;
        const double w11 = fx * fy;

        velocity_out[0] = w00 * field_.u[id00] + w10 * field_.u[id10] +
                          w01 * field_.u[id01] + w11 * field_.u[id11];
        velocity_out[1] = w00 * field_.v[id00] + w10 * field_.v[id10] +
                          w01 * field_.v[id01] + w11 * field_.v[id11];
        velocity_out[2] = 0.0; // 2D field: no out-of-plane velocity.
        pressure_out = w00 * field_.p[id00] + w10 * field_.p[id10] +
                       w01 * field_.p[id01] + w11 * field_.p[id11];
        return true;
    }

    double density() const override { return rho_; }
    double viscosity() const override { return mu_; }

private:
    fluid::fdm::FDMFieldData field_;
    double rho_;
    double mu_;
    double p_ref_;
};

} // namespace

// ── Factory functions ─────────────────────────────────────────────

std::unique_ptr<IFlowField3D> make_uniform_field(const UniformFieldConfig& config)
{
    return std::make_unique<UniformFieldImpl>(config);
}

// The factory signature has no status channel, so an invalid config is
// reported by returning nullptr (never thrown); callers must check the
// returned pointer.  A valid config requires dims >= 2 per axis, spacing > 0
// per axis, and velocity / pressure arrays sized exactly for the grid.
std::unique_ptr<IFlowField3D> make_structured_grid_field(const StructuredGridConfig& config)
{
    for (int32_t d : config.dims)
        if (d < 2)
            return nullptr;
    for (double s : config.spacing)
        if (s <= 0.0)
            return nullptr;

    const int64_t nx = config.dims[0];
    const int64_t ny = config.dims[1];
    const int64_t nz = config.dims[2];
    const std::size_t node_count = static_cast<std::size_t>(nx * ny * nz);

    if (config.velocity.size() != 3 * node_count)
        return nullptr;
    if (config.pressure.size() != node_count)
        return nullptr;

    return std::make_unique<StructuredGridFieldImpl>(config);
}

std::unique_ptr<IFlowField3D> make_fdm_field_adapter(const fluid::fdm::FDMFieldData& field,
                                                     double rho, double mu, double p_ref)
{
    return std::make_unique<FdmFieldAdapter>(field, rho, mu, p_ref);
}

// ── Surface sampling ──────────────────────────────────────────────

fluid::forces::SurfaceFlow sample_flow(const IFlowField3D& field,
                                       std::span<const std::array<double, 3>> points,
                                       std::span<const std::array<double, 3>> normals,
                                       std::span<const double> areas,
                                       std::span<const int32_t> element_index,
                                       double p_ref)
{
    // Convenience sampler for parallel user arrays; a length mismatch is not
    // an error here — return an invalid SurfaceFlow that the caller rejects.
    // An empty point set is likewise rejected (valid() requires n > 0).
    const std::size_t n = points.size();
    if (n == 0 || normals.size() != n || areas.size() != n ||
        element_index.size() != n)
        return {};

    fluid::forces::SurfaceFlow flow;
    flow.points.reserve(n);
    flow.normals.reserve(n);
    flow.velocity.reserve(n);
    flow.shear_traction.reserve(n);
    flow.pressure.reserve(n);
    flow.area.reserve(n);
    flow.element_index.reserve(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        std::array<double, 3> velocity{0.0, 0.0, 0.0};
        double pressure = p_ref;
        // Out-of-bounds samples are not fatal: the point is kept with zero
        // velocity at p_ref instead of rejecting the whole surface.
        field.sample(points[i], velocity, pressure);

        flow.points.push_back(points[i]);
        flow.normals.push_back(normals[i]);
        flow.velocity.push_back(velocity);
        flow.shear_traction.push_back({0.0, 0.0, 0.0}); // no wall model
        flow.pressure.push_back(pressure);
        flow.area.push_back(areas[i]);
        flow.element_index.push_back(element_index[i]);
    }

    flow.density = field.density();
    flow.viscosity = field.viscosity();
    flow.p_ref = p_ref;
    return flow;
}

fluid::forces::SurfaceFlow sample_flow(const IFlowField3D& field,
                                       const fluid::forces::BladeSurface& surface,
                                       double p_ref)
{
    return sample_flow(field, surface.points, surface.normals, surface.areas,
                       surface.element_index, p_ref);
}

} // namespace exd::physics::coupling