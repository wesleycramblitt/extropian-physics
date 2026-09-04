// Field-channel adapters: wrap a structured grid as the coupling
// field interfaces (scalar/vector flow fields).  The interpolation
// combination arithmetic is shared: mesh::interp::trilinear.

#include <exd/engine/coupling/field_channels.hpp>
#include <exd/engine/mesh/interpolation.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace exd::engine::coupling
{

// local alias: the combination arithmetic is shared (mesh::interp::trilinear)
inline double trilinear_scalar(const std::vector<double>& values,
                               std::size_t base, std::size_t stride,
                               int nx, int ny,
                               double fx, double fy, double fz)
{
    return exd::engine::mesh::interp::trilinear(values, base, stride, nx, ny, fx, fy, fz);
}

namespace {

// Map a world-space query onto the local cell fractions and the lower-node
// flat index of the containing cell.  Returns false when the query lies
// outside the grid (any component < 0 or > dims-1).
inline bool cell_fractions(const std::array<double, 3>& p,
                           const std::array<double, 3>& origin,
                           const std::array<double, 3>& spacing,
                           const std::array<int32_t, 3>& dims,
                           double& fx, double& fy, double& fz,
                           std::size_t& node_base)
{
    const int nx = dims[0];
    const int ny = dims[1];
    const int nz = dims[2];

    // Node-space coordinates; clamp-free bounds check — any component outside
    // [0, dims-1] is outside the field.
    double n[3];
    for (int a = 0; a < 3; ++a)
    {
        n[a] = (p[a] - origin[a]) / spacing[a];
        if (n[a] < 0.0 || n[a] > static_cast<double>(dims[a] - 1))
            return false;
    }

    int i = static_cast<int>(n[0]);
    int j = static_cast<int>(n[1]);
    int k = static_cast<int>(n[2]);

    // Exactly on the last node: clamp to the outermost cell so the local
    // fraction is 1 and the sample stays a pure node value.
    i = std::min(i, nx - 2);
    j = std::min(j, ny - 2);
    k = std::min(k, nz - 2);

    fx = n[0] - static_cast<double>(i);
    fy = n[1] - static_cast<double>(j);
    fz = n[2] - static_cast<double>(k);

    node_base = static_cast<std::size_t>(i + nx * (j + ny * k));
    return true;
}

// ── ScalarGridChannel ─────────────────────────────────────────────
// Node-centered scalar field on a regular grid, trilinear sampling.

class ScalarGridChannel final : public IScalarField3D
{
public:
    explicit ScalarGridChannel(StructuredScalarGrid grid)
        : grid_(std::move(grid))
    {
    }

    const char* name() const { return "scalar_grid"; }

    bool sample(const std::array<double, 3>& p, double& value_out) const override
    {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        std::size_t base = 0;
        if (!cell_fractions(p, grid_.origin, grid_.spacing, grid_.dims,
                            fx, fy, fz, base))
            return false;

        value_out = trilinear_scalar(grid_.values, base, 1,
                                     grid_.dims[0], grid_.dims[1], fx, fy, fz);
        return true;
    }

private:
    StructuredScalarGrid grid_;
};

// ── VectorGridChannel ─────────────────────────────────────────────
// Node-centered vector field on a regular grid.  Values are laid out flat as
// 3·nx·ny·nz with node index idx(i,j,k) = i + nx·(j + ny·k) and a
// per-component start offset of 3·idx.

class VectorGridChannel final : public IVectorField3D
{
public:
    explicit VectorGridChannel(StructuredVectorGrid grid)
        : grid_(std::move(grid))
    {
    }

    const char* name() const { return "vector_grid"; }

    bool sample(const std::array<double, 3>& p,
                std::array<double, 3>& value_out) const override
    {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        std::size_t base = 0;
        if (!cell_fractions(p, grid_.origin, grid_.spacing, grid_.dims,
                            fx, fy, fz, base))
            return false;

        for (int c = 0; c < 3; ++c)
            value_out[c] = trilinear_scalar(grid_.values,
                                            3 * base + static_cast<std::size_t>(c),
                                            3, grid_.dims[0], grid_.dims[1],
                                            fx, fy, fz);
        return true;
    }

private:
    StructuredVectorGrid grid_;
};

// ── Fluid field adapters ──────────────────────────────────────────
// A fluid flow field carries velocity + pressure together; the adapters
// expose each component as a standalone generic channel.  The adapters hold
// a pointer to the underlying field — the caller must keep the field alive.

class VelocityFieldAdapter final : public IVectorField3D
{
public:
    explicit VelocityFieldAdapter(const IFlowField3D& field)
        : field_(&field)
    {
    }

    const char* name() const { return "velocity_adapter"; }

    bool sample(const std::array<double, 3>& p,
                std::array<double, 3>& value_out) const override
    {
        double pressure_dummy = 0.0;
        // Velocity flows through; pressure is discarded.
        return field_->sample(p, value_out, pressure_dummy);
    }

private:
    const IFlowField3D* field_;  // non-owning: caller must keep the field alive
};

class PressureFieldAdapter final : public IScalarField3D
{
public:
    explicit PressureFieldAdapter(const IFlowField3D& field)
        : field_(&field)
    {
    }

    const char* name() const { return "pressure_adapter"; }

    bool sample(const std::array<double, 3>& p, double& value_out) const override
    {
        std::array<double, 3> velocity_dummy{0.0, 0.0, 0.0};
        // Pressure flows through; velocity is discarded.
        return field_->sample(p, velocity_dummy, value_out);
    }

private:
    const IFlowField3D* field_;  // non-owning: caller must keep the field alive
};

} // anonymous namespace

// ── Factories ─────────────────────────────────────────────────────

// The factory signatures have no status channel, so an invalid config is
// reported by returning nullptr (never thrown).  A valid config requires
// dims >= 2 per axis, spacing > 0 per axis, and a value array sized exactly
// for the grid.  The grid is copied into the channel.

std::unique_ptr<IScalarField3D> make_scalar_grid_field(const StructuredScalarGrid& grid)
{
    for (int32_t d : grid.dims)
        if (d < 2)
            return nullptr;
    for (double s : grid.spacing)
        if (s <= 0.0)
            return nullptr;

    const int64_t nx = grid.dims[0];
    const int64_t ny = grid.dims[1];
    const int64_t nz = grid.dims[2];
    const std::size_t node_count = static_cast<std::size_t>(nx * ny * nz);

    if (grid.values.size() != node_count)
        return nullptr;

    return std::make_unique<ScalarGridChannel>(grid);
}

std::unique_ptr<IVectorField3D> make_vector_grid_field(const StructuredVectorGrid& grid)
{
    for (int32_t d : grid.dims)
        if (d < 2)
            return nullptr;
    for (double s : grid.spacing)
        if (s <= 0.0)
            return nullptr;

    const int64_t nx = grid.dims[0];
    const int64_t ny = grid.dims[1];
    const int64_t nz = grid.dims[2];
    const std::size_t node_count = static_cast<std::size_t>(nx * ny * nz);

    if (grid.values.size() != 3 * node_count)
        return nullptr;

    return std::make_unique<VectorGridChannel>(grid);
}

std::unique_ptr<IVectorField3D> make_velocity_field_adapter(const IFlowField3D& field)
{
    return std::make_unique<VelocityFieldAdapter>(field);
}

std::unique_ptr<IScalarField3D> make_pressure_field_adapter(const IFlowField3D& field)
{
    return std::make_unique<PressureFieldAdapter>(field);
}

} // namespace exd::engine::coupling