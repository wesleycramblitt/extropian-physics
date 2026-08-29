// Field-channel tests: scalar/vector structured-grid trilinear channels and
// the fluid-field adapters that expose velocity / pressure as generic
// channels.

#include <exd/physics/coupling/field_channels.hpp>
#include <exd/physics/coupling/field_sampler.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <vector>

using namespace exd::physics;
using namespace exd::physics::coupling;

namespace
{

// ── Fixtures ──────────────────────────────────────────────────────

// 3×3×3 node-centered grid, origin (0,0,0), spacing 1: values = x + y + z
// with values[i + nx*(j + ny*k)] = x_i + y_j + z_k (nodes at (i, j, k)).
StructuredScalarGrid make_linear_scalar_grid()
{
    constexpr int nx = 3;
    constexpr int ny = 3;
    constexpr int nz = 3;

    StructuredScalarGrid grid;
    grid.origin = {0.0, 0.0, 0.0};
    grid.spacing = {1.0, 1.0, 1.0};
    grid.dims = {nx, ny, nz};
    grid.values.assign(nx * ny * nz, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i + nx * (j + ny * k));
                grid.values[idx] = static_cast<double>(i + j + k);
            }
    return grid;
}

// Same grid but with a vector field [vx, vy, vz] = [x, y, z] per node.
StructuredVectorGrid make_linear_vector_grid()
{
    constexpr int nx = 3;
    constexpr int ny = 3;
    constexpr int nz = 3;

    StructuredVectorGrid grid;
    grid.origin = {0.0, 0.0, 0.0};
    grid.spacing = {1.0, 1.0, 1.0};
    grid.dims = {nx, ny, nz};
    grid.values.assign(3 * nx * ny * nz, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i + nx * (j + ny * k));
                grid.values[3 * idx + 0] = static_cast<double>(i);
                grid.values[3 * idx + 1] = static_cast<double>(j);
                grid.values[3 * idx + 2] = static_cast<double>(k);
            }
    return grid;
}

} // anonymous namespace

// ── Scalar grid channel ───────────────────────────────────────────

TEST_CASE("scalar grid channel reproduces a linear scalar field")
{
    auto field = make_scalar_grid_field(make_linear_scalar_grid());
    REQUIRE(field != nullptr);

    double value = 0.0;

    // Interior point: trilinear interpolation reproduces the linear field.
    CHECK(field->sample({1.5, 0.0, 0.0}, value));
    CHECK(value == doctest::Approx(1.5));

    CHECK(field->sample({0.25, 1.75, 0.5}, value));
    CHECK(value == doctest::Approx(2.5));

    // Exactly on a node: exact node value.
    CHECK(field->sample({2.0, 1.0, 0.0}, value));
    CHECK(value == doctest::Approx(3.0));

    // Exactly on the last node along an axis: pure node value.
    CHECK(field->sample({2.0, 2.0, 2.0}, value));
    CHECK(value == doctest::Approx(6.0));

    // Out of bounds.
    CHECK_FALSE(field->sample({3.1, 0.0, 0.0}, value));
    CHECK_FALSE(field->sample({0.0, -0.1, 0.0}, value));
    CHECK_FALSE(field->sample({0.0, 0.0, 2.5}, value));
}

TEST_CASE("scalar grid factory rejects invalid configs")
{
    CHECK(make_scalar_grid_field(StructuredScalarGrid{}) == nullptr);

    auto grid = make_linear_scalar_grid();
    grid.dims = {1, 3, 3};
    CHECK(make_scalar_grid_field(grid) == nullptr);

    grid = make_linear_scalar_grid();
    grid.spacing = {0.0, 1.0, 1.0};
    CHECK(make_scalar_grid_field(grid) == nullptr);

    grid = make_linear_scalar_grid();
    grid.values.resize(grid.values.size() - 1);
    CHECK(make_scalar_grid_field(grid) == nullptr);
}

// ── Vector grid channel ───────────────────────────────────────────

TEST_CASE("vector grid channel reproduces a linear vector field")
{
    auto field = make_vector_grid_field(make_linear_vector_grid());
    REQUIRE(field != nullptr);

    std::array<double, 3> vec{0.0, 0.0, 0.0};

    // Cell center: trilinear interpolation of [x, y, z] is exact.
    CHECK(field->sample({0.5, 0.5, 0.5}, vec));
    CHECK(vec[0] == doctest::Approx(0.5));
    CHECK(vec[1] == doctest::Approx(0.5));
    CHECK(vec[2] == doctest::Approx(0.5));

    // Interior mixed point.
    CHECK(field->sample({1.75, 0.25, 1.5}, vec));
    CHECK(vec[0] == doctest::Approx(1.75));
    CHECK(vec[1] == doctest::Approx(0.25));
    CHECK(vec[2] == doctest::Approx(1.5));

    // Exactly on a corner node.
    CHECK(field->sample({2.0, 1.0, 0.0}, vec));
    CHECK(vec[0] == doctest::Approx(2.0));
    CHECK(vec[1] == doctest::Approx(1.0));
    CHECK(vec[2] == doctest::Approx(0.0));

    // Out of bounds.
    CHECK_FALSE(field->sample({3.1, 0.0, 0.0}, vec));
    CHECK_FALSE(field->sample({0.0, 0.0, 3.0}, vec));
}

TEST_CASE("vector grid factory rejects invalid configs")
{
    CHECK(make_vector_grid_field(StructuredVectorGrid{}) == nullptr);

    auto grid = make_linear_vector_grid();
    grid.dims = {2, 2, 1};
    CHECK(make_vector_grid_field(grid) == nullptr);

    grid = make_linear_vector_grid();
    grid.values.resize(grid.values.size() - 1);
    CHECK(make_vector_grid_field(grid) == nullptr);

    grid = make_linear_vector_grid();
    grid.spacing = {1.0, -1.0, 1.0};
    CHECK(make_vector_grid_field(grid) == nullptr);
}

// ── Fluid adapters ────────────────────────────────────────────────

TEST_CASE("velocity adapter exposes the velocity of a uniform field")
{
    UniformFieldConfig cfg;
    cfg.velocity = {1.0, 2.0, 3.0};
    cfg.p_ref = 101325.0;
    auto fluid = make_uniform_field(cfg);
    REQUIRE(fluid != nullptr);

    auto adapter = make_velocity_field_adapter(*fluid);
    REQUIRE(adapter != nullptr);

    std::array<double, 3> vel{0.0, 0.0, 0.0};
    // Uniform fields answer every query inside bounds.
    CHECK(adapter->sample({-17.0, 0.5, 900.0}, vel));
    CHECK(vel[0] == doctest::Approx(1.0));
    CHECK(vel[1] == doctest::Approx(2.0));
    CHECK(vel[2] == doctest::Approx(3.0));
}

TEST_CASE("pressure adapter exposes the pressure of a uniform field")
{
    UniformFieldConfig cfg;
    cfg.p_ref = 424242.0;
    auto fluid = make_uniform_field(cfg);
    REQUIRE(fluid != nullptr);

    auto adapter = make_pressure_field_adapter(*fluid);
    REQUIRE(adapter != nullptr);

    double pressure = 0.0;
    CHECK(adapter->sample({0.0, 0.0, 0.0}, pressure));
    CHECK(pressure == doctest::Approx(424242.0));
}