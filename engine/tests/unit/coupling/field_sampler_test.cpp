// Field-sampler unit tests: uniform field, structured-grid trilinear
// interpolation, 2D FDM adapter, and the surface sampling helpers.

#include <exd/engine/coupling/field_sampler.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <vector>

using namespace exd::engine;
using namespace exd::engine::coupling;

namespace
{

// ── Fixtures ──────────────────────────────────────────────────────

// 3×3×3 grid, origin (0,0,0), spacing 1: u = x, v = y, w = 0, p = 100000.
StructuredGridConfig make_linear_grid()
{
    StructuredGridConfig cfg;
    cfg.dims = {3, 3, 3};
    cfg.origin = {0.0, 0.0, 0.0};
    cfg.spacing = {1.0, 1.0, 1.0};
    cfg.rho = 1.2;
    cfg.mu = 1.8e-5;

    constexpr int nx = 3;
    constexpr int ny = 3;
    constexpr int nz = 3;
    cfg.velocity.assign(3 * nx * ny * nz, 0.0);
    cfg.pressure.assign(nx * ny * nz, 100000.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i + nx * (j + ny * k));
                cfg.velocity[3 * idx + 0] = static_cast<double>(i);
                cfg.velocity[3 * idx + 1] = static_cast<double>(j);
                cfg.velocity[3 * idx + 2] = 0.0;
                cfg.pressure[idx] = 100000.0;
            }
    return cfg;
}

// Same grid with a piecewise cell-wise x-velocity: ~0 in cell [0,1] and ~10
// in cell [1,2]; node values 0/5/10 so interpolation blends the two cells.
StructuredGridConfig make_stepped_grid()
{
    StructuredGridConfig cfg;
    cfg.dims = {3, 3, 3};
    cfg.origin = {0.0, 0.0, 0.0};
    cfg.spacing = {1.0, 1.0, 1.0};

    constexpr int nx = 3;
    constexpr int ny = 3;
    constexpr int nz = 3;
    cfg.velocity.assign(3 * nx * ny * nz, 0.0);
    cfg.pressure.assign(nx * ny * nz, 1000.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i + nx * (j + ny * k));
                cfg.velocity[3 * idx + 0] = 5.0 * static_cast<double>(i);
            }
    return cfg;
}

// 4×4 cell-centered 2D FDM field: x = y = {0, 0.5, 1, 1.5}, u = 2x, v = -y,
// p = 1000·(x + y) at cell centers (row-major index i + j·nx).
exd::engine::physics::fluid::fdm::FDMFieldData make_fdm_field()
{
    exd::engine::physics::fluid::fdm::FDMFieldData f;
    f.nx = 4;
    f.ny = 4;
    f.x = {0.0, 0.5, 1.0, 1.5};
    f.y = {0.0, 0.5, 1.0, 1.5};
    f.u.assign(16, 0.0);
    f.v.assign(16, 0.0);
    f.p.assign(16, 0.0);
    for (int j = 0; j < f.ny; ++j)
        for (int i = 0; i < f.nx; ++i)
        {
            const std::size_t id = f.index(i, j);
            f.u[id] = 2.0 * f.x[i];
            f.v[id] = -1.0 * f.y[j];
            f.p[id] = 1000.0 * (f.x[i] + f.y[j]);
        }
    return f;
}

} // anonymous namespace

// ── Uniform field ─────────────────────────────────────────────────

TEST_CASE("uniform field samples constant velocity and pressure everywhere")
{
    UniformFieldConfig cfg;
    cfg.velocity = {5.0, -1.0, 2.5};
    cfg.rho = 1.0;
    cfg.mu = 2.0e-5;
    cfg.p_ref = 90000.0;

    auto field = make_uniform_field(cfg);
    REQUIRE(field != nullptr);

    std::array<double, 3> v{0.0, 0.0, 0.0};
    double p = 0.0;
    CHECK(field->sample({123.4, -5.6, 7.8}, v, p));
    CHECK(v[0] == doctest::Approx(5.0));
    CHECK(v[1] == doctest::Approx(-1.0));
    CHECK(v[2] == doctest::Approx(2.5));
    CHECK(p == doctest::Approx(90000.0));
    CHECK(field->density() == doctest::Approx(1.0));
    CHECK(field->viscosity() == doctest::Approx(2.0e-5));
}

// ── Structured grid ───────────────────────────────────────────────

TEST_CASE("structured grid trilinear interpolation reproduces a linear field")
{
    auto field = make_structured_grid_field(make_linear_grid());
    REQUIRE(field != nullptr);

    std::array<double, 3> v{0.0, 0.0, 0.0};
    double p = 0.0;

    // Interior point: trilinear reproduces the linear field exactly.
    CHECK(field->sample({1.5, 1.5, 1.5}, v, p));
    CHECK(v[0] == doctest::Approx(1.5));
    CHECK(v[1] == doctest::Approx(1.5));
    CHECK(v[2] == doctest::Approx(0.0));
    CHECK(p == doctest::Approx(100000.0));

    CHECK(field->sample({0.25, 1.75, 0.5}, v, p));
    CHECK(v[0] == doctest::Approx(0.25));
    CHECK(v[1] == doctest::Approx(1.75));
    CHECK(v[2] == doctest::Approx(0.0));

    // Exactly on a node: exact node value.
    CHECK(field->sample({1.0, 2.0, 0.0}, v, p));
    CHECK(v[0] == doctest::Approx(1.0));
    CHECK(v[1] == doctest::Approx(2.0));
    CHECK(v[2] == doctest::Approx(0.0));

    // Out of bounds.
    CHECK_FALSE(field->sample({3.1, 0.0, 0.0}, v, p));
    CHECK_FALSE(field->sample({0.0, -0.1, 0.0}, v, p));
    CHECK_FALSE(field->sample({0.0, 0.0, 2.5}, v, p));

    CHECK(field->density() == doctest::Approx(1.2));
    CHECK(field->viscosity() == doctest::Approx(1.8e-5));
}

TEST_CASE("structured grid interpolates between piecewise cell values")
{
    auto field = make_structured_grid_field(make_stepped_grid());
    REQUIRE(field != nullptr);

    std::array<double, 3> v{0.0, 0.0, 0.0};
    double p = 0.0;

    // Cell centers: interpolation lands between adjacent cell values (0/10).
    CHECK(field->sample({0.5, 0.5, 0.5}, v, p));
    CHECK(v[0] == doctest::Approx(2.5));
    CHECK(v[1] == doctest::Approx(0.0));
    CHECK(v[2] == doctest::Approx(0.0));

    CHECK(field->sample({1.5, 0.5, 0.5}, v, p));
    CHECK(v[0] == doctest::Approx(7.5));
    CHECK(v[1] == doctest::Approx(0.0));

    // Boundary samples on the domain corners succeed.
    CHECK(field->sample({0.0, 0.0, 0.0}, v, p));
    CHECK(v[0] == doctest::Approx(0.0));
    CHECK(field->sample({2.0, 2.0, 2.0}, v, p));
    CHECK(v[0] == doctest::Approx(10.0));
}

TEST_CASE("structured grid factory rejects invalid configs")
{
    CHECK(make_structured_grid_field(StructuredGridConfig{}) == nullptr);

    auto cfg = make_linear_grid();
    cfg.dims = {1, 3, 3};
    CHECK(make_structured_grid_field(cfg) == nullptr);

    cfg = make_linear_grid();
    cfg.spacing = {0.0, 1.0, 1.0};
    CHECK(make_structured_grid_field(cfg) == nullptr);

    cfg = make_linear_grid();
    cfg.velocity.resize(cfg.velocity.size() - 1);
    CHECK(make_structured_grid_field(cfg) == nullptr);

    cfg = make_linear_grid();
    cfg.pressure.clear();
    CHECK(make_structured_grid_field(cfg) == nullptr);
}

// ── FDM adapter ───────────────────────────────────────────────────

TEST_CASE("fdm field adapter bilinearly samples in-plane cell-centered data")
{
    auto field = make_fdm_field_adapter(make_fdm_field(), 1.2, 1.8e-5, 101325.0);
    REQUIRE(field != nullptr);

    std::array<double, 3> v{0.0, 0.0, 0.0};
    double p = 0.0;

    // Query on the x[1] = 0.5 / y[1] = 0.5 node column → exact values.
    CHECK(field->sample({0.5, 0.5, 0.0}, v, p));
    CHECK(v[0] == doctest::Approx(1.0));
    CHECK(v[1] == doctest::Approx(-0.5));
    CHECK(v[2] == doctest::Approx(0.0));
    CHECK(p == doctest::Approx(1000.0));

    // Midway between nodes.
    CHECK(field->sample({0.25, 0.25, 0.0}, v, p));
    CHECK(v[0] == doctest::Approx(0.5));
    CHECK(v[1] == doctest::Approx(-0.25));
    CHECK(v[2] == doctest::Approx(0.0));
    CHECK(p == doctest::Approx(500.0));

    // The z coordinate is ignored entirely.
    CHECK(field->sample({0.5, 0.5, 99.0}, v, p));
    CHECK(v[0] == doctest::Approx(1.0));
    CHECK(v[1] == doctest::Approx(-0.5));
    CHECK(v[2] == doctest::Approx(0.0));

    // Exactly on the last node coordinate of both axes.
    CHECK(field->sample({1.5, 1.5, 0.0}, v, p));
    CHECK(v[0] == doctest::Approx(3.0));
    CHECK(v[1] == doctest::Approx(-1.5));
    CHECK(v[2] == doctest::Approx(0.0));
    CHECK(p == doctest::Approx(3000.0));

    // Out of bounds.
    CHECK_FALSE(field->sample({-0.01, 0.5, 0.0}, v, p));
    CHECK_FALSE(field->sample({1.6, 0.5, 0.0}, v, p));
    CHECK_FALSE(field->sample({0.5, -0.01, 0.0}, v, p));
    CHECK_FALSE(field->sample({0.5, 2.0, 0.0}, v, p));

    CHECK(field->density() == doctest::Approx(1.2));
    CHECK(field->viscosity() == doctest::Approx(1.8e-5));
}

// ── Surface sampling ──────────────────────────────────────────────

TEST_CASE("sample_flow maps a uniform field onto a blade surface")
{
    UniformFieldConfig ucfg;
    ucfg.velocity = {4.0, 3.0, 2.0};
    ucfg.rho = 1.0;
    ucfg.mu = 2.0e-5;
    ucfg.p_ref = 95000.0;
    auto uniform = make_uniform_field(ucfg);
    REQUIRE(uniform != nullptr);

    exd::engine::physics::fluid::forces::BladeSurface surface;
    surface.points = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    surface.normals = {{0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
    surface.areas = {0.25, 0.5};
    surface.element_index = {2, 5};

    const double p_ref = 95000.0;
    auto flow = sample_flow(*uniform, surface, p_ref);
    REQUIRE(flow.valid());
    REQUIRE(flow.points.size() == 2);
    CHECK(flow.normals.size() == 2);
    CHECK(flow.velocity[0][0] == doctest::Approx(4.0));
    CHECK(flow.velocity[0][1] == doctest::Approx(3.0));
    CHECK(flow.velocity[1][2] == doctest::Approx(2.0));
    CHECK(flow.shear_traction[0][0] == doctest::Approx(0.0));
    CHECK(flow.shear_traction[1][2] == doctest::Approx(0.0));
    CHECK(flow.pressure[0] == doctest::Approx(p_ref));
    CHECK(flow.pressure[1] == doctest::Approx(p_ref));
    CHECK(flow.area[0] == doctest::Approx(0.25));
    CHECK(flow.area[1] == doctest::Approx(0.5));
    CHECK(flow.element_index[0] == 2);
    CHECK(flow.element_index[1] == 5);
    CHECK(flow.density == doctest::Approx(1.0));
    CHECK(flow.viscosity == doctest::Approx(2.0e-5));
    CHECK(flow.p_ref == doctest::Approx(p_ref));
}

TEST_CASE("sample_flow substitutes fallback values outside the field extent")
{
    auto field = make_structured_grid_field(make_linear_grid());
    REQUIRE(field != nullptr);

    const std::array<std::array<double, 3>, 2> points = {
        {{1.0, 1.0, 1.0}, {5.0, 5.0, 5.0}}};
    const std::array<std::array<double, 3>, 2> normals = {
        {{0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}}};
    const std::array<double, 2> areas = {0.1, 0.2};
    const std::array<int32_t, 2> element_index = {0, 1};

    auto flow = sample_flow(*field, points, normals, areas, element_index, 98000.0);
    REQUIRE(flow.valid());
    REQUIRE(flow.points.size() == 2);

    // In-bounds point sampled from the field (field pressure 100000).
    CHECK(flow.velocity[0][0] == doctest::Approx(1.0));
    CHECK(flow.velocity[0][1] == doctest::Approx(1.0));
    CHECK(flow.velocity[0][2] == doctest::Approx(0.0));
    CHECK(flow.pressure[0] == doctest::Approx(100000.0));

    // Out-of-bounds point: zero velocity, reference pressure, still present.
    CHECK(flow.velocity[1][0] == doctest::Approx(0.0));
    CHECK(flow.velocity[1][1] == doctest::Approx(0.0));
    CHECK(flow.velocity[1][2] == doctest::Approx(0.0));
    CHECK(flow.pressure[1] == doctest::Approx(98000.0));
    CHECK(flow.shear_traction[1][0] == doctest::Approx(0.0));
    CHECK(flow.shear_traction[1][2] == doctest::Approx(0.0));
    CHECK(flow.area[1] == doctest::Approx(0.2));
    CHECK(flow.element_index[1] == 1);
}

TEST_CASE("sample_flow rejects mismatched array lengths")
{
    auto uniform = make_uniform_field(UniformFieldConfig{});
    REQUIRE(uniform != nullptr);

    const std::array<std::array<double, 3>, 2> points = {{{0, 0, 0}, {1, 1, 1}}};
    const std::array<std::array<double, 3>, 1> normals = {{{0, 0, 1}}};
    const std::array<double, 2> areas = {1.0, 1.0};
    const std::array<int32_t, 2> element_index = {0, 1};

    auto flow = sample_flow(*uniform, points, normals, areas, element_index, 101325.0);
    CHECK_FALSE(flow.valid());
}