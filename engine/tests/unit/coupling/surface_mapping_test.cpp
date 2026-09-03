// Surface-mapping tests: nearest-neighbor point transfer (transfer_nearest)
// and trilinear structured-grid sampling (transfer_trilinear).

#include <exd/engine/coupling/surface_mapping.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace exd::engine;
using namespace exd::engine::coupling;

namespace
{

// 5x5x5 node-centered grid, origin (0,0,0), spacing 0.5:
// value = x + y + z (a trilinear field reproduced exactly by interpolation).
StructuredScalarGrid make_linear_grid()
{
    constexpr int32_t nx = 5;
    constexpr int32_t ny = 5;
    constexpr int32_t nz = 5;

    StructuredScalarGrid grid;
    grid.origin = {0.0, 0.0, 0.0};
    grid.spacing = {0.5, 0.5, 0.5};
    grid.dims = {nx, ny, nz};
    grid.values.assign(nx * ny * nz, 0.0);
    for (int32_t k = 0; k < nz; ++k)
        for (int32_t j = 0; j < ny; ++j)
            for (int32_t i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i + nx * (j + ny * k));
                grid.values[idx] = 0.5 * static_cast<double>(i + j + k);
            }
    return grid;
}

} // anonymous namespace

// -- Nearest transfer ---------------------------------------------

TEST_CASE("nearest transfer")
{
    // Corners of a square in the z = 0 plane.
    const std::vector<std::array<double, 3>> corners = {
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0},
        {1.0, 1.0, 0.0},
        {-1.0, 1.0, 0.0},
    };
    const std::vector<double> corner_values = {10.0, 20.0, 30.0, 40.0};

    SUBCASE("target at the center picks the first nearest corner (tie)")
    {
        const std::vector<std::array<double, 3>> target_points = {{0.0, 0.0, 0.0}};
        std::vector<double> target_values(1, -1.0);
        ModelStatus status;
        REQUIRE(transfer_nearest(corner_values, corners, target_points,
                                 target_values, status));
        REQUIRE(status.ok);
        CHECK(target_values[0] == doctest::Approx(10.0));
    }

    SUBCASE("target at an exact source point copies its value")
    {
        const std::vector<std::array<double, 3>> target_points = {corners[2]};
        std::vector<double> target_values(1, -1.0);
        ModelStatus status;
        REQUIRE(transfer_nearest(corner_values, corners, target_points,
                                 target_values, status));
        REQUIRE(status.ok);
        CHECK(target_values[0] == doctest::Approx(30.0));
    }

    SUBCASE("source size mismatch is an error")
    {
        const std::vector<std::array<double, 3>> target_points = {{0.0, 0.0, 0.0}};
        std::vector<double> target_values(1, -1.0);
        std::vector<double> bad_source = {10.0, 20.0, 30.0}; // 3 vs 4 source points
        ModelStatus status;
        CHECK_FALSE(transfer_nearest(bad_source, corners, target_points,
                                     target_values, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }
}

// -- Trilinear grid transfer --------------------------------------

TEST_CASE("trilinear on aligned grids")
{
    const StructuredScalarGrid source = make_linear_grid();

    SUBCASE("half-spacing offset target matches the analytic linear field")
    {
        StructuredScalarGrid target;
        target.origin = {0.25, 0.25, 0.25}; // squarely inside the source box [0,2]^3
        target.spacing = {0.5, 0.5, 0.5};
        target.dims = {3, 3, 3};
        std::vector<double> out(3 * 3 * 3, 0.0);

        ModelStatus status;
        REQUIRE(transfer_trilinear(source, target, out, status));
        REQUIRE(status.ok);
        CHECK(status.warnings.empty());

        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t j = 0; j < 3; ++j)
                for (std::size_t i = 0; i < 3; ++i)
                {
                    const double x = 0.25 + 0.5 * static_cast<double>(i);
                    const double y = 0.25 + 0.5 * static_cast<double>(j);
                    const double z = 0.25 + 0.5 * static_cast<double>(k);
                    const std::size_t idx = i + 3 * (j + 3 * k);
                    // Trilinear reproduces linear fields exactly.
                    CHECK(out[idx] == doctest::Approx(x + y + z).epsilon(1e-9));
                }
    }

    SUBCASE("out-of-bounds target nodes are NaN with a warning")
    {
        StructuredScalarGrid target;
        target.origin = {1.8, 1.8, 1.8}; // nodes 1.8 / 2.0 / 2.2 per axis
        target.spacing = {0.2, 0.2, 0.2};
        target.dims = {3, 3, 3};
        std::vector<double> out(3 * 3 * 3, 0.0);

        ModelStatus status;
        REQUIRE(transfer_trilinear(source, target, out, status));
        REQUIRE(status.ok);
        CHECK_FALSE(status.warnings.empty());

        // Nodes on the 2.2 planes are outside the source box -> NaN.
        std::size_t nan_count = 0;
        for (double value : out)
        {
            if (std::isnan(value))
                ++nan_count;
        }
        CHECK(nan_count > 0);

        // Node (2.0, 2.0, 2.0) is exactly the source corner -> in bounds.
        const std::size_t corner_index = 1u + 3u * (1u + 3u * 1u);
        CHECK(out[corner_index] == doctest::Approx(6.0));
    }
}

TEST_CASE("trilinear reproduces node values")
{
    const StructuredScalarGrid source = make_linear_grid();

    // Target nodes coincide with source nodes on the sub-box
    // [0.5, 1.5]^3 (spacing 0.5, aligned origins).
    StructuredScalarGrid target;
    target.origin = {0.5, 0.5, 0.5};
    target.spacing = {0.5, 0.5, 0.5};
    target.dims = {3, 3, 3};

    std::vector<double> out(3 * 3 * 3, 0.0);
    ModelStatus status;
    REQUIRE(transfer_trilinear(source, target, out, status));
    REQUIRE(status.ok);

    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t i = 0; i < 3; ++i)
            {
                const double x = 0.5 + 0.5 * static_cast<double>(i);
                const double y = 0.5 + 0.5 * static_cast<double>(j);
                const double z = 0.5 + 0.5 * static_cast<double>(k);
                const std::size_t idx = i + 3 * (j + 3 * k);
                CHECK(out[idx] == doctest::Approx(x + y + z));
            }
}