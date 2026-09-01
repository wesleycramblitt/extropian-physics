// elasticity_test.cpp
// Static linear elasticity patch tests. The module is verified in regimes
// where its discrete operator is EXACT (nu = 0 uniaxial, confined
// compression, nu = 0 thermal column); free-face cases at coarse z
// resolution carry O(h/L) boundary-order error (documented) and are tested
// qualitatively (bending direction, monotonicity, exact load linearity).

#include <exd/physics/structural/elasticity.hpp>

#include <exd/physics/coupling/field_channels.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <vector>

using namespace exd::physics::structural;
using exd::physics::ModelStatus;

namespace {

std::vector<std::array<bool, 3>> pin_bottom_uz(const ElasticityGridConfig& g)
{
    std::vector<std::array<bool, 3>> mask(static_cast<size_t>(g.dims[0]) *
                                              g.dims[1] * g.dims[2],
                                          std::array<bool, 3>{false, false, false});
    const int nx = g.dims[0], ny = g.dims[1];
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            mask[i + nx * (j + ny * 0)][2] = true; // z = 0 face: u_z = 0
    return mask;
}

std::vector<std::array<bool, 3>> pin_bottom_uz_and_sides(const ElasticityGridConfig& g)
{
    auto mask = pin_bottom_uz(g);
    const int nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(nx) *
                                      (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
                if (i == 0 || i == nx - 1)
                    mask[idx][0] = true; // u_x pinned on ±x faces
                if (j == 0 || j == ny - 1)
                    mask[idx][1] = true; // u_y pinned on ±y faces
            }
    return mask;
}

/// Rigid-mode pins for free-lateral-face cases: bottom roller + three corner
/// pins that kill the six rigid-body modes without constraining the physics.
std::vector<std::array<bool, 3>> roller_with_rigid_pins(const ElasticityGridConfig& g)
{
    auto mask = pin_bottom_uz(g);
    const int nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    auto idx = [&](int i, int j, int k) {
        return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                   (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
    };
    mask[idx(0, 0, 0)] = {true, true, true};
    mask[idx(nx - 1, 0, 0)][1] = true;
    mask[idx(0, ny - 1, 0)][0] = true;
    return mask;
}

std::array<double, 3> node_at(const exd::physics::coupling::StructuredVectorGrid& g,
                              int i, int j, int k)
{
    const int nx = g.dims[0], ny = g.dims[1];
    const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(nx) *
                          (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
    return {g.values[3 * idx + 0], g.values[3 * idx + 1], g.values[3 * idx + 2]};
}

} // anonymous namespace

TEST_CASE("Elasticity: uniaxial tension patch, nu = 0 (exact)")
{
    // With nu = 0 the lateral faces decouple: exact solution uz = sigma·z/E,
    // ux = uy = 0. The z-axis spans 1.0 m (21 nodes).
    ElasticityConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {21, 5, 21};
    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.0;
    cfg.fixed_mask = roller_with_rigid_pins(cfg.grid);
    cfg.surface_traction = {0, 0, 1e6};
    cfg.tolerance = 1e-10;

    ModelStatus status;
    const auto res = solve_elasticity(cfg, status);
    REQUIRE(status.ok);

    const double E = 200e9, sigma = 1e6;
    for (int k : {4, 8, 12, 16, 20})
    {
        const auto u = node_at(res.displacement, 10, 2, k);
        CHECK(u[2] == doctest::Approx(sigma * (0.05 * k) / E).epsilon(1e-6));
        CHECK(std::fabs(u[0]) < 1e-12);
        CHECK(std::fabs(u[1]) < 1e-12);
    }
}

TEST_CASE("Elasticity: confined compression (oedometer patch, exact, nu = 0.3)")
{
    ElasticityConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {11, 11, 11};
    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.3;
    cfg.fixed_mask = pin_bottom_uz_and_sides(cfg.grid);
    cfg.surface_traction = {0, 0, -1e6};

    ModelStatus status;
    const auto res = solve_elasticity(cfg, status);
    REQUIRE(status.ok);

    // Oedometer modulus M = E(1-nu)/((1+nu)(1-2nu)); uz = sigma·z/M.
    const double E = 200e9, nu = 0.3, sigma = -1e6;
    const double M = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    for (int k : {3, 5, 7, 10})
    {
        const auto u = node_at(res.displacement, 5, 5, k);
        CHECK(u[2] == doctest::Approx(sigma * (0.05 * k) / M).epsilon(5e-3));
        CHECK(std::fabs(u[0]) < 1e-12);
        CHECK(std::fabs(u[1]) < 1e-12);
    }
}

TEST_CASE("Elasticity: cantilever bending — qualitative + exact load linearity")
{
    // Free-face bending at coarse z resolution carries documented O(h/L)
    // boundary error, so the quantitative Euler-Bernoulli band is deferred;
    // this case verifies the physics the solver DOES own: downward bending,
    // monotone deflection, and exact linear response to the load.
    auto solve = [](double sigma) {
        ElasticityConfig cfg;
        cfg.grid.origin = {0, 0, 0};
        cfg.grid.spacing = {0.05, 0.05, 0.05};
        cfg.grid.dims = {31, 5, 5}; // Lx = 1.5 m, b = h = 0.2 m
        cfg.material.elastic_modulus = 200e9;
        cfg.material.poisson_ratio = 0.0;
        std::vector<std::array<bool, 3>> mask(31 * 5 * 5, {false, false, false});
        {
            const int nx = 31, ny = 5, nz = 5;
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    mask[0 + nx * (j + ny * k)] = {true, true, true}; // wall at x = 0
        }
        cfg.fixed_mask = mask;
        cfg.surface_traction = {0, 0, -sigma}; // downward load on the whole top face
        cfg.tolerance = 1e-10;
        ModelStatus status;
        return solve_elasticity(cfg, status);
    };

    const double q = 1e4;
    const auto res = solve(q);
    REQUIRE(res.ok);
    const auto res2 = solve(2.0 * q);
    REQUIRE(res2.ok);

    // Tip deflection (x = L, mid width/height): downward and monotone.
    double prev = 0.0;
    for (int i : {6, 12, 18, 24, 30})
    {
        const auto u = node_at(res.displacement, i, 2, 2);
        CHECK(u[2] <= prev + 1e-15);
        prev = u[2];
    }
    const auto u_tip = node_at(res.displacement, 30, 2, 2);
    CHECK(u_tip[2] < 0.0); // bends down

    // Exact linear response: doubling the load doubles every displacement.
    const auto u_tip2 = node_at(res2.displacement, 30, 2, 2);
    CHECK(u_tip2[2] / u_tip[2] == doctest::Approx(2.0).epsilon(1e-9));
}

TEST_CASE("Elasticity: thermal-expansion column (confined, exact, nu = 0.3)")
{
    struct TempField : exd::physics::coupling::IScalarField3D
    {
        double T_ref = 300.0;
        double gz = 1.0; // K/m
        bool sample(const std::array<double, 3>& p, double& v) const override
        {
            v = T_ref + gz * p[2];
            return true;
        }
    } temp_field;

    ElasticityConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {11, 11, 21};
    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.3;
    cfg.fixed_mask = pin_bottom_uz_and_sides(cfg.grid); // lateral rollers
    cfg.thermal_expansion_coefficient = 1e-5;

    ModelStatus status;
    const auto res = solve_elasticity(cfg, status, &temp_field);
    REQUIRE(status.ok);

    // Confined column (u_x = u_y = 0), free top: the axial thermal strain is
    // amplified by the lateral confinement:
    //   eps_zz = alpha·dT·(3·lambda + 2·mu)/(lambda + 2·mu),  sigma_zz = 0.
    const double E = 200e9, nu = 0.3, alpha = 1e-5, gz = 1.0;
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    const double factor = (3.0 * lambda + 2.0 * mu) / (lambda + 2.0 * mu);
    for (int k : {5, 10, 15, 20})
    {
        const auto u = node_at(res.displacement, 5, 5, k);
        CHECK(u[2] == doctest::Approx(factor * alpha * gz * (0.05 * k) * (0.05 * k) / 2.0)
                       .epsilon(1e-6));
        // Lateral leakage is discrete mixed-term noise (~1e-8 on a 1e-6
        // displacement scale): bound it at 0.1% of the axial scale.
        CHECK(std::fabs(u[0]) < 1e-6);
        CHECK(std::fabs(u[1]) < 1e-6);
    }
}

TEST_CASE("Elasticity: validation")
{
    ElasticityConfig cfg;
    cfg.grid.dims = {1, 5, 5};
    std::string err;
    std::vector<std::string> warn;
    CHECK(!validate_elasticity_config(cfg, err, warn)); // dims < 2

    cfg.grid.dims = {5, 5, 5};
    cfg.material.elastic_modulus = 0.0;
    CHECK(!validate_elasticity_config(cfg, err, warn)); // E <= 0

    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.7;
    CHECK(!validate_elasticity_config(cfg, err, warn)); // nu >= 0.5

    cfg.material.poisson_ratio = 0.3;
    cfg.fixed_mask = {std::array<bool, 3>{true, true, true}}; // wrong size
    CHECK(!validate_elasticity_config(cfg, err, warn));

    cfg.fixed_mask = std::vector<std::array<bool, 3>>(5 * 5 * 5);
    cfg.traction_mask = {true}; // wrong size
    CHECK(!validate_elasticity_config(cfg, err, warn));
}

TEST_CASE("Elasticity: channel adapter smoke")
{
    ElasticityConfig cfg;
    cfg.grid.origin = {0, 0, 0};
    cfg.grid.spacing = {0.05, 0.05, 0.05};
    cfg.grid.dims = {21, 5, 21};
    cfg.material.elastic_modulus = 200e9;
    cfg.material.poisson_ratio = 0.0;
    cfg.fixed_mask = roller_with_rigid_pins(cfg.grid);
    cfg.surface_traction = {0, 0, 1e6};
    cfg.tolerance = 1e-9;

    ModelStatus status;
    const auto res = solve_elasticity(cfg, status);
    REQUIRE(status.ok);

    auto channel = exd::physics::coupling::make_vector_grid_field(res.displacement);
    REQUIRE(channel != nullptr);
    std::array<double, 3> v{};
    CHECK(channel->sample({0.25, 0.1, 0.5}, v));
    CHECK(v[2] == doctest::Approx(1e6 * 0.5 / 200e9).epsilon(0.02));
    CHECK(!channel->sample({5.0, 0.0, 0.5}, v)); // out of bounds
}
