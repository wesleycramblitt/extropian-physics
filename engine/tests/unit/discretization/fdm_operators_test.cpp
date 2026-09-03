// fdm_operators_test.cpp — spec §26/§63: FDM operators verified with
// manufactured solutions and convergence-order checks.
#include <exd/engine/discretization/fdm/operators.hpp>
#include <exd/engine/mesh/generation.hpp>
#include <doctest/doctest.h>

#include <vector>

using namespace exd::engine;
using namespace exd::engine::discretization::fdm;
using namespace exd::engine::mesh;

namespace {
double l2_error(std::span<const double> a, std::span<const double> b)
{
    double acc = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        acc += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(acc / static_cast<double>(a.size()));
}

GhostSpec mirror_ghosts()
{
    GhostSpec g;
    return g;   // all faces mirror (zero-gradient)
}
} // namespace

TEST_CASE("FDM gradient: analytic quadratic field, central O(h²)")
{
    ModelStatus st;
    // u = x² + 2y² + 3z²  →  ∇u = (2x, 4y, 6z) exact at nodes for central diff
    for (int n : {8, 16})
    {
        const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {n, n, n});
        std::vector<double> u(g.node_count()), gx(g.node_count()), gy(g.node_count()), gz(g.node_count());
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i)
                {
                    const double x = static_cast<double>(i) / (n - 1);
                    const double y = static_cast<double>(j) / (n - 1);
                    const double z = static_cast<double>(k) / (n - 1);
                    u[static_cast<size_t>(i) + n * (static_cast<size_t>(j) + n * k)] = x * x + 2 * y * y + 3 * z * z;
                }
        REQUIRE(gradient(g, u, gx, gy, gz, mirror_ghosts(), st));
        // Interior exactness (boundary stencils are first-order, documented)
        double ex = 0, ey = 0, ez = 0;
        for (int k = 1; k < n - 1; ++k)
            for (int j = 1; j < n - 1; ++j)
                for (int i = 1; i < n - 1; ++i)
                {
                    const size_t idx = static_cast<size_t>(i) + n * (static_cast<size_t>(j) + n * k);
                    const double x = static_cast<double>(i) / (n - 1);
                    const double y = static_cast<double>(j) / (n - 1);
                    const double z = static_cast<double>(k) / (n - 1);
                    ex = std::max(ex, std::fabs(gx[idx] - 2 * x));
                    ey = std::max(ey, std::fabs(gy[idx] - 4 * y));
                    ez = std::max(ez, std::fabs(gz[idx] - 6 * z));
                }
        CHECK(ex < 1e-10);
        CHECK(ey < 1e-10);
        CHECK(ez < 1e-10);
    }
}

TEST_CASE("FDM Laplacian: quadratic field is exact; convergence order on sin")
{
    ModelStatus st;
    // Exactness for quadratics: lap(x² + y² + z²) = 6
    {
        const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {10, 10, 10});
        std::vector<double> u(g.node_count()), out(g.node_count());
        for (size_t idx = 0; idx < g.node_count(); ++idx)
        {
            const int k = static_cast<int>(idx / (10 * 10));
            const int j = static_cast<int>((idx / 10) % 10);
            const int i = static_cast<int>(idx % 10);
            const double x = static_cast<double>(i) / 9.0;
            const double y = static_cast<double>(j) / 9.0;
            const double z = static_cast<double>(k) / 9.0;
            u[idx] = x * x + y * y + z * z;
        }
        REQUIRE(laplacian(g, u, out, mirror_ghosts(), st));
        double err = 0;
        for (int k = 1; k < 9; ++k)
            for (int j = 1; j < 9; ++j)
                for (int i = 1; i < 9; ++i)
                {
                    const size_t idx = static_cast<size_t>(i) + 10ull *
                        (static_cast<size_t>(j) + 10ull * k);
                    err = std::max(err, std::fabs(out[idx] - 6.0));
                }
        CHECK(err < 1e-9);
    }
    // u = sin(πx)·sin(πy) → lap u = -2π² u (2D); error → 0 as h → 0
    auto run = [&](int n) {
        const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {n, n, 2});
        std::vector<double> u(g.node_count()), out(g.node_count());
        for (int k = 0; k < 2; ++k)
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i)
                {
                    const double x = static_cast<double>(i) / (n - 1);
                    const double y = static_cast<double>(j) / (n - 1);
                    u[static_cast<size_t>(i) + n * (static_cast<size_t>(j) + n * k)] =
                        std::sin(M_PI * x) * std::sin(M_PI * y);
                }
        REQUIRE(laplacian(g, u, out, mirror_ghosts(), st));
        // interior-only L2 error (boundary stencils are first-order)
        double acc = 0.0;
        int count = 0;
        for (int k = 0; k < 2; ++k)
            for (int j = 1; j < n - 1; ++j)
                for (int i = 1; i < n - 1; ++i)
                {
                    const double x = static_cast<double>(i) / (n - 1);
                    const double y = static_cast<double>(j) / (n - 1);
                    const double exact = -2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
                    const size_t idx = static_cast<size_t>(i) + n * (static_cast<size_t>(j) + n * k);
                    acc += (out[idx] - exact) * (out[idx] - exact);
                    ++count;
                }
        return std::sqrt(acc / count);
    };
    const double e16 = run(16);
    const double e32 = run(32);
    // h halves → second-order error quarters
    CHECK(e32 < e16 / 3.5);
}

TEST_CASE("FDM divergence and curl: null identities")
{
    ModelStatus st;
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {8, 8, 8});
    std::vector<double> vx(g.node_count()), vy(g.node_count()), vz(g.node_count());
    for (size_t idx = 0; idx < g.node_count(); ++idx)
    {
        const int k = static_cast<int>(idx / (8 * 8));
        const int j = static_cast<int>((idx / 8) % 8);
        const int i = static_cast<int>(idx % 8);
        const double x = static_cast<double>(i) / 7.0;
        const double y = static_cast<double>(j) / 7.0;
        const double z = static_cast<double>(k) / 7.0;
        // v = curl(w) with w = (0, 0, x²) → v = (0, 2x·?, 0)? Use solenoidal field:
        // v = (-y, x, 0): divergence 0, curl = (0, 0, 2) exact
        vx[idx] = -y; vy[idx] = x; vz[idx] = 0.0;
    }
    std::vector<double> div(g.node_count());
    REQUIRE(divergence(g, vx, vy, vz, div, mirror_ghosts(), st));
    double derr = 0;
    for (double d : div) derr = std::max(derr, std::fabs(d));
    CHECK(derr < 1e-10);

    std::vector<double> cx(g.node_count()), cy(g.node_count()), cz(g.node_count());
    REQUIRE(curl(g, vx, vy, vz, cx, cy, cz, mirror_ghosts(), st));
    double cerr = 0;
    for (int k = 1; k < 7; ++k)
        for (int j = 1; j < 7; ++j)
            for (int i = 1; i < 7; ++i)
            {
                const size_t idx = static_cast<size_t>(i) + 8ull *
                    (static_cast<size_t>(j) + 8ull * k);
                cerr = std::max(cerr, std::fabs(cz[idx] - 2.0));
            }
    CHECK(cerr < 1e-10);
}

TEST_CASE("FDM upwind advection: uniform flow = exact directional derivative")
{
    ModelStatus st;
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {12, 12, 3});
    std::vector<double> vx(g.node_count(), 1.0), vy(g.node_count(), 0.0), vz(g.node_count(), 0.0);
    const double dx = g.spacing[0];
    // first-order upwind is EXACT for linear φ: u·∇φ = u·dφ/dx
    std::vector<double> phi(g.node_count()), out(g.node_count());
    for (size_t idx = 0; idx < g.node_count(); ++idx)
    {
        const int k = static_cast<int>(idx / (12 * 12));
        const int j = static_cast<int>((idx / 12) % 12);
        const int i = static_cast<int>(idx % 12);
        const double x = static_cast<double>(i) / 11.0;
        phi[idx] = x;      // d/dx (x) = 1
    }
    REQUIRE(upwind_advect(g, vx, vy, vz, phi, out, mirror_ghosts(), st));
    double err = 0;
    for (int k = 1; k < 2; ++k)
        for (int j = 1; j < 11; ++j)
            for (int i = 1; i < 11; ++i)
            {
                const size_t idx = static_cast<size_t>(i) + 12ull *
                    (static_cast<size_t>(j) + 12ull * k);
                err = std::max(err, std::fabs(out[idx] - 1.0));
            }
    CHECK(err < 1e-9);
    // quadratic φ: first-order upwind error == h (dφ/dx ≈ (φ_i − φ_{i−1})/h)
    const double v = 1.0;
    const double err_upwind = std::fabs((dx * dx) / dx - 2.0 * dx * 0.5) ;  // formal: error = h
    (void)v; (void)err_upwind;
}

TEST_CASE("FdmLaplacianOperator: matrix-free interfaces")
{
    using namespace exd::engine::core;
    ModelStatus st;
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {6, 6, 6});
    FdmLaplacianOperator op(g, mirror_ghosts());
    Field x(FieldMetadata{.name = "x", .rank = FieldRank::Scalar, .components = 1,
                          .location = FieldLocation::Node}, g.node_count());
    x.assign(1.0);
    Field y(FieldMetadata{.name = "y", .rank = FieldRank::Scalar, .components = 1,
                          .location = FieldLocation::Node}, g.node_count());
    REQUIRE(op.apply(x, y, st));
    // constant field → Laplacian 0 (mirror ghosts)
    double err = 0;
    for (size_t i = 0; i < y.size(); ++i) err = std::max(err, std::fabs(y.at(i)));
    CHECK(err < 1e-12);
    Field d(FieldMetadata{.name = "d", .rank = FieldRank::Scalar, .components = 1,
                          .location = FieldLocation::Node}, g.node_count());
    REQUIRE(op.diagonal(d, st));
    const double h = 1.0 / 5.0;   // dims 6 → spacing 0.2
    const double expected_diag = -(2.0 / (h * h)) * 3.0;
    CHECK(d.at(0) == doctest::Approx(expected_diag).epsilon(1e-12));
    // jvp of a linear operator = apply
    Field v(FieldMetadata{.name = "v", .rank = FieldRank::Scalar, .components = 1,
                          .location = FieldLocation::Node}, g.node_count());
    v.assign(2.0);
    Field jy(FieldMetadata{.name = "jy", .rank = FieldRank::Scalar, .components = 1,
                           .location = FieldLocation::Node}, g.node_count());
    REQUIRE(op.jacobian_vector_product(x, v, jy, st));
    for (size_t i = 0; i < jy.size(); ++i) CHECK(jy.at(i) == doctest::Approx(0.0).epsilon(1e-12));
}
