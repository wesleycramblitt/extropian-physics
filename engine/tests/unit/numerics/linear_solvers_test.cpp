// linear_solvers_test.cpp — spec §35/§36: CG/GMRES/BiCGSTAB on the
// matrix-free operator interface; Poisson verified against analytic.
#include <exd/engine/discretization/fdm/operators.hpp>
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/numerics/bicgstab.hpp>
#include <exd/engine/numerics/cg.hpp>
#include <exd/engine/numerics/gmres.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::core;
using namespace exd::engine::discretization::fdm;
using namespace exd::engine::mesh;
using namespace exd::engine::numerics;

namespace {

/// Dirichlet ghosts on all faces with the analytic boundary values of the
/// manufactured solution u = sin(πx)sin(πy)sin(πz) (zero on the box walls).
GhostSpec dirichlet_zero()
{
    GhostSpec g;
    for (auto& f : g.faces)
    {
        f.dirichlet = true;
        f.value = 0.0;
    }
    return g;
}

Field make_field(const std::string& name, size_t n)
{
    return Field(FieldMetadata{.name = name, .rank = FieldRank::Scalar,
                               .components = 1, .location = FieldLocation::Node}, n);
}

} // namespace

TEST_CASE("CG: Poisson with Dirichlet BC matches analytic solution")
{
    ModelStatus st;
    const int n = 16;
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {n, n, n});
    FdmLaplacianOperator lap(g, dirichlet_zero());

    const size_t N = g.node_count();
    Field rhs = make_field("rhs", N);
    Field sol = make_field("sol", N);
    sol.assign(0.0);

    // −Δu = f with u = sin(πx)sin(πy)sin(πz): f = 3π² u.
    // The FDM Laplacian is Δ (negative definite); CG needs the SPD Poisson
    // operator −Δ, composed with NegatedOperator (§35 operator composition).
    NegatedOperator poisson(lap);
    for (size_t idx = 0; idx < N; ++idx)
    {
        const int k = static_cast<int>(idx / (n * n));
        const int j = static_cast<int>((idx / n) % n);
        const int i = static_cast<int>(idx % n);
        const double x = static_cast<double>(i) / (n - 1);
        const double y = static_cast<double>(j) / (n - 1);
        const double z = static_cast<double>(k) / (n - 1);
        rhs.set(idx, 3.0 * M_PI * M_PI * std::sin(M_PI * x) *
                     std::sin(M_PI * y) * std::sin(M_PI * z));
    }

    IterativeSolverConfig cfg;
    cfg.tolerance = 1e-8;
    cfg.max_iterations = 500;
    auto rep = solve_cg(poisson, rhs, sol, cfg, st);
    REQUIRE(st.ok);
    CHECK(rep.converged);

    // verify against analytic solution
    double l2 = 0.0;
    for (size_t idx = 0; idx < N; ++idx)
    {
        const int k = static_cast<int>(idx / (n * n));
        const int j = static_cast<int>((idx / n) % n);
        const int i = static_cast<int>(idx % n);
        const double x = static_cast<double>(i) / (n - 1);
        const double y = static_cast<double>(j) / (n - 1);
        const double z = static_cast<double>(k) / (n - 1);
        const double exact = std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
        l2 += (sol.at(idx) - exact) * (sol.at(idx) - exact);
    }
    l2 = std::sqrt(l2 / static_cast<double>(N));
    // 2nd-order discretization on n=16: h=1/15 → error ~ (π²h²/6)·scale
    CHECK(l2 < 2e-2);
}

TEST_CASE("GMRES and BiCGSTAB: SPD Poisson agree with CG")
{
    ModelStatus st;
    const int n = 10;
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {n, n, n});
    FdmLaplacianOperator lap(g, dirichlet_zero());
    NegatedOperator poisson(lap);      // SPD: −Δ
    const size_t N = g.node_count();

    Field rhs = make_field("rhs", N);
    rhs.assign(0.0);
    for (int i = 0; i < n; ++i)
        rhs.set(static_cast<size_t>(i), 1.0);   // line source

    Field x_cg = make_field("x_cg", N);
    Field x_gmres = make_field("x_gmres", N);
    Field x_bicg = make_field("x_bicg", N);

    IterativeSolverConfig cfg;
    cfg.tolerance = 1e-9;
    cfg.max_iterations = 300;
    auto rep_cg = solve_cg(poisson, rhs, x_cg, cfg, st);
    REQUIRE(st.ok);
    REQUIRE(rep_cg.converged);

    GmresConfig gcfg;
    gcfg.tolerance = 1e-9;
    gcfg.max_iterations = 300;
    gcfg.restart = 20;
    auto rep_gmres = solve_gmres(poisson, rhs, x_gmres, gcfg, st);
    REQUIRE(st.ok);
    REQUIRE(rep_gmres.converged);

    auto rep_bicg = solve_bicgstab(poisson, rhs, x_bicg, cfg, st);
    REQUIRE(st.ok);
    REQUIRE(rep_bicg.converged);

    // all three agree
    for (size_t i = 0; i < N; ++i)
    {
        CHECK(x_gmres.at(i) == doctest::Approx(x_cg.at(i)).epsilon(1e-6));
        CHECK(x_bicg.at(i) == doctest::Approx(x_cg.at(i)).epsilon(1e-6));
    }
}

TEST_CASE("Field algebra: dot/axpy/norm")
{
    const size_t N = 4;
    Field a = make_field("a", N);
    Field b = make_field("b", N);
    for (size_t i = 0; i < N; ++i)
    {
        a.set(i, static_cast<double>(i) + 1.0);
        b.set(i, 2.0);
    }
    CHECK(dot(a.data(), b.data()) == doctest::Approx(20.0).epsilon(1e-12));
    CHECK(norm2(a.data()) == doctest::Approx(std::sqrt(30.0)).epsilon(1e-12));
    axpy(1.0, a.data(), b.data());
    CHECK(b.at(0) == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(b.at(3) == doctest::Approx(6.0).epsilon(1e-12));
}
