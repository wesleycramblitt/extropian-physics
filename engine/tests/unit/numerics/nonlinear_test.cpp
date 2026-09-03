// nonlinear_test.cpp — spec §37: fixed point + Newton-Krylov.
#include <exd/engine/numerics/nonlinear.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine::core;
using namespace exd::engine::numerics;

namespace {
Field scalar_field(std::string name, double v0)
{
    Field f(FieldMetadata{.name = name, .rank = FieldRank::Scalar, .components = 1,
                          .location = FieldLocation::Global}, 1);
    f.set(0, v0);
    return f;
}
} // namespace

TEST_CASE("Fixed point: x = cos(x) converges to the Dottie number")
{
    ModelStatus st;
    Field x = scalar_field("x", 0.5);
    NonlinearSolverConfig cfg;
    cfg.tolerance = 1e-12;
    cfg.max_iterations = 100;
    auto g = [&](const Field& in, Field& out, ModelStatus&) -> bool {
        out.set(0, std::cos(in.at(0)));
        return true;
    };
    auto rep = solve_fixed_point(g, x, cfg, st);
    CHECK(rep.converged);
    // Dottie number ≈ 0.7390851332
    CHECK(x.at(0) == doctest::Approx(0.7390851332151607).epsilon(1e-8));
}

TEST_CASE("Newton-Krylov: x² - 2 = 0 → sqrt(2)")
{
    ModelStatus st;
    Field x = scalar_field("x", 1.0);
    NonlinearSolverConfig cfg;
    cfg.tolerance = 1e-11;
    cfg.max_iterations = 20;
    auto residual = [&](const Field& in, Field& out, ModelStatus&) -> bool {
        out.set(0, in.at(0) * in.at(0) - 2.0);
        return true;
    };
    auto rep = solve_newton(residual, x, cfg, st);
    CHECK(st.ok);
    CHECK(rep.converged);
    CHECK(x.at(0) == doctest::Approx(std::sqrt(2.0)).epsilon(1e-9));
}

TEST_CASE("Newton-Krylov: nonlinear system x²+y²=1, x=0.5")
{
    ModelStatus st;
    Field x(FieldMetadata{.name = "xy", .rank = FieldRank::Vector, .components = 2,
                          .location = FieldLocation::Global}, 2);
    x.set(0, 0.8);
    x.set(1, 0.8);
    NonlinearSolverConfig cfg;
    cfg.tolerance = 1e-10;
    cfg.max_iterations = 30;
    auto residual = [](const Field& in, Field& out, ModelStatus&) -> bool {
        out.set(0, in.at(0) * in.at(0) + in.at(1) * in.at(1) - 1.0);
        out.set(1, in.at(0) - 0.5);
        return true;
    };
    auto rep = solve_newton(residual, x, cfg, st);
    CHECK(st.ok);
    CHECK(rep.converged);
    CHECK(x.at(0) == doctest::Approx(0.5).epsilon(1e-8));
    CHECK(x.at(1) == doctest::Approx(std::sqrt(0.75)).epsilon(1e-8));
}
