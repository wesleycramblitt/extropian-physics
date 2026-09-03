// reactor_test.cpp
// 0D batch reactor: first-order decay, reversible equilibrium, second-order
// kinetics, Arrhenius scaling, conservation and determinism.

#include <exd/engine/physics/reaction/reactor.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <vector>

using namespace exd::engine::physics::reaction;
using exd::engine::core::ModelStatus;

TEST_CASE("Chemistry: first-order decay is exact and mass-conserving")
{
    ChemistryConfig cfg;
    cfg.species = {{"A", 1.0}, {"B", 0.0}};
    cfg.reactions = {{{0}, {1.0}, {1}, {1.0}, /*k=*/1.0, 0.0, 0.0}};
    cfg.dt = 1e-3;
    cfg.end_time = 1.0;

    ModelStatus status;
    const auto res = solve_chemistry(cfg, status);
    REQUIRE(status.ok);

    CHECK(res.final_concentrations[0] == doctest::Approx(std::exp(-1.0)).epsilon(1e-6));
    CHECK(res.final_concentrations[1] == doctest::Approx(1.0 - std::exp(-1.0)).epsilon(1e-6));
    // Mass conservation: c_A + c_B == 1.
    CHECK(res.final_concentrations[0] + res.final_concentrations[1] ==
          doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("Chemistry: reversible reaction reaches equilibrium monotonically")
{
    ChemistryConfig cfg;
    cfg.species = {{"A", 1.0}, {"B", 0.0}};
    cfg.reactions = {// A -> B
                     {{0}, {1.0}, {1}, {1.0}, /*kf=*/2.0, 0.0, 0.0},
                     // B -> A
                     {{1}, {1.0}, {0}, {1.0}, /*kb=*/1.0, 0.0, 0.0}};
    cfg.dt = 1e-3;
    cfg.end_time = 10.0;

    ModelStatus status;
    const auto res = solve_chemistry(cfg, status);
    REQUIRE(status.ok);

    // K = kf/kb = 2 -> c_A_eq = 1/3, c_B_eq = 2/3.
    const double cA_eq = 1.0 / 3.0;
    CHECK(res.final_concentrations[0] == doctest::Approx(cA_eq).epsilon(5e-3));
    CHECK(res.final_concentrations[1] == doctest::Approx(2.0 / 3.0).epsilon(5e-3));

    // Monotone approach: |c_A - c_Aeq| strictly decreases over the last 20%.
    const size_t start = static_cast<size_t>(res.history.size() * 0.8);
    double prev = 1e9;
    for (size_t i = start; i < res.history.size(); ++i)
    {
        const double err = std::fabs(res.history[i][0] - cA_eq);
        CHECK(err <= prev + 1e-12);
        prev = err;
    }
}

TEST_CASE("Chemistry: second-order kinetics match the analytic solution")
{
    ChemistryConfig cfg;
    cfg.species = {{"A", 1.0}, {"B", 0.0}};
    // 2A -> B
    cfg.reactions = {{{0}, {2.0}, {1}, {1.0}, /*k=*/1.0, 0.0, 0.0}};
    cfg.dt = 1e-3;
    cfg.end_time = 1.0;

    ModelStatus status;
    const auto res = solve_chemistry(cfg, status);
    REQUIRE(status.ok);

    // dc_A/dt = -2·k·c_A^2 -> c_A(t) = 1/(1 + 2·k·t).
    CHECK(res.final_concentrations[0] == doctest::Approx(1.0 / 3.0).epsilon(1e-4));
    // Conservation of A units: c_A + 2·c_B == 1.
    CHECK(res.final_concentrations[0] + 2.0 * res.final_concentrations[1] ==
          doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("Chemistry: Arrhenius temperature scaling matches the analytic ratio")
{
    const double Ea = 30000.0; // J/mol
    const double R = 8.314;
    const double A_pre = 1e3;

    const auto run_at = [&](double T, double end, double dt) {
        ChemistryConfig cfg;
        cfg.species = {{"A", 1.0}, {"B", 0.0}};
        cfg.reactions = {{{0}, {1.0}, {1}, {1.0}, 0.0, Ea, A_pre}};
        cfg.temperature = T;
        cfg.dt = dt;
        cfg.end_time = end;
        ModelStatus status;
        const auto res = solve_chemistry(cfg, status);
        REQUIRE(status.ok);
        return res;
    };

    // k = A·exp(-Ea/(R·T)); half-life t_half = ln(2)/k.
    const auto half_life = [&](const auto& res) {
        // Find the first history time where c_A <= 0.5 (linear interpolation).
        const auto& h = res.history;
        const auto& tt = res.time_history;
        for (size_t i = 1; i < h.size(); ++i)
            if (h[i][0] <= 0.5)
            {
                const double f = (0.5 - h[i - 1][0]) / (h[i][0] - h[i - 1][0]);
                return tt[i - 1] + f * (tt[i] - tt[i - 1]);
            }
        return 1e9;
    };

    const auto res300 = run_at(300.0, 300.0, 0.05);
    const auto res600 = run_at(600.0, 1.5, 5e-4);

    const double t300 = half_life(res300);
    const double t600 = half_life(res600);
    REQUIRE(t300 < 200.0);
    REQUIRE(t600 > 1e-4);

    const double ratio_analytic =
        std::exp((-Ea / (R * 600.0)) - (-Ea / (R * 300.0))); // k600/k300
    const double ratio_measured = t300 / t600;               // = k600/k300
    CHECK(ratio_measured == doctest::Approx(ratio_analytic).epsilon(0.10));
}

TEST_CASE("Chemistry: validation")
{
    ChemistryConfig cfg; // species empty
    std::string err;
    std::vector<std::string> warn;
    CHECK(!validate_chemistry_config(cfg, err, warn));

    cfg.species = {{"A", 1.0}};
    cfg.reactions = {{{5}, {1.0}, {0}, {1.0}, 1.0, 0.0, 0.0}}; // index out of range
    CHECK(!validate_chemistry_config(cfg, err, warn));

    cfg.reactions = {{{0}, {-1.0}, {0}, {1.0}, 1.0, 0.0, 0.0}}; // negative stoich
    CHECK(!validate_chemistry_config(cfg, err, warn));

    std::string err2;
    cfg.reactions.clear();
    cfg.dt = 0.0;
    CHECK(!validate_chemistry_config(cfg, err2, warn));
}

TEST_CASE("Chemistry: determinism and mass conservation for a balanced reaction")
{
    ChemistryConfig cfg;
    cfg.species = {{"A", 1.0}, {"B", 1.0}, {"C", 0.0}};
    // A + B -> 2C: atom-balanced AND mole-count-conserving (2 -> 2 molecules).
    cfg.reactions = {{{0, 1}, {1.0, 1.0}, {2}, {2.0}, /*k=*/1.0, 0.0, 0.0}};
    cfg.dt = 1e-3;
    cfg.end_time = 2.0;

    ModelStatus status;
    const auto res = solve_chemistry(cfg, status);
    REQUIRE(status.ok);

    const auto res2 = solve_chemistry(cfg, status);
    REQUIRE(status.ok);
    for (size_t i = 0; i < res.final_concentrations.size(); ++i)
        CHECK(res.final_concentrations[i] == res2.final_concentrations[i]);

    for (const auto& c : res.final_concentrations)
        CHECK(std::isfinite(c));
    // Mole-count conserving reaction: total concentration drift is
    // rounding-level (the metric is |sum_new - sum_old|/(1 + sum_old)).
    CHECK(res.max_total_mass_error < 1e-9);
    // Species never go negative.
    for (const auto& snap : res.history)
        for (double c : snap)
            CHECK(c >= 0.0);
}
