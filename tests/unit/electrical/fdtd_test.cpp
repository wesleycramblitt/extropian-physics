// Unit tests for the 3D FDTD (Yee lattice) time-domain Maxwell solver.
// Verifies: allocation/validation, auto-dt stability, pulse speed ~ c,
// sign-reversed echo from the low wall, energy conservation, determinism.
#include <exd/physics/electrical/fdtd.hpp>
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace exd::physics::electrical;

namespace {

constexpr double EPS0 = 8.8541878128e-12;  // F/m
constexpr double MU0 = 1.25663706212e-6;   // H/m

/// 1D-like default config: 200x8x8 cells, auto dt, gaussian soft source
/// (t0 = 30 steps, sigma = 6 steps) at x-plane 8.
FdtdConfig base_config() {
    FdtdConfig c;
    c.dims = {200, 8, 8};
    c.spacing = {1e-3, 1e-3, 1e-3};
    c.source_plane_index = 8;
    c.source_amplitude = 1.0;
    c.source_t0 = 30.0;
    c.source_sigma = 6.0;
    c.max_steps = 400;
    return c;  // dt = 0.0 -> auto (courant_factor x CFL)
}

/// The solver's auto time step for this config (mirrors time_step() in
/// fdtd.cpp, which is not exposed).
double auto_dt(const FdtdConfig& c) {
    const double c0 = 1.0 / std::sqrt(EPS0 * MU0);
    const double dx = c.spacing[0];
    const double dy = c.spacing[1];
    const double dz = c.spacing[2];
    const double cfl = 1.0 / (c0 * std::sqrt(1.0 / (dx * dx)
                                            + 1.0 / (dy * dy)
                                            + 1.0 / (dz * dz)));
    return c.courant_factor * cfl;
}

/// Drives init + step manually, recording per-step Ez at two x-planes
/// (middle of the y/z cross section) alongside the diagnostics history.
struct ProbedRun {
    bool valid = false;
    std::string error;
    FdtdField field;
    std::vector<FdtdStepResult> history;
    std::vector<double> ez_at_a;  // Ez at probe plane a
    std::vector<double> ez_at_b;  // Ez at probe plane b
};

ProbedRun run_probed(const FdtdConfig& config, int probe_a_x, int probe_b_x) {
    ProbedRun run;
    exd::physics::ModelStatus status;
    FdtdField field;
    if (!init_fdtd_field(config, field, status)) {
        run.error = status.error;
        return run;
    }
    const std::size_t nx = static_cast<std::size_t>(config.dims[0]);
    const std::size_t ny = static_cast<std::size_t>(config.dims[1]);
    const std::size_t nz = static_cast<std::size_t>(config.dims[2]);
    const std::size_t jj = ny / 2;
    const std::size_t kk = nz / 2;
    const std::size_t a = static_cast<std::size_t>(probe_a_x) + nx * (jj + ny * kk);
    const std::size_t b = static_cast<std::size_t>(probe_b_x) + nx * (jj + ny * kk);
    for (int step = 0; step < config.max_steps; ++step) {
        FdtdStepResult step_result;
        if (!step_fdtd(config, field, step_result, status)) {
            run.error = status.error;
            return run;
        }
        run.history.push_back(step_result);
        run.ez_at_a.push_back(field.ez[a]);
        run.ez_at_b.push_back(field.ez[b]);
    }
    run.valid = true;
    run.field = std::move(field);
    return run;
}

/// First step where |Ez| exceeds `threshold` (pulse-front arrival).
int first_crossing(const std::vector<double>& v, double threshold) {
    for (std::size_t i = 0; i < v.size(); ++i)
        if (std::abs(v[i]) > threshold)
            return static_cast<int>(i);
    return -1;
}

/// First local maximum of `v` above `threshold` (index out-param).
double first_positive_peak(const std::vector<double>& v, double threshold,
                           int& index) {
    for (std::size_t i = 1; i + 1 < v.size(); ++i) {
        if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > threshold) {
            index = static_cast<int>(i);
            return v[i];
        }
    }
    index = -1;
    return 0.0;
}

} // anonymous namespace

TEST_CASE("FDTD: init_fdtd_field allocates a valid, sized, zeroed field") {
    auto config = base_config();
    FdtdField field;
    exd::physics::ModelStatus status;

    REQUIRE(init_fdtd_field(config, field, status));
    CHECK(status.ok);
    CHECK(field.valid);
    CHECK(field.dims == config.dims);

    const std::size_t n = static_cast<std::size_t>(config.dims[0])
                        * static_cast<std::size_t>(config.dims[1])
                        * static_cast<std::size_t>(config.dims[2]);
    CHECK(field.ex.size() == n);
    CHECK(field.ey.size() == n);
    CHECK(field.ez.size() == n);
    CHECK(field.hx.size() == n);
    CHECK(field.hy.size() == n);
    CHECK(field.hz.size() == n);

    const bool zeroed = std::all_of(field.ex.begin(), field.ex.end(),
                                    [](double x) { return x == 0.0; }) &&
                        std::all_of(field.hz.begin(), field.hz.end(),
                                    [](double x) { return x == 0.0; });
    CHECK(zeroed);

    CHECK(field.step == 0);
    CHECK(field.t == 0.0);
}

TEST_CASE("FDTD: auto dt run is stable and finite (1D-like, 400 steps)") {
    auto config = base_config();
    auto run = run_probed(config, 8, 28);

    REQUIRE(run.valid);
    for (const auto& h : run.history) {
        CHECK(std::isfinite(h.max_e));
        CHECK(std::isfinite(h.max_h));
        CHECK(std::isfinite(h.energy));
        CHECK(h.max_e >= 0.0);
        CHECK(h.max_e <= 10.0);
    }

    // Final field: all six arrays finite (no NaN/inf seeding).
    for (const auto* v : {&run.field.ex, &run.field.ey, &run.field.ez,
                          &run.field.hx, &run.field.hy, &run.field.hz}) {
        for (double x : *v)
            CHECK(std::isfinite(x));
    }

    // Time bookkeeping: t = max_steps * auto dt, step = max_steps.
    CHECK(run.field.t == doctest::Approx(config.max_steps * auto_dt(config))
                             .epsilon(1e-12));
    CHECK(run.field.step == static_cast<std::uint64_t>(config.max_steps));
    CHECK(run.history.size() == static_cast<std::size_t>(config.max_steps));
}

TEST_CASE("FDTD: gaussian pulse travels 20 cells at ~c (within 15%)") {
    auto config = base_config();
    auto run = run_probed(config, 8, 28);

    REQUIRE(run.valid);
    const double c0 = 1.0 / std::sqrt(EPS0 * MU0);
    const double dx = config.spacing[0];
    const double dt = auto_dt(config);

    // Propagation delay for 20 cells (x=8 -> x=28) at the speed of light.
    const double expected_steps = 20.0 * dx / (c0 * dt);

    // Pulse-front arrival at each probe (first crossing of 10% amplitude;
    // argmax picks up later re-circulation/surface-wave lobes in this
    // compressed 200x8x8 box, so we measure the first arrival instead).
    const double threshold = 0.1 * config.source_amplitude;
    const int arrival_a = first_crossing(run.ez_at_a, threshold);
    const int arrival_b = first_crossing(run.ez_at_b, threshold);

    REQUIRE(arrival_a >= 0);
    REQUIRE(arrival_b >= 0);

    const int delay = arrival_b - arrival_a;
    CHECK(delay == doctest::Approx(expected_steps).epsilon(0.15));
}

TEST_CASE("FDTD: pulse reflects off the low wall with a sign-reversed echo") {
    // Probe sits behind the source plane (x=2, source at x=10): first a
    // positive pulse passing on its way to the x- wall, then a later,
    // sign-reversed echo.
    auto config = base_config();
    config.dims = {80, 8, 8};
    config.source_plane_index = 10;
    config.max_steps = 150;
    auto run = run_probed(config, 2, 2);

    REQUIRE(run.valid);

    int pos_step = -1;
    const double pos_peak = first_positive_peak(
        run.ez_at_a, 0.2 * config.source_amplitude, pos_step);

    double neg_peak = 0.0;
    int neg_step = -1;
    for (std::size_t i = 0; i < run.ez_at_a.size(); ++i) {
        if (run.ez_at_a[i] < neg_peak) {
            neg_peak = run.ez_at_a[i];
            neg_step = static_cast<int>(i);
        }
    }

    REQUIRE(pos_step >= 0);
    REQUIRE(neg_step >= 0);
    CHECK(pos_peak > 0.0);
    CHECK(std::abs(neg_peak) > 0.15 * config.source_amplitude);
    CHECK(neg_step > pos_step);
}

TEST_CASE("FDTD: total electromagnetic energy is conserved after the source decays") {
    auto config = base_config();
    auto run = run_probed(config, 8, 28);

    REQUIRE(run.valid);
    REQUIRE(run.history.size() > 350);

    // The gaussian is fully decayed by step ~ t0 + 5*sigma (60); the wave
    // is still inside the box at steps 150 and 350 (lossless PEC walls).
    const double e150 = run.history[150].energy;
    const double e350 = run.history[350].energy;

    CHECK(e150 > 0.0);
    CHECK(e350 / e150 == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("FDTD: invalid config is rejected") {
    FdtdField field;
    exd::physics::ModelStatus status;

    SUBCASE("dims too small") {
        FdtdConfig c = base_config();
        c.dims = {2, 2, 2};
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }
    SUBCASE("spacing zero") {
        FdtdConfig c = base_config();
        c.spacing = {0.0, 0.0, 0.0};
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }
    SUBCASE("eps_r zero") {
        FdtdConfig c = base_config();
        c.eps_r = 0.0;
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }
    SUBCASE("mu_r zero") {
        FdtdConfig c = base_config();
        c.mu_r = 0.0;
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }
    SUBCASE("courant_factor at limit") {
        FdtdConfig c = base_config();
        c.courant_factor = 1.0;
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }
    SUBCASE("source plane out of range") {
        FdtdConfig c = base_config();
        c.source_plane_index = c.dims[0] + 10;
        CHECK_FALSE(init_fdtd_field(c, field, status));
    }

    // Every failure must carry an ok=false status and a message.
    FdtdConfig bad;
    bad.dims = {2, 2, 2};
    CHECK_FALSE(init_fdtd_field(bad, field, status));
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
}

TEST_CASE("FDTD: stepping an invalid or mismatched field is rejected") {
    auto config = base_config();
    FdtdStepResult step_result;
    exd::physics::ModelStatus status;

    SUBCASE("default-constructed (invalid) field") {
        FdtdField field;
        CHECK_FALSE(step_fdtd(config, field, step_result, status));
    }
    SUBCASE("field dims do not match config") {
        FdtdField field;
        REQUIRE(init_fdtd_field(config, field, status));
        field.dims = {20, 8, 8};
        CHECK_FALSE(step_fdtd(config, field, step_result, status));
    }
    SUBCASE("field arrays undersized") {
        FdtdField field;
        REQUIRE(init_fdtd_field(config, field, status));
        field.ez.pop_back();
        CHECK_FALSE(step_fdtd(config, field, step_result, status));
    }

    FdtdField field;
    CHECK_FALSE(step_fdtd(config, field, step_result, status));
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
}

TEST_CASE("FDTD: run_fdtd surfaces init failures and returns its history") {
    SUBCASE("invalid config yields invalid result with error") {
        FdtdConfig c = base_config();
        c.dims = {2, 2, 2};
        auto result = run_fdtd(c);
        CHECK_FALSE(result.valid);
        CHECK_FALSE(result.error.empty());
    }

    SUBCASE("valid config returns field plus per-step history") {
        auto config = base_config();
        config.max_steps = 20;
        auto result = run_fdtd(config);
        REQUIRE(result.valid);
        CHECK(result.field.valid);
        CHECK(result.history.size() == 20);
        CHECK(result.field.step == 20);
    }
}

TEST_CASE("FDTD: two identical runs produce identical energy histories") {
    auto config = base_config();
    auto run_a = run_probed(config, 8, 28);
    auto run_b = run_probed(config, 8, 28);

    REQUIRE(run_a.valid);
    REQUIRE(run_b.valid);
    REQUIRE(run_a.history.size() == run_b.history.size());

    const bool equal = std::equal(
        run_a.history.begin(), run_a.history.end(), run_b.history.begin(),
        [](const FdtdStepResult& x, const FdtdStepResult& y) {
            return std::abs(x.energy - y.energy) <= 1e-12;
        });
    CHECK(equal);
}