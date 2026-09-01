// wave_solver_test.cpp
// Phase I acoustics domain: 1D plane-wave speed, 3D box-mode frequencies,
// leapfrog energy conservation, the CFL stability guard, validation and the
// grid-channel adapter smoke test.

#include <exd/physics/acoustics/wave_solver.hpp>
#include <exd/physics/coupling/field_channels.hpp>   // test-only: wrap the result grid

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace exd::physics;
using namespace exd::physics::acoustics;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Flat index: i + nx*(j + ny*k).
std::size_t flat_index(const WaveGridConfig& g, int i, int j, int k)
{
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    return static_cast<std::size_t>(i + nx * (j + ny * k));
}

// The box-mode seed the solver uses (dims==2 axes contribute constant 1).
double mode_shape(const WaveGridConfig& g, const std::array<int32_t, 3>& mode,
                  int i, int j, int k)
{
    const int nx = g.dims[0];
    const int ny = g.dims[1];
    const int nz = g.dims[2];
    const double Lx = g.spacing[0] * (nx - 1);
    const double Ly = g.spacing[1] * (ny - 1);
    const double Lz = g.spacing[2] * (nz - 1);
    const double fx = (nx > 2) ? std::sin(mode[0] * kPi * (g.spacing[0] * i) / Lx) : 1.0;
    const double fy = (ny > 2) ? std::sin(mode[1] * kPi * (g.spacing[1] * j) / Ly) : 1.0;
    const double fz = (nz > 2) ? std::sin(mode[2] * kPi * (g.spacing[2] * k) / Lz) : 1.0;
    return fx * fy * fz;
}

// Sign-change times of a probe signal, linearly interpolated between samples.
std::vector<double> crossing_times(const std::vector<double>& hist, double dt)
{
    std::vector<double> times;
    int last_sign = 0;
    for (std::size_t i = 0; i + 1 < hist.size(); ++i)
    {
        const double v0 = hist[i];
        const double v1 = hist[i + 1];
        if (v0 != 0.0)
            last_sign = (v0 > 0.0) ? 1 : -1;
        if (v1 != 0.0)
        {
            const int s1 = (v1 > 0.0) ? 1 : -1;
            if (last_sign != 0 && last_sign != s1)
            {
                const double t = (static_cast<double>(i) + (0.0 - v0) / (v1 - v0)) * dt;
                times.push_back(t);
                last_sign = s1;
            }
        }
    }
    return times;
}

bool all_finite(const std::vector<double>& v)
{
    for (double x : v)
        if (!std::isfinite(x))
            return false;
    return true;
}

// Absolute Pearson correlation of the final field with the mode shape.
double mode_correlation(const WaveResult& result, const WaveConfig& config)
{
    WaveGridConfig g;   // same layout as the result's StructuredScalarGrid
    g.origin = result.pressure_final.origin;
    g.spacing = result.pressure_final.spacing;
    g.dims = result.pressure_final.dims;
    double num = 0.0, den_p = 0.0, den_m = 0.0;
    for (int k = 0; k < g.dims[2]; ++k)
        for (int j = 0; j < g.dims[1]; ++j)
            for (int i = 0; i < g.dims[0]; ++i)
            {
                const std::size_t I = flat_index(g, i, j, k);
                const double p = result.pressure_final.values[I];
                const double m = mode_shape(g, config.initial_mode, i, j, k);
                num += p * m;
                den_p += p * p;
                den_m += m * m;
            }
    return std::abs(num) / std::sqrt(den_p * den_m);
}

} // anonymous namespace

TEST_CASE("acoustics: 1D plane wave speed")
{
    // Slab 10 m long along x (dx = 0.05, nx = 201), thin y/z axes (dims = 2).
    // Pressure-release walls along x, symmetry planes along y/z: the l = 1
    // mode is a plane wave with period T = 2L/(c*l).
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {201, 2, 2};
    config.sound_speed = 343.0;
    config.initial_mode = {1, 1, 1};
    config.amplitude = 1.0;

    // Probe at node 50 of 200 (25% from the left wall, x = 2.5 m).
    config.probe_index = static_cast<int32_t>(flat_index(config.grid, 50, 0, 0));
    config.max_steps = 1000;   // ~= 1.15 periods

    ModelStatus status;
    const WaveResult result = solve_wave(config, status);
    REQUIRE(result.ok);
    REQUIRE(status.ok);

    const double T_an = 2.0 * 10.0 / 343.0;          // ~= 0.058309 s
    const double T_quarter = T_an / 4.0;

    const std::vector<double> crossings = crossing_times(result.probe_history, result.dt_used);
    REQUIRE(crossings.size() >= 2);

    // First zero crossing ~= T/4.
    CHECK(std::abs(crossings[0] - T_quarter) < 0.01 * T_an);
    // Periodic time from the first two crossings ~= T.
    const double T_meas = 2.0 * (crossings[1] - crossings[0]);
    CHECK(std::abs(T_meas - T_an) < 0.01 * T_an);

    // The first maximum after the first crossing returns at t ~= T.
    std::size_t max_index = 0;
    for (std::size_t i = 1; i < result.probe_history.size(); ++i)
        if (result.probe_history[i] > result.probe_history[max_index])
            max_index = i;
    const double t_max = static_cast<double>(max_index) * result.dt_used;
    CHECK(std::abs(t_max - T_an) < 0.01 * T_an);

    // Amplitude decay < 2% over one period (probe returns to its t0 value).
    const double initial_probe = result.probe_history.front();
    const double returned_probe = result.probe_history[max_index];
    CHECK(std::abs(returned_probe / initial_probe - 1.0) < 0.02);

    // And the same check via the field energy over exactly one numerical
    // period (the crossing-interval period measured above).
    config.max_steps = static_cast<uint64_t>(std::llround(T_meas / result.dt_used));
    ModelStatus status2;
    const WaveResult one_period = solve_wave(config, status2);
    REQUIRE(one_period.ok);
    CHECK(one_period.energy_final / one_period.energy_initial > 0.98);
    CHECK(one_period.energy_final / one_period.energy_initial < 1.02);
}

TEST_CASE("acoustics: 3D box eigenfrequencies")
{
    // Unit box 21^3 (dx = 0.05), c = 1000, all walls pressure-release.
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.amplitude = 1.0;
    config.max_steps = 2000;

    // ---- Mode (1,1,1): f = c/2*sqrt(3) ~= 866.0 Hz.  Probe at the center.
    config.initial_mode = {1, 1, 1};
    config.probe_index = static_cast<int32_t>(flat_index(config.grid, 10, 10, 10));
    ModelStatus status;
    const WaveResult res111 = solve_wave(config, status);
    REQUIRE(res111.ok);

    const double t_total = static_cast<double>(res111.probe_history.size() - 1) * res111.dt_used;
    const auto crossings111 = crossing_times(res111.probe_history, res111.dt_used);
    REQUIRE(crossings111.size() >= 4);
    const double f111 = static_cast<double>(crossings111.size()) / (2.0 * t_total);
    const double f_an_111 = 1000.0 / 2.0 * std::sqrt(3.0);
    CHECK(std::abs(f111 - f_an_111) / f_an_111 < 0.02);

    // The final field is the standing (1,1,1) mode.
    CHECK(mode_correlation(res111, config) > 0.9);

    // ---- Mode (1,1,2): f = c/2*sqrt(6) ~= 1224.7 Hz in a unit box (the
    // sqrt(5) in the task brief is incorrect for (1,1,2); sqrt(5)
    // corresponds to a 2D mode).  Probe at the z = L/4 antinode (the
    // z = L/2 plane is a nodal plane for n = 2).
    config.initial_mode = {1, 1, 2};
    config.probe_index = static_cast<int32_t>(flat_index(config.grid, 10, 10, 5));
    ModelStatus status2;
    const WaveResult res112 = solve_wave(config, status2);
    REQUIRE(res112.ok);

    const double t_total2 = static_cast<double>(res112.probe_history.size() - 1) * res112.dt_used;
    const auto crossings112 = crossing_times(res112.probe_history, res112.dt_used);
    REQUIRE(crossings112.size() >= 4);
    const double f112 = static_cast<double>(crossings112.size()) / (2.0 * t_total2);
    const double f_an_112 = 1000.0 / 2.0 * std::sqrt(6.0);
    CHECK(std::abs(f112 - f_an_112) / f_an_112 < 0.02);

    CHECK(mode_correlation(res112, config) > 0.9);
}

TEST_CASE("acoustics: energy conservation (non-decaying)")
{
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.initial_mode = {1, 1, 1};
    config.amplitude = 1.0;
    config.probe_index = static_cast<int32_t>(flat_index(config.grid, 10, 10, 10));

    // First, measure the numerical period from a 24-period run.
    config.max_steps = 1200;
    ModelStatus status;
    const WaveResult measure = solve_wave(config, status);
    REQUIRE(measure.ok);
    const double t_total = static_cast<double>(measure.probe_history.size() - 1) * measure.dt_used;
    const auto crossings = crossing_times(measure.probe_history, measure.dt_used);
    REQUIRE(crossings.size() >= 4);
    const double f_num = static_cast<double>(crossings.size()) / (2.0 * t_total);
    const double T_num = 1.0 / f_num;

    // Then run exactly 10 numerical periods: the leapfrog scheme has no
    // amplitude loss for the eigenmodes, so sum(p^2)*cell_volume returns to
    // its initial value to within a few percent.
    config.max_steps = static_cast<uint64_t>(std::llround(10.0 * T_num / measure.dt_used));
    CHECK(config.max_steps > 0);

    ModelStatus status2;
    const WaveResult result = solve_wave(config, status2);
    REQUIRE(result.ok);
    const double ratio = result.energy_final / result.energy_initial;
    CHECK(ratio > 0.97);
    CHECK(ratio < 1.03);
}

TEST_CASE("acoustics: CFL guard clamps and warns")
{
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.initial_mode = {1, 1, 1};
    config.probe_index = static_cast<int32_t>(flat_index(config.grid, 10, 10, 10));
    config.max_steps = 200;

    // The exact von Neumann bound for this grid.
    const double cfl_limit = 1.0 /
        (config.sound_speed *
         std::sqrt(1.0 / (0.05 * 0.05) + 1.0 / (0.05 * 0.05) + 1.0 / (0.05 * 0.05)));

    config.dt = 2.0 * cfl_limit;   // over the limit: the guard must engage
    ModelStatus status;
    const WaveResult result = solve_wave(config, status);
    REQUIRE(result.ok);
    REQUIRE(status.ok);

    bool warned = false;
    for (const std::string& w : status.warnings)
        if (w.find("CFL violation") != std::string::npos)
            warned = true;
    CHECK(warned);
    CHECK(result.dt_used == doctest::Approx(cfl_limit).epsilon(1e-12));

    // Clamped scheme stays finite and bounded (a genuinely unstable run
    // would blow up).
    CHECK(all_finite(result.pressure_final.values));
    CHECK(all_finite(result.probe_history));
    CHECK(result.max_pressure < 100.0 * config.amplitude);
}

TEST_CASE("acoustics: validation reports errors, never throws")
{
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;

    std::string error;
    std::vector<std::string> warnings;

    // dims < 2 per axis
    config.grid.dims = {1, 3, 3};
    CHECK_FALSE(validate_wave_config(config, error, warnings));
    CHECK(error.find("dims") != std::string::npos);
    ModelStatus status;
    const WaveResult bad1 = solve_wave(config, status);
    CHECK_FALSE(bad1.ok);
    CHECK_FALSE(status.ok);
    CHECK(bad1.pressure_final.values.empty());

    // sound_speed <= 0
    config = WaveConfig{};
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 0.0;
    CHECK_FALSE(validate_wave_config(config, error, warnings));
    const WaveResult bad2 = solve_wave(config, status);
    CHECK_FALSE(bad2.ok);
    CHECK(error.find("sound speed") != std::string::npos);

    // probe_index out of range
    config = WaveConfig{};
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.probe_index = 21 * 21 * 21;   // one past the end
    CHECK_FALSE(validate_wave_config(config, error, warnings));
    CHECK(error.find("probe_index") != std::string::npos);

    // dt < 0
    config = WaveConfig{};
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.dt = -1.0;
    CHECK_FALSE(validate_wave_config(config, error, warnings));
}

TEST_CASE("acoustics: channel adapter smoke test")
{
    // Short run: the standing wave is ~ cos(omega*2*dt) ~= 0.97 near its
    // peak, so the antinode sample is ~= amplitude and the nodal-plane
    // sample ~= 0.
    WaveConfig config;
    config.grid.spacing = {0.05, 0.05, 0.05};
    config.grid.dims = {21, 21, 21};
    config.sound_speed = 1000.0;
    config.initial_mode = {1, 1, 1};
    config.amplitude = 1.5;
    config.max_steps = 2;

    ModelStatus status;
    const WaveResult result = solve_wave(config, status);
    REQUIRE(result.ok);

    auto channel = coupling::make_scalar_grid_field(result.pressure_final);
    REQUIRE(channel != nullptr);

    double value = 0.0;
    // Antinode at the box center.
    CHECK(channel->sample({0.5, 0.5, 0.5}, value));
    CHECK(value > 0.9 * config.amplitude);
    CHECK(value <= config.amplitude + 1e-12);

    // A node (the z = 0 pressure-release wall) samples ~= 0.
    CHECK(channel->sample({0.5, 0.5, 0.0}, value));
    CHECK(std::abs(value) < 0.05 * config.amplitude);

    // Out of bounds is handled gracefully by the channel.
    CHECK_FALSE(channel->sample({2.0, 0.5, 0.5}, value));
}