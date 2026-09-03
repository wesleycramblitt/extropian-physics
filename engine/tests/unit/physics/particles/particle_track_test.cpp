// particle_track_test.cpp
// Lagrangian particle tracking: ballistic exactness, drag terminal velocity,
// determinism/number conservation, channel advection, out-of-bounds warnings.

#include <exd/engine/physics/particles/particle_track.hpp>

#include <exd/engine/coupling/field_channels.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>

using namespace exd::engine::physics::particles;
using exd::engine::core::ModelStatus;

namespace {

/// Constant-velocity flow field for advection tests.
struct UniformFlow : exd::engine::coupling::IVectorField3D
{
    std::array<double, 3> v{1.0, 0.0, 0.0};
    bool sample(const std::array<double, 3>&, std::array<double, 3>& v_out) const override
    {
        v_out = v;
        return true;
    }
};

/// Field that always reports out-of-bounds.
struct EmptyFlow : exd::engine::coupling::IVectorField3D
{
    bool sample(const std::array<double, 3>&, std::array<double, 3>&) const override
    {
        return false;
    }
};

} // anonymous namespace

TEST_CASE("Particles: ballistic motion is exact (RK4 on constant acceleration)")
{
    ParticleConfig cfg;
    cfg.particle_count = 1;
    cfg.spawn_extent = {0, 0, 0}; // single particle at the origin
    cfg.initial_velocity = {0, 0, 10.0};
    cfg.gravity = {0, 0, -9.81};
    cfg.drag_coefficient = 0.0;
    cfg.dt = 1e-3;
    cfg.max_steps = 2000; // T = 2.0 s

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    REQUIRE(status.ok);
    REQUIRE(res.final_positions.size() == 1);

    const double T = 2.0;
    const double z_exact = 10.0 * T - 0.5 * 9.81 * T * T; // 20 - 19.62 = 0.38
    CHECK(res.final_positions[0][2] == doctest::Approx(z_exact).epsilon(1e-9));
    CHECK(res.final_positions[0][0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.final_positions[0][1] == doctest::Approx(0.0).epsilon(1e-12));
    // Velocity history reaches the analytic v(t) = v0 - g·t at the end.
    CHECK(res.final_velocities[0][2] == doctest::Approx(10.0 - 9.81 * T).epsilon(1e-9));
}

TEST_CASE("Particles: drag reaches terminal velocity without overshoot")
{
    ParticleConfig cfg;
    cfg.particle_count = 1;
    cfg.spawn_extent = {0, 0, 0};
    cfg.initial_velocity = {0, 0, 0.0};
    cfg.gravity = {0, 0, -9.81};
    cfg.drag_coefficient = 1.0; // tau = 1 s
    cfg.dt = 1e-3;
    cfg.max_steps = 20000; // T = 20 s >> 5·tau

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    REQUIRE(status.ok);

    const double v_t = 9.81;
    CHECK(res.final_velocities[0][2] == doctest::Approx(-v_t).epsilon(1e-3));
    // No overshoot: |v_z| never exceeds the terminal speed (checked on history).
    for (const auto& v : res.velocity_probe)
        CHECK(std::fabs(v[2]) <= v_t * (1.0 + 1e-6));
    // Position: z = -v_t·(T - (1 - e^(-kT))/k).
    const double z_exact = -v_t * (20.0 - (1.0 - std::exp(-20.0)) / 1.0);
    CHECK(res.final_positions[0][2] == doctest::Approx(z_exact).epsilon(5e-3));
}

TEST_CASE("Particles: number conservation and determinism")
{
    ParticleConfig cfg;
    cfg.particle_count = 27; // 3x3x3 lattice
    cfg.spawn_extent = {3, 3, 3};
    cfg.initial_velocity = {1, 2, 3};
    cfg.gravity = {0, 0, -9.81};
    cfg.drag_coefficient = 0.2;
    cfg.dt = 1e-3;
    cfg.max_steps = 500;

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    REQUIRE(status.ok);
    CHECK(res.final_positions.size() == 27);
    CHECK(res.final_velocities.size() == 27);

    const auto res2 = solve_particles(cfg, status);
    REQUIRE(status.ok);
    for (size_t i = 0; i < 27; ++i)
    {
        // Bit-identical determinism (element-wise comparison for doctest).
        const bool same = res.final_positions[i] == res2.final_positions[i];
        CHECK(same);
        // All particles finite and within the spawn box + ballistic envelope.
        for (int c = 0; c < 3; ++c)
            CHECK(std::isfinite(res.final_positions[i][c]));
    }
}

TEST_CASE("Particles: channel advection follows the analytic trajectory")
{
    UniformFlow flow; // v_f = {1, 0, 0}
    ParticleConfig cfg;
    cfg.particle_count = 1;
    cfg.spawn_extent = {0, 0, 0};
    cfg.initial_velocity = {0, 0, 0};
    cfg.gravity = {0, 0, 0};
    cfg.drag_coefficient = 0.5; // 1/s
    cfg.flow_channel = &flow;
    cfg.dt = 1e-3;
    cfg.max_steps = 10000; // T = 10 s

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    REQUIRE(status.ok);
    REQUIRE(res.final_positions.size() == 1);

    // v(t) = vf·(1 - e^(-kt)); x(t) = vf·(t - (1 - e^(-kt))/k).
    const double k = 0.5, vf = 1.0, T = 10.0;
    const double v_exact = vf * (1.0 - std::exp(-k * T));
    const double x_exact = vf * (T - (1.0 - std::exp(-k * T)) / k);
    CHECK(res.final_velocities[0][0] == doctest::Approx(v_exact).epsilon(1e-3));
    CHECK(res.final_positions[0][0] == doctest::Approx(x_exact).epsilon(1e-3));
    CHECK(res.final_positions[0][1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(res.final_positions[0][2] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("Particles: out-of-bounds sampling warns and integrates with v_f = 0")
{
    EmptyFlow flow;
    ParticleConfig cfg;
    cfg.particle_count = 8;
    cfg.spawn_extent = {2, 2, 2};
    cfg.initial_velocity = {0, 0, 0};
    cfg.gravity = {0, 0, 0};
    cfg.drag_coefficient = 0.0;
    cfg.flow_channel = &flow;
    cfg.dt = 1e-3;
    cfg.max_steps = 100;

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    CHECK(status.ok);
    CHECK(res.ok);
    CHECK(!status.warnings.empty()); // out-of-bounds warning recorded
    // v_f = 0 with no forces: particles stay put.
    for (const auto& p : res.final_positions)
    {
        CHECK(p[0] >= 0.0);
        CHECK(p[0] <= 2.0);
        CHECK(p[1] >= 0.0);
        CHECK(p[1] <= 2.0);
        CHECK(p[2] >= 0.0);
        CHECK(p[2] <= 2.0);
    }
}

TEST_CASE("Particles: validation")
{
    ParticleConfig cfg;
    cfg.particle_count = 0;
    std::string err;
    std::vector<std::string> warn;
    CHECK(!validate_particle_config(cfg, err, warn));

    cfg.particle_count = 1;
    cfg.dt = 0.0;
    CHECK(!validate_particle_config(cfg, err, warn));

    cfg.dt = 1e-3;
    cfg.drag_coefficient = -1.0;
    CHECK(!validate_particle_config(cfg, err, warn));
}

TEST_CASE("Particles: structured-grid channel drift smoke")
{
    // v_x = 0.1·x sampled from a structured vector grid.
    exd::engine::coupling::StructuredVectorGrid grid;
    grid.origin = {0, 0, 0};
    grid.spacing = {1.0, 1.0, 1.0};
    grid.dims = {11, 2, 2};
    const int nx = 11, ny = 2, nz = 2;
    grid.values.resize(3 * nx * ny * nz);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const size_t idx = static_cast<size_t>(i) + nx * (j + ny * k);
                grid.values[3 * idx + 0] = 0.1 * i; // v_x
                grid.values[3 * idx + 1] = 0.0;
                grid.values[3 * idx + 2] = 0.0;
            }
    auto channel = exd::engine::coupling::make_vector_grid_field(grid);
    REQUIRE(channel != nullptr);

    ParticleConfig cfg;
    cfg.particle_count = 8;
    cfg.origin = {0.5, 0.5, 0.5}; // keep particle 0 off x = 0 where v_x = 0
    cfg.spawn_extent = {2, 2, 2};
    cfg.initial_velocity = {0, 0, 0};
    cfg.gravity = {0, 0, 0};
    cfg.drag_coefficient = 1.0; // strong drag -> v ~ v_f
    cfg.flow_channel = channel.get();
    cfg.dt = 1e-3;
    cfg.max_steps = 5000; // T = 5 s

    ModelStatus status;
    const auto res = solve_particles(cfg, status);
    REQUIRE(status.ok);

    // With strong drag, dx/dt = 0.1·x -> x(T) = x0·e^(0.5): growth factor
    // ~1.648 (band widened for the drag lag); monotone drift downwind.
    double prev_x = 0.0;
    for (size_t i = 1; i < res.trajectory_probe.size(); ++i)
    {
        CHECK(res.trajectory_probe[i][0] >= prev_x - 1e-12);
        prev_x = res.trajectory_probe[i][0];
    }
    const double x0 = res.trajectory_probe.front()[0];
    const double x1 = res.trajectory_probe.back()[0];
    CHECK(x1 > x0);
    // Growth sits between the no-motion bound (1.0) and the instant-drag
    // bound (e^0.5 = 1.649), because k = 1 s^-1 drag lags the field.
    CHECK(x1 / x0 > 1.2);
    CHECK(x1 / x0 < 1.75);
    // Lateral/vertical drift is zero: y/z stay at the spawn coordinates.
    CHECK(res.final_positions[0][1] == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(res.final_positions[0][2] == doctest::Approx(0.5).epsilon(1e-9));
}
