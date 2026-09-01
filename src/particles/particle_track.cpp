// particle_track.cpp
// Deterministic Lagrangian particle tracking over sampled fields.  Spawn is a
// fixed lattice (largest product decomposition <= requested count), dynamics
// is integrated with the shared solver::integrate_step() integrator over one
// flat [particle states...] vector.

#include <exd/physics/particles/particle_track.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace exd::physics::particles {

namespace {

/// Largest a*b*c <= count with a <= b <= c, balanced (prefer the largest a,
/// then largest b among ties, so 27 -> (3,3,3) and 100 -> (4,5,5)).
std::array<int32_t, 3> lattice_dims(uint64_t count)
{
    int32_t ba = 1, bb = 1, bc = 1;
    uint64_t best = 1;
    for (uint64_t a = 1; a * a * a <= count; ++a)
    {
        for (uint64_t b = a; a * b * b <= count; ++b)
        {
            const uint64_t c = count / (a * b);
            if (c < b) continue;
            const uint64_t p = a * b * c;
            if (p > best || (p == best && a > static_cast<uint64_t>(ba))
                || (p == best && a == static_cast<uint64_t>(ba)
                    && b > static_cast<uint64_t>(bb)))
            {
                best = p;
                ba = static_cast<int32_t>(a);
                bb = static_cast<int32_t>(b);
                bc = static_cast<int32_t>(c);
            }
        }
    }
    return {ba, bb, bc};
}

} // namespace

bool validate_particle_config(const ParticleConfig& config,
                              std::string& error,
                              std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.particle_count < 1)
    {
        error = "particles: particle_count must be >= 1";
        return false;
    }
    if (!(config.dt > 0.0))
    {
        error = "particles: dt must be > 0";
        return false;
    }
    if (config.drag_coefficient < 0.0)
    {
        error = "particles: drag_coefficient must be >= 0";
        return false;
    }
    if (config.max_steps == 0)
    {
        error = "particles: max_steps must be > 0";
        return false;
    }
    if (config.history_interval == 0)
    {
        error = "particles: history_interval must be >= 1";
        return false;
    }
    for (int c = 0; c < 3; ++c)
    {
        if (!std::isfinite(config.gravity[static_cast<std::size_t>(c)])
            || !std::isfinite(config.initial_velocity[static_cast<std::size_t>(c)])
            || !std::isfinite(config.origin[static_cast<std::size_t>(c)])
            || !std::isfinite(config.spawn_extent[static_cast<std::size_t>(c)]))
        {
            error = "particles: origin, spawn_extent, initial_velocity and gravity must be finite";
            return false;
        }
        if (config.spawn_extent[static_cast<std::size_t>(c)] < 0.0)
        {
            error = "particles: spawn_extent components must be >= 0";
            return false;
        }
    }

    if (config.flow_channel == nullptr && config.drag_coefficient > 0.0)
        warnings.push_back("particles: drag without a flow channel; v_f = 0 everywhere");
    return true;
}

ParticleResult solve_particles(const ParticleConfig& config, ModelStatus& status)
{
    ParticleResult res;
    std::string error;
    std::vector<std::string> warnings;
    if (!validate_particle_config(config, error, warnings))
    {
        res.status = ModelStatus{false, error, warnings};
        status = res.status;
        return res;
    }

    const std::array<int32_t, 3> ld = lattice_dims(config.particle_count);
    const int64_t np_i = static_cast<int64_t>(ld[0]) * ld[1] * ld[2];
    const std::size_t np = static_cast<std::size_t>(np_i);

    // ── deterministic lattice spawn, row-major i + nx*(j + ny*k) ──
    std::vector<std::array<double, 3>> spawn_pos(np);
    for (int k = 0; k < ld[2]; ++k)
        for (int j = 0; j < ld[1]; ++j)
            for (int i = 0; i < ld[0]; ++i)
            {
                const std::size_t n = static_cast<std::size_t>(i)
                                      + static_cast<std::size_t>(ld[0])
                                          * (static_cast<std::size_t>(j)
                                             + static_cast<std::size_t>(ld[1])
                                                 * static_cast<std::size_t>(k));
                std::array<double, 3> p;
                p[0] = config.origin[0]
                       + (ld[0] > 1
                              ? static_cast<double>(i) / static_cast<double>(ld[0] - 1)
                              : 0.5)
                             * config.spawn_extent[0];
                p[1] = config.origin[1]
                       + (ld[1] > 1
                              ? static_cast<double>(j) / static_cast<double>(ld[1] - 1)
                              : 0.5)
                             * config.spawn_extent[1];
                p[2] = config.origin[2]
                       + (ld[2] > 1
                              ? static_cast<double>(k) / static_cast<double>(ld[2] - 1)
                              : 0.5)
                             * config.spawn_extent[2];
                spawn_pos[n] = p;
            }

    // ── one flat state vector: [pos p0, vel p0, pos p1, vel p1, ...] ──
    std::vector<double> state(6 * np, 0.0);
    for (std::size_t n = 0; n < np; ++n)
    {
        state[6 * n + 0] = spawn_pos[n][0];
        state[6 * n + 1] = spawn_pos[n][1];
        state[6 * n + 2] = spawn_pos[n][2];
        state[6 * n + 3] = config.initial_velocity[0];
        state[6 * n + 4] = config.initial_velocity[1];
        state[6 * n + 5] = config.initial_velocity[2];
    }

    bool warned_oob = false;
    const double k = config.drag_coefficient;
    solver::DerivativeFn deriv = [&](std::span<const double> st,
                                     std::span<double> dst, double)
    {
        for (std::size_t n = 0; n < np; ++n)
        {
            const std::array<double, 3> pos = {st[6 * n + 0], st[6 * n + 1], st[6 * n + 2]};
            const std::array<double, 3> vel = {st[6 * n + 3], st[6 * n + 4], st[6 * n + 5]};
            std::array<double, 3> vf = {0.0, 0.0, 0.0};
            if (config.flow_channel != nullptr)
            {
                std::array<double, 3> s = {0.0, 0.0, 0.0};
                if (config.flow_channel->sample(pos, s)) vf = s;
                else if (!warned_oob) warned_oob = true;
            }
            dst[6 * n + 0] = vel[0];
            dst[6 * n + 1] = vel[1];
            dst[6 * n + 2] = vel[2];
            dst[6 * n + 3] = config.gravity[0] - k * (vel[0] - vf[0]);
            dst[6 * n + 4] = config.gravity[1] - k * (vel[1] - vf[1]);
            dst[6 * n + 5] = config.gravity[2] - k * (vel[2] - vf[2]);
        }
    };

    auto record_probe = [&](const std::vector<double>& st, double t)
    {
        res.trajectory_probe.push_back({st[0], st[1], st[2]});
        res.velocity_probe.push_back({st[3], st[4], st[5]});
        res.time_history.push_back(t);
    };

    double t = 0.0;
    record_probe(state, t);
    for (uint64_t step = 0; step < config.max_steps; ++step)
    {
        exd::physics::ModelStatus ist;
        double dt_used = 0.0;
        const bool ok = solver::integrate_step(config.integration, t, config.dt,
                                               std::span<double>(state), deriv,
                                               ist, &dt_used);
        if (!ok || !ist.ok)
        {
            const std::string err = "particles: integration failed: "
                                    + (ist.error.empty() ? "integrator step rejected"
                                                         : ist.error);
            res.status = ModelStatus{false, err, warnings};
            status = res.status;
            return res;
        }
        t += (dt_used > 0.0) ? dt_used : config.dt;
        if ((step + 1) % config.history_interval == 0) record_probe(state, t);
    }
    record_probe(state, t);

    for (std::size_t i = 0; i < 6 * np; ++i)
    {
        if (!std::isfinite(state[i]))
        {
            res.status = ModelStatus{false, "particles: non-finite state after integration", warnings};
            status = res.status;
            return res;
        }
    }

    if (warned_oob)
        warnings.push_back("particles: flow channel sample out of bounds; v_f treated as 0 there");

    res.final_positions.resize(np);
    res.final_velocities.resize(np);
    double speed_sum = 0.0;
    for (std::size_t n = 0; n < np; ++n)
    {
        res.final_positions[n] = {state[6 * n + 0], state[6 * n + 1], state[6 * n + 2]};
        res.final_velocities[n] = {state[6 * n + 3], state[6 * n + 4], state[6 * n + 5]};
        const std::array<double, 3> v = res.final_velocities[n];
        speed_sum += std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    res.mean_speed = speed_sum / static_cast<double>(np);

    res.ok = true;
    res.status = ModelStatus{true, "", warnings};
    status = res.status;
    return res;
}

} // namespace exd::physics::particles