#include <exd/engine/presets/multiphysics/aeroacoustics.hpp>

#include <cmath>

namespace exd::engine::presets::multiphysics {

AeroacousticsResult run_aeroacoustics(const AeroacousticsConfig& config)
{
    AeroacousticsResult result;
    core::ModelStatus& status = result.status;

    // ── flow: steady duct (presets/cfd) ──
    auto flow_cfg = presets::cfd::make_channel_flow(config.flow, status);
    if (!status.ok) return result;
    result.flow_result = exd::engine::physics::fluid::fdm3::solve_fdm3(flow_cfg);
    if (!result.flow_result.valid)
    {
        status.ok = false;
        status.error = "aeroacoustics: fdm3 flow failed: " + result.flow_result.error;
        return result;
    }

    // ── cross-section mean velocity (mass-conservation value = inlet) ──
    const auto& f = result.flow_result.field;
    double u_sum = 0.0;
    for (size_t i = 0; i < static_cast<size_t>(f.nx); ++i)
        for (size_t j = 0; j < static_cast<size_t>(f.ny); ++j)
            for (size_t k = 0; k < static_cast<size_t>(f.nz); ++k)
                u_sum += f.u[f.index(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k))];
    const size_t n_cells = static_cast<size_t>(f.nx * f.ny * f.nz);
    result.mean_velocity = {u_sum / static_cast<double>(n_cells), 0.0, 0.0};

    // ── acoustics on the same lattice with the mean flow ──
    exd::engine::physics::acoustics::WaveConfig wave;
    wave.grid.origin = {0.0, 0.0, 0.0};
    wave.grid.spacing = flow_cfg.lx / flow_cfg.nx == 0.0
                            ? std::array<double, 3>{0.05, 0.05, 0.05}
                            : std::array<double, 3>{flow_cfg.lx / flow_cfg.nx,
                                                    flow_cfg.ly / flow_cfg.ny,
                                                    flow_cfg.lz / flow_cfg.nz};
    wave.grid.dims = {flow_cfg.nx, flow_cfg.ny, flow_cfg.nz};
    wave.sound_speed = config.sound_speed;
    wave.mean_flow = config.include_flow ? result.mean_velocity
                                         : std::array<double, 3>{0.0, 0.0, 0.0};
    wave.max_steps = config.max_steps;
    wave.initial_pressure.assign(static_cast<size_t>(wave.grid.dims[0]) *
                                     wave.grid.dims[1] * wave.grid.dims[2], 0.0);
    const double hx = wave.grid.spacing[0];
    for (int32_t k = 0; k < wave.grid.dims[2]; ++k)
        for (int32_t j = 0; j < wave.grid.dims[1]; ++j)
            for (int32_t i = 0; i < wave.grid.dims[0]; ++i)
            {
                const double x = static_cast<double>(i) * hx;
                const double r = (x - config.pulse_center[0]) / config.pulse_width;
                wave.initial_pressure[static_cast<size_t>(i) +
                                      static_cast<size_t>(wave.grid.dims[0]) *
                                          (static_cast<size_t>(j) +
                                           static_cast<size_t>(wave.grid.dims[1]) * k)] =
                    config.amplitude * std::exp(-0.5 * r * r);
            }
    if (config.probe_index >= 0)
        wave.probe_index = config.probe_index;
    else
    {
        const size_t wc = static_cast<size_t>(wave.grid.dims[0]) * wave.grid.dims[1] *
                          wave.grid.dims[2];
        wave.probe_index = static_cast<int32_t>(wc / 2);
    }

    result.wave_result = exd::engine::physics::acoustics::solve_wave(wave, status);
    if (!result.wave_result.ok)
    {
        status.ok = false;
        status.error = "aeroacoustics: wave solve failed";
        return result;
    }
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::multiphysics
