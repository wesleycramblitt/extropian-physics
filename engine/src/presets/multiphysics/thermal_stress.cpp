#include <exd/engine/presets/multiphysics/thermal_stress.hpp>

#include <exd/engine/coupling/field_channels.hpp>

namespace exd::engine::presets::multiphysics {

ThermalStressResult run_thermal_stress(const ThermalStressConfig& config)
{
    ThermalStressResult result;
    core::ModelStatus& status = result.status;

    // ── thermal: steady conduction slab, linear T ──
    exd::engine::physics::thermal::ThermalConfig tcfg;
    tcfg.grid.origin = {0, 0, 0};
    tcfg.grid.spacing = {config.length / (config.nx - 1),
                         config.length / (config.nx - 1),
                         config.length / (config.nx - 1)};
    tcfg.grid.dims = {config.nx, config.ny, config.nz};
    tcfg.material.conductivity = config.conductivity;
    tcfg.material.density = 1.0;
    tcfg.material.specific_heat = 1.0;
    tcfg.boundary_values = {config.t_right, config.t_left, 0, 0, 0, 0};
    for (int f = 0; f < 6; ++f)
        tcfg.boundary_kind[static_cast<size_t>(f)] =
            exd::engine::physics::thermal::ThermalBoundaryKind::Insulated;
    tcfg.boundary_kind[0] = exd::engine::physics::thermal::ThermalBoundaryKind::FixedValue;
    tcfg.boundary_kind[1] = exd::engine::physics::thermal::ThermalBoundaryKind::FixedValue;
    tcfg.tolerance = 1e-8;
    result.thermal = exd::engine::physics::thermal::solve_thermal(tcfg, status);
    if (!result.thermal.ok)
    {
        status.ok = false;
        status.error = "thermal_stress: thermal solve failed";
        return result;
    }

    // temperature channel over the solved field
    auto temp_channel = exd::engine::coupling::make_scalar_grid_field(result.thermal.temperature);
    if (!temp_channel)
    {
        status.ok = false;
        status.error = "thermal_stress: could not build temperature channel";
        return result;
    }

    // ── structural: free bar, axial roller pin at the x− center node ──
    exd::engine::physics::structural::ElasticityConfig ecfg;
    ecfg.grid.origin = {0, 0, 0};
    ecfg.grid.spacing = tcfg.grid.spacing;
    ecfg.grid.dims = {config.nx, config.ny, config.nz};
    ecfg.material.elastic_modulus = config.elastic_modulus;
    ecfg.material.poisson_ratio = config.poisson_ratio;
    ecfg.thermal_expansion_coefficient = config.thermal_expansion;
    ecfg.tolerance = config.tolerance;
    ecfg.max_iterations = config.max_iterations;

    const size_t N = static_cast<size_t>(config.nx) * config.ny * config.nz;
    const auto idx3 = [&](int i, int j, int k) {
        return static_cast<size_t>(i) + static_cast<size_t>(config.nx) *
                   (static_cast<size_t>(j) + static_cast<size_t>(config.ny) * k);
    };
    ecfg.fixed_mask.assign(N, {false, false, false});
    // rigid-mode killer (translation + the three rotations) while leaving the
    // free axial expansion unconstrained: full pin at (0,0,0), rollers at
    // (nx-1,0,0) u_y/u_z and (0,ny-1,0) u_z.
    const size_t pin0 = idx3(0, 0, 0);
    ecfg.fixed_mask[pin0] = {true, true, true};
    const size_t pin1 = idx3(config.nx - 1, 0, 0);
    ecfg.fixed_mask[pin1][1] = true;
    ecfg.fixed_mask[pin1][2] = true;
    const size_t pin2 = idx3(0, config.ny - 1, 0);
    ecfg.fixed_mask[pin2][2] = true;

    result.structural = exd::engine::physics::structural::solve_elasticity(
        ecfg, status, temp_channel.get());
    if (!result.structural.ok)
    {
        status.ok = false;
        status.error = "thermal_stress: structural solve failed";
        return result;
    }

    // tip displacement on the center line (i = nx-1, j = ny/2, k = nz/2)
    const size_t tip = static_cast<size_t>(config.nx - 1) +
                       static_cast<size_t>(config.nx) *
                           (static_cast<size_t>(config.ny / 2) +
                            static_cast<size_t>(config.ny) * (config.nz / 2));
    result.measured_tip_displacement = result.structural.displacement.values[3 * tip + 0];
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::multiphysics
