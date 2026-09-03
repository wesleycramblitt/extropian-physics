#include <exd/engine/presets/multiphysics/joule_heating.hpp>

#include <exd/engine/coupling/field_channels.hpp>

namespace exd::engine::presets::multiphysics {

JouleHeatingResult run_joule_heating(const JouleHeatingConfig& config)
{
    JouleHeatingResult result;
    core::ModelStatus& status = result.status;

    // ── electrostatics: potential between the x faces ──
    exd::engine::physics::electromagnetics::StaticFieldConfig ecfg;
    ecfg.mode = exd::engine::physics::electromagnetics::StaticFieldMode::Electrostatic;
    ecfg.dims = {config.nx, config.ny, config.nz};
    ecfg.spacing = {config.spacing, config.spacing, config.spacing};
    // Parallel-plate capacitor: the plates ARE the x faces; the side walls
    // are NEUMANN (zero normal gradient), so the field is the EXACT linear
    // bridge with E = −V/L uniform — the discrete solution the module's
    // Neumann-face support guarantees (see static_fields_test).
    ecfg.face_values = {0.0, config.voltage, 0.0, 0.0, 0.0, 0.0};
    for (int a = 0; a < 6; ++a)
        ecfg.face_kind[static_cast<size_t>(a)] = (a < 2)
            ? exd::engine::physics::electromagnetics::FaceKind::Dirichlet
            : exd::engine::physics::electromagnetics::FaceKind::Neumann;
    const double gap_d = (config.nx - 1) * config.spacing;   // full x extent
    ecfg.sor_omega = 1.8;
    ecfg.tolerance = 1e-10;
    ecfg.max_iterations = 200000;
    result.electric = exd::engine::physics::electromagnetics::solve_static_field(ecfg);
    if (!result.electric.ok)
    {
        status.ok = false;
        status.error = "joule_heating: static field failed: " + result.electric.error;
        return result;
    }

    // ── per-node heat source grid q = σ|E|² on the thermal lattice ──
    // E = −∇φ computed directly from the solved POTENTIAL field (central
    // differences on the interior grid; zero-gradient beyond the faces).
    const auto& pot = result.electric.potential;
    exd::engine::coupling::StructuredScalarGrid source_grid;
    source_grid.origin = pot.origin;
    source_grid.spacing = pot.spacing;
    source_grid.dims = pot.dims;
    source_grid.values.assign(static_cast<size_t>(config.nx) * config.ny * config.nz, 0.0);
    auto phi_at = [&](int32_t i, int32_t j, int32_t k) -> double {
        const int32_t im = (i < 0) ? 0 : (i >= pot.dims[0] ? pot.dims[0] - 1 : i);
        const int32_t jm = (j < 0) ? 0 : (j >= pot.dims[1] ? pot.dims[1] - 1 : j);
        const int32_t km = (k < 0) ? 0 : (k >= pot.dims[2] ? pot.dims[2] - 1 : k);
        return pot.values[static_cast<size_t>(im) +
                          static_cast<size_t>(pot.dims[0]) *
                              (static_cast<size_t>(jm) +
                               static_cast<size_t>(pot.dims[1]) * km)];
    };
    double power_sum = 0.0;
    const double cell_vol = config.spacing * config.spacing * config.spacing;
    for (int32_t k = 0; k < pot.dims[2]; ++k)
        for (int32_t j = 0; j < pot.dims[1]; ++j)
            for (int32_t i = 0; i < pot.dims[0]; ++i)
            {
                const double ex = -(phi_at(i + 1, j, k) - phi_at(i - 1, j, k)) /
                                  (2.0 * config.spacing);
                const double ey = -(phi_at(i, j + 1, k) - phi_at(i, j - 1, k)) /
                                  (2.0 * config.spacing);
                const double ez = -(phi_at(i, j, k + 1) - phi_at(i, j, k - 1)) /
                                  (2.0 * config.spacing);
                const double e2 = ex * ex + ey * ey + ez * ez;
                const double q = config.conductivity * e2;
                source_grid.values[static_cast<size_t>(i) +
                                   static_cast<size_t>(pot.dims[0]) *
                                       (static_cast<size_t>(j) +
                                        static_cast<size_t>(pot.dims[1]) * k)] = q;
                // cell convention (matching the thermal module's integral):
                // q sampled at the cell MIN corner over (nx−1)(ny−1)(nz−1)
                // cells
                if (i < pot.dims[0] - 1 && j < pot.dims[1] - 1 && k < pot.dims[2] - 1)
                    power_sum += q * cell_vol;
            }
    result.total_power = power_sum;
    // parallel-plate anchor: P = σ·V²·A/d with A the full cross-section and
    // d the uniform-gap separation (plates inset one cell from each face)
    const double d = gap_d;
    const double area = (config.ny - 1) * (config.nz - 1) * config.spacing * config.spacing;
    result.analytic_power = config.conductivity * config.voltage * config.voltage *
                            area / d;
    auto source_channel = exd::engine::coupling::make_scalar_grid_field(source_grid);
    if (!source_channel)
    {
        status.ok = false;
        status.error = "joule_heating: could not build q channel";
        return result;
    }

    // ── thermal: fixed T at both x faces, insulated elsewhere, q from the
    //    source channel ──
    exd::engine::physics::thermal::ThermalConfig tcfg;
    // SAME lattice as the electrostatic field (origin/spacing/dims) so the
    // q source channel aligns node-for-node
    tcfg.grid.origin = result.electric.potential.origin;
    tcfg.grid.spacing = {config.spacing, config.spacing, config.spacing};
    tcfg.grid.dims = {config.nx, config.ny, config.nz};
    tcfg.material.conductivity = config.thermal_conductivity;
    tcfg.material.density = 1.0;
    tcfg.material.specific_heat = 1.0;
    tcfg.boundary_values = {config.t_wall, config.t_wall, 0, 0, 0, 0};
    for (int f = 0; f < 6; ++f)
        tcfg.boundary_kind[static_cast<size_t>(f)] =
            exd::engine::physics::thermal::ThermalBoundaryKind::Insulated;
    tcfg.boundary_kind[0] = exd::engine::physics::thermal::ThermalBoundaryKind::FixedValue;
    tcfg.boundary_kind[1] = exd::engine::physics::thermal::ThermalBoundaryKind::FixedValue;
    tcfg.source_channel = source_channel.get();
    tcfg.tolerance = config.tolerance;
    result.thermal = exd::engine::physics::thermal::solve_thermal(tcfg, status);
    if (!result.thermal.ok)
    {
        status.ok = false;
        status.error = "joule_heating: thermal solve failed";
        return result;
    }

    // mid temperature and center heat rate (the local-q parabola anchor)
    const size_t mid = static_cast<size_t>(config.nx / 2) +
                       static_cast<size_t>(config.nx) *
                           (static_cast<size_t>(config.ny / 2) +
                            static_cast<size_t>(config.ny) * (config.nz / 2));
    result.mid_temperature = result.thermal.temperature.values[mid];
    result.center_source = source_grid.values[mid];
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::multiphysics
