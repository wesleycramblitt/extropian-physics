#include <exd/engine/presets/multiphysics/species_in_flow.hpp>

#include <exd/engine/coupling/field_channels.hpp>

#include <cmath>

namespace exd::engine::presets::multiphysics {

SpeciesInFlowResult run_species_in_flow(const SpeciesInFlowConfig& config)
{
    SpeciesInFlowResult result;
    core::ModelStatus& status = result.status;

    // ── steady duct flow ──
    auto flow_cfg = presets::cfd::make_channel_flow(config.flow, status);
    if (!status.ok) return result;
    result.flow_result = exd::engine::physics::fluid::fdm3::solve_fdm3(flow_cfg);
    if (!result.flow_result.valid)
    {
        status.ok = false;
        status.error = "species_in_flow: fdm3 failed: " + result.flow_result.error;
        return result;
    }

    // ── velocity channel over the solved field (cell-center lattice) ──
    const auto& f = result.flow_result.field;
    exd::engine::coupling::StructuredVectorGrid vel;
    vel.origin = {0.025, 0.025, 0.025};
    vel.spacing = {0.05, 0.05, 0.05};
    vel.dims = {f.nx, f.ny, f.nz};
    vel.values.resize(3ull * static_cast<size_t>(f.nx * f.ny * f.nz));
    for (int k = 0; k < f.nz; ++k)
        for (int j = 0; j < f.ny; ++j)
            for (int i = 0; i < f.nx; ++i)
            {
                const size_t idx = f.index(i, j, k);
                vel.values[3 * idx + 0] = f.u[idx];
                vel.values[3 * idx + 1] = f.v[idx];
                vel.values[3 * idx + 2] = f.w[idx];
            }
    auto velocity_channel = exd::engine::coupling::make_vector_grid_field(vel);
    if (!velocity_channel)
    {
        status.ok = false;
        status.error = "species_in_flow: could not build velocity channel";
        return result;
    }

    // ── species: advection-decay on the same lattice ──
    exd::engine::physics::species::SpeciesConfig scfg;
    scfg.grid.origin = {0.025, 0.025, 0.025};
    scfg.grid.spacing = {0.05, 0.05, 0.05};
    scfg.grid.dims = {f.nx, f.ny, f.nz};
    scfg.species = {"c"};
    scfg.diffusivity = {0.0};
    scfg.decay_rate = {config.decay_rate};
    scfg.initial_concentration = {{0.0}};
    scfg.velocity_channel = velocity_channel.get();
    scfg.dt = config.species_dt;
    scfg.max_steps = config.max_steps;
    scfg.steady = true;
    scfg.steady_tolerance = config.steady_tolerance;
    scfg.boundary_faces = {
        {exd::engine::mesh::BoundaryId::XNeg, true, {config.inlet_concentration}},
    };
    result.species = exd::engine::physics::species::solve_species(scfg);
    if (!result.species.ok)
    {
        status.ok = false;
        status.error = "species_in_flow: species solve failed";
        return result;
    }

    // outlet plane (cell-center lattice: last plane = nx-1): FLUX-weighted
    // concentration (the slow near-wall lanes decay more; the flux weighting
    // reflects what actually leaves the domain)
    const size_t plane = static_cast<size_t>(f.nx - 1);
    double c_flux = 0.0, u_flux = 0.0;
    for (int k = 0; k < f.nz; ++k)
        for (int j = 0; j < f.ny; ++j)
        {
            const size_t idx = plane + static_cast<size_t>(f.nx) *
                                   (static_cast<size_t>(j) + static_cast<size_t>(f.ny) * k);
            c_flux += f.u[idx] * result.species.concentration[0].values[idx];
            u_flux += f.u[idx];
        }
    result.outlet_concentration = u_flux > 0.0 ? c_flux / u_flux : 0.0;
    result.analytic_outlet = config.inlet_concentration *
                             std::exp(-config.decay_rate * config.flow.length /
                                      config.flow.inlet_velocity);
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::multiphysics
