// Species transport — operator-split advection/diffusion/reaction on
// structured grids, composed from the core runtime (spec §34-35).

#include <exd/engine/physics/species/species_solver.hpp>

#include <exd/engine/discretization/fdm/operators.hpp>
#include <exd/engine/numerics/cg.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <cmath>
#include <span>

namespace exd::engine::physics::species {

using namespace exd::engine::core;
using exd::engine::coupling::IScalarField3D;
using exd::engine::coupling::IVectorField3D;
using exd::engine::discretization::fdm::DiffusionStepOperator;
using exd::engine::discretization::fdm::GhostSpec;
using exd::engine::discretization::fdm::FdmLaplacianOperator;
using exd::engine::numerics::IterativeSolverConfig;
using exd::engine::numerics::solve_cg;

namespace {

bool is_mirror(const GhostSpec& gs)
{
    for (auto& f : gs.faces)
        if (f.dirichlet) return false;
    return true;
}

GhostSpec make_ghosts(const SpeciesConfig& config)
{
    GhostSpec gs;
    for (auto& face : config.boundary_faces)
    {
        const size_t idx = static_cast<size_t>(face.face);
        gs.faces[idx].dirichlet = face.fixed;
        gs.faces[idx].value = face.value.empty() ? 0.0 : face.value.front();
    }
    return gs;
}

/// Pin Dirichlet face nodes to their fixed values (per species; uses the
/// first species value when a face spec carries a single value).
void pin_boundary(const SpeciesConfig& config, mesh::StructuredScalarGrid& c, size_t s)
{
    for (auto& face : config.boundary_faces)
    {
        if (!face.fixed) continue;
        const double v = face.value.empty() ? 0.0 :
                         face.value[std::min(s, face.value.size() - 1)];
        const int32_t nx = c.dims[0], ny = c.dims[1], nz = c.dims[2];
        switch (face.face)
        {
        case mesh::BoundaryId::XNeg:
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    c.values[static_cast<size_t>(nx) *
                             (static_cast<size_t>(j) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::XPos:
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    c.values[static_cast<size_t>(nx - 1) +
                             static_cast<size_t>(nx) * (static_cast<size_t>(j) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::YNeg:
            for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                    c.values[static_cast<size_t>(i) +
                             static_cast<size_t>(nx) * static_cast<size_t>(ny) * k] = v;
            break;
        case mesh::BoundaryId::YPos:
            for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                    c.values[static_cast<size_t>(i) +
                             static_cast<size_t>(nx) * (static_cast<size_t>(ny - 1) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::ZNeg:
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    c.values[static_cast<size_t>(i) + static_cast<size_t>(nx) * static_cast<size_t>(j)] = v;
            break;
        case mesh::BoundaryId::ZPos:
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    c.values[static_cast<size_t>(i) + static_cast<size_t>(nx) *
                             (static_cast<size_t>(j) + static_cast<size_t>(ny) * (static_cast<size_t>(nz - 1)))] = v;
            break;
        }
    }
}

} // anonymous namespace

bool validate_species_config(const SpeciesConfig& config, ModelStatus& status)
{
    if (!config.grid.validate(status)) return false;
    const size_t ns = config.species.size();
    if (ns == 0)
    {
        status.ok = false;
        status.error = "species: at least one species required";
        return false;
    }
    if (config.diffusivity.size() != ns || config.decay_rate.size() != ns ||
        config.initial_concentration.size() != ns)
    {
        status.ok = false;
        status.error = "species: diffusivity/decay_rate/initial_concentration must match species count";
        return false;
    }
    for (size_t s = 0; s < ns; ++s)
    {
        const auto& ic = config.initial_concentration[s];
        if (!ic.empty() && ic.size() != 1 && ic.size() != config.grid.node_count())
        {
            status.ok = false;
            status.error = "species: initial_concentration[" + std::to_string(s) +
                           "] must be 1 value (uniform) or the node count";
            return false;
        }
    }
    for (auto& conv : config.conversions)
    {
        if (conv.from >= ns || conv.to >= ns)
        {
            status.ok = false;
            status.error = "species: conversion index out of range";
            return false;
        }
        if (conv.from == conv.to)
        {
            status.ok = false;
            status.error = "species: conversion must be between distinct species";
            return false;
        }
    }
    for (auto& fv : config.boundary_faces)
    {
        if (!fv.value.empty() && fv.value.size() != ns)
        {
            status.ok = false;
            status.error = "species: boundary face value count must match species count (or be 1)";
            return false;
        }
    }
    return true;
}

bool step_species(const SpeciesConfig& config,
                  std::vector<mesh::StructuredScalarGrid>& state,
                  double dt, double& max_change_out,
                  ModelStatus& status)
{
    const size_t ns = config.species.size();
    if (state.size() != ns || !state[0].validate(status)) return false;

    const auto& g = config.grid;
    const size_t N = g.node_count();
    const GhostSpec gs = make_ghosts(config);
    const bool pure_local = is_mirror(gs);   // no boundary flux when all mirror

    // ── velocity at nodes: channel sample or body velocity ──
    std::vector<double> vx(N, config.body_velocity[0]);
    std::vector<double> vy(N, config.body_velocity[1]);
    std::vector<double> vz(N, config.body_velocity[2]);
    if (config.velocity_channel)
    {
        double umax = 0.0;
        for (int k = 0; k < g.dims[2]; ++k)
            for (int j = 0; j < g.dims[1]; ++j)
                for (int i = 0; i < g.dims[0]; ++i)
                {
                    std::array<double, 3> v{};
                    if (!config.velocity_channel->sample(g.node_coords(i, j, k), v))
                    {
                        status.ok = false;
                        status.error = "species: velocity channel out of bounds at a node";
                        return false;
                    }
                    const size_t idx = static_cast<size_t>(i) +
                        static_cast<size_t>(g.dims[0]) *
                        (static_cast<size_t>(j) + static_cast<size_t>(g.dims[1]) * k);
                    vx[idx] = v[0]; vy[idx] = v[1]; vz[idx] = v[2];
                    umax = std::max(umax, std::max(std::fabs(v[0]),
                             std::max(std::fabs(v[1]), std::fabs(v[2]))));
                }
        const double h_min = std::min({g.spacing[0], g.spacing[1], g.spacing[2]});
        if (umax > 0.0 && dt > 0.5 * h_min / umax)
            status.warnings.push_back("species: dt exceeds advective CFL 0.5·h/u; step clamped");
        if (umax > 0.0)
            dt = std::min(dt, 0.5 * h_min / umax);
    }

    max_change_out = 0.0;
    std::vector<double> adv(N);
    for (size_t s = 0; s < ns; ++s)
    {
        auto& c = state[s];
        const double D = config.diffusivity[s];
        const double k = config.decay_rate[s];

        // 1) explicit upwind advection (first-order, monotone)
        if (config.velocity_channel || config.body_velocity != std::array<double, 3>{0, 0, 0})
        {
            if (!exd::engine::discretization::fdm::upwind_advect(g, vx, vy, vz, c.values, adv, gs, status))
                return false;
            for (size_t i = 0; i < N; ++i)
                c.values[i] -= dt * adv[i];
        }

        // 2) implicit diffusion (I − dt·D·Δ) c_new = c_adv, CG on SPD −Δ
        if (D > 0.0)
        {
            Field rhs(FieldMetadata{
                .name = "rhs", .rank = FieldRank::Scalar, .components = 1,
                .location = FieldLocation::Node}, N);
            for (size_t i = 0; i < N; ++i) rhs.data()[i] = c.values[i];
            Field sol(FieldMetadata{
                .name = "sol", .rank = FieldRank::Scalar, .components = 1,
                .location = FieldLocation::Node}, N);
            std::copy(c.values.begin(), c.values.end(), sol.data().begin());
            DiffusionStepOperator op(g, gs, dt * D);
            IterativeSolverConfig cfg;
            cfg.tolerance = 1e-12;
            cfg.max_iterations = 800;
            auto rep = solve_cg(op, rhs, sol, cfg, status);
            if (!rep.converged && status.ok)
            {
                status.warnings.push_back("species: diffusion CG did not fully converge (" +
                                          std::to_string(rep.iterations) + " iters)");
            }
            if (!status.ok) return false;
            std::copy(sol.data().begin(), sol.data().end(), c.values.begin());
        }

        // 3) exact first-order decay
        if (k > 0.0)
        {
            const double f = std::exp(-k * dt);
            for (auto& v : c.values) v *= f;
        }

        // Dirichlet pins
        if (!pure_local) pin_boundary(config, c, s);
    }
    // pairs A → B (exact exponential pair; applied once per step, not per
    // species)
    for (auto& conv : config.conversions)
    {
        const double f = std::exp(-conv.rate * dt);
        for (size_t i = 0; i < N; ++i)
        {
            const double d = state[conv.from].values[i] * (1.0 - f);
            state[conv.from].values[i] -= d;
            state[conv.to].values[i] += d;
        }
    }
    (void)pure_local;
    return true;
}

SpeciesResult solve_species(const SpeciesConfig& config)
{
    SpeciesResult result;
    ModelStatus& status = result.status;
    if (!validate_species_config(config, status))
    {
        result.ok = false;
        return result;
    }
    const size_t ns = config.species.size();
    const size_t N = config.grid.node_count();

    result.concentration.resize(ns);
    for (size_t s = 0; s < ns; ++s)
    {
        auto& c = result.concentration[s];
        static_cast<mesh::StructuredGrid&>(c) = config.grid;   // base lattice
        const auto& ic = config.initial_concentration[s];
        if (ic.size() == 1)
            c.values.assign(N, ic.front());      // uniform initial value
        else if (ic.size() == N)
            c.values = ic;                        // per-node initial field
        else
        {
            status.ok = false;
            status.error = "species: initial_concentration must be 1 value (uniform) or the node count";
            result.ok = false;
            return result;
        }
        // apply Dirichlet pins to the initial state
        pin_boundary(config, c, s);
    }

    const uint64_t max_steps = config.max_steps;
    for (uint64_t it = 0; it < max_steps; ++it)
    {
        double change = 0.0;
        double dt = config.dt;
        // measure pre-step state for the change norm
        std::vector<std::vector<double>> before(ns);
        for (size_t s = 0; s < ns; ++s) before[s] = result.concentration[s].values;
        if (!step_species(config, result.concentration, dt, change, status))
        {
            result.ok = false;
            return result;
        }
        result.steps = it + 1;
        result.time += dt;
        result.max_change = 0.0;
        for (size_t s = 0; s < ns; ++s)
            for (size_t i = 0; i < N; ++i)
                result.max_change = std::max(result.max_change,
                    std::fabs(result.concentration[s].values[i] - before[s][i]));
        if (config.steady && result.max_change < config.steady_tolerance)
            break;
        if (result.time > 1e18) break;   // safety
    }
    if (result.steps >= max_steps && config.steady && result.max_change >= config.steady_tolerance)
        status.warnings.push_back("species: steady tolerance not reached within max_steps");
    if (!config.steady && result.steps >= max_steps)
        status.warnings.push_back("species: max_steps reached");

    result.total_mass.resize(ns, 0.0);
    for (size_t s = 0; s < ns; ++s)
        for (double v : result.concentration[s].values)
            result.total_mass[s] += v;
    result.ok = status.ok;
    return result;
}

std::unique_ptr<IScalarField3D> make_concentration_channel(
    const std::vector<mesh::StructuredScalarGrid>& state, size_t s)
{
    if (s >= state.size()) return nullptr;
    return exd::engine::coupling::make_scalar_grid_field(state[s]);
}

} // namespace exd::engine::physics::species
