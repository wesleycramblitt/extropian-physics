// Porous media — implicit pressure diffusion composed from the core runtime
// (spec §34-35): DiffusionStepOperator + CG.

#include <exd/engine/physics/porous/porous_solver.hpp>

#include <exd/engine/discretization/fdm/operators.hpp>
#include <exd/engine/numerics/cg.hpp>

#include <cmath>

namespace exd::engine::physics::porous {

using namespace exd::engine::core;
using exd::engine::discretization::fdm::DiffusionStepOperator;
using exd::engine::discretization::fdm::FdmLaplacianOperator;
using exd::engine::discretization::fdm::GhostSpec;
using exd::engine::numerics::IterativeSolverConfig;
using exd::engine::numerics::NegatedOperator;
using exd::engine::numerics::ScaledOperator;
using exd::engine::numerics::solve_cg;

namespace {

/// Linear part L(x) = A(x) − A(0) of an AFFINE operator A (the
/// eliminated-Dirichlet operators are affine in the fixed face values:
/// A(x) = M·x + c).  Iterative solvers require the linear part.
class AffineToLinear final : public IOperator
{
public:
    AffineToLinear(const IOperator& inner, size_t n)
        : inner_(inner), fixed_(FieldMetadata{
              .name = "affine_const", .rank = FieldRank::Scalar, .components = 1,
              .location = FieldLocation::Node}, n)
    {
        zero_ = fixed_;
        zero_.assign(0.0);
        inner_.apply(zero_, fixed_, fixed_status_);   // fixed_ = A(0)
    }
    const OperatorInfo& info() const override { return inner_.info(); }
    bool apply(const Field& in, Field& out, ModelStatus& status) const override
    {
        if (!inner_.apply(in, out, status)) return false;
        for (size_t i = 0; i < out.size(); ++i)
            out.data()[i] -= fixed_.data()[i];   // L(x) = A(x) − A(0)
        return true;
    }
    bool apply_transpose(const Field& in, Field& out, ModelStatus& status) const override
    {
        if (status.ok /* inner transpose unsupported */) {}
        status.ok = false;
        status.error = "affine-to-linear: transpose not supported";
        return false;
    }
    bool diagonal(Field& out, ModelStatus& status) const override
    {
        return inner_.diagonal(out, status);   // A(0) contributes no diagonal
    }
    bool jacobian_vector_product(const Field& x, const Field& v,
                                 Field& out, ModelStatus& status) const override
    {
        return apply(v, out, status);   // linear part: J(x)·v = L(v)
    }

private:
    const IOperator& inner_;
    mutable Field zero_;
    Field fixed_;
    ModelStatus fixed_status_;
};

GhostSpec make_ghosts(const PorousConfig& config)
{
    GhostSpec gs;
    for (auto& face : config.boundary_faces)
    {
        const size_t idx = static_cast<size_t>(face.face);
        gs.faces[idx].dirichlet = face.fixed;
        gs.faces[idx].value = face.value;
    }
    return gs;
}

void pin_faces(const PorousConfig& config, mesh::StructuredScalarGrid& p)
{
    const int32_t nx = p.dims[0], ny = p.dims[1], nz = p.dims[2];
    for (auto& face : config.boundary_faces)
    {
        if (!face.fixed) continue;
        const double v = face.value;
        switch (face.face)
        {
        case mesh::BoundaryId::XNeg:
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    p.values[static_cast<size_t>(nx) *
                             (static_cast<size_t>(j) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::XPos:
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    p.values[static_cast<size_t>(nx - 1) +
                             static_cast<size_t>(nx) *
                                 (static_cast<size_t>(j) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::YNeg:
            for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                    p.values[static_cast<size_t>(i) +
                             static_cast<size_t>(nx) * static_cast<size_t>(ny) * k] = v;
            break;
        case mesh::BoundaryId::YPos:
            for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                    p.values[static_cast<size_t>(i) +
                             static_cast<size_t>(nx) *
                                 (static_cast<size_t>(ny - 1) + static_cast<size_t>(ny) * k)] = v;
            break;
        case mesh::BoundaryId::ZNeg:
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    p.values[static_cast<size_t>(i) +
                             static_cast<size_t>(nx) * static_cast<size_t>(j)] = v;
            break;
        case mesh::BoundaryId::ZPos:
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    p.values[static_cast<size_t>(i) + static_cast<size_t>(nx) *
                             (static_cast<size_t>(j) +
                              static_cast<size_t>(ny) * (static_cast<size_t>(nz - 1)))] = v;
            break;
        }
    }
}

} // namespace

bool validate_porous_config(const PorousConfig& config, ModelStatus& status)
{
    if (!config.grid.validate(status)) return false;
    if (!(config.porosity > 0.0) || !(config.compressibility > 0.0) ||
        !(config.viscosity > 0.0) || !(config.permeability > 0.0))
    {
        status.ok = false;
        status.error = "porous: permeability/viscosity/porosity/compressibility must be > 0";
        return false;
    }
    if (!config.source_rate.empty() && config.source_rate.size() != config.grid.node_count())
    {
        status.ok = false;
        status.error = "porous: source_rate must be empty or the node count";
        return false;
    }
    return true;
}

PorousResult solve_porous(const PorousConfig& config)
{
    PorousResult result;
    ModelStatus& status = result.status;
    if (!validate_porous_config(config, status))
    {
        result.ok = false;
        return result;
    }

    const size_t N = config.grid.node_count();
    if (!config.initial_pressure_field.empty() &&
        config.initial_pressure_field.size() != N)
    {
        status.ok = false;
        status.error = "porous: initial_pressure_field must be empty or the node count";
        result.ok = false;
        return result;
    }
    static_cast<mesh::StructuredGrid&>(result.pressure) = config.grid;
    if (!config.initial_pressure_field.empty())
        result.pressure.values = config.initial_pressure_field;
    else
        result.pressure.values.assign(N, config.initial_pressure);
    pin_faces(config, result.pressure);

    const double K = hydraulic_diffusivity(config);
    const GhostSpec gs = make_ghosts(config);
    bool has_source = !config.source_rate.empty();
    const double cell_vol = config.grid.cell_volume();

    // ── steady: DIRECT one-shot solve −K·Δp = q with the boundary values
    //    transferred to the RHS (eliminated-Dirichlet operator, SPD, spec
    //    §35 matrix-free).  No pseudo-transient. ──
    if (config.steady)
    {
        FdmLaplacianOperator lap(config.grid, gs);
        NegatedOperator neg(lap);
        ScaledOperator A(neg, K);          // A = −K·Δ (eliminated Dirichlet)
        // A is AFFINE: A(x) = M·x + c with the fixed face values inside c.
        // The steady state satisfies A(p0 + w) = 0 → M·w = −A(p0), where
        // p0 carries the face values.  Solve the LINEAR part L = A − A(0)
        // with b = −A(p0) (evaluated with the affine operator).
        AffineToLinear L(A, N);
        mesh::StructuredScalarGrid p0;
        static_cast<mesh::StructuredGrid&>(p0) = config.grid;
        p0.values.assign(N, 0.0);
        pin_faces(config, p0);             // faces = g, interior = 0
        Field p0_f(FieldMetadata{
            .name = "p0", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node}, N);
        std::copy(p0.values.begin(), p0.values.end(), p0_f.data().begin());
        Field rhs(FieldMetadata{
            .name = "rhs", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node}, N);
        if (!A.apply(p0_f, rhs, status)) return result;
        for (size_t i = 0; i < N; ++i) rhs.data()[i] *= -1.0;   // b = −A(p0)
        if (has_source)
            for (size_t i = 0; i < N; ++i)
                rhs.data()[i] += config.source_rate[i];
        Field sol(FieldMetadata{
            .name = "sol", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node}, N);
        sol.assign(0.0);                   // w starts at 0
        IterativeSolverConfig cfg;
        cfg.tolerance = 1e-10;
        cfg.max_iterations = 20000;
        auto rep = solve_cg(L, rhs, sol, cfg, status);
        if (!rep.converged && status.ok)
            status.warnings.push_back("porous: steady CG did not fully converge");
        if (!status.ok)
        {
            result.ok = false;
            return result;
        }
        for (size_t i = 0; i < N; ++i)
            result.pressure.values[i] = p0_f.data()[i] + sol.data()[i];
        pin_faces(config, result.pressure);
        result.steps = 1;
        result.max_change = rep.final_residual;
        result.time = 0.0;
        result.total_mass = 0.0;
        for (double p : result.pressure.values)
            result.total_mass += p * config.porosity * config.compressibility * cell_vol;
        result.ok = true;
        result.status = ModelStatus{true, "", status.warnings};
        return result;
    }

    std::vector<double> before(N);
    for (uint64_t it = 0; it < config.max_steps; ++it)
    {
        before = result.pressure.values;
        // implicit pressure diffusion: (I − dt·K·Δ)p^{n+1} = p^n + dt·q
        Field rhs(FieldMetadata{
            .name = "rhs", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node}, N);
        for (size_t i = 0; i < N; ++i)
        {
            rhs.data()[i] = result.pressure.values[i];
            if (has_source) rhs.data()[i] += config.dt * config.source_rate[i];
        }
        Field sol(FieldMetadata{
            .name = "sol", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node}, N);
        std::copy(result.pressure.values.begin(), result.pressure.values.end(),
                  sol.data().begin());
        DiffusionStepOperator op(config.grid, gs, config.dt * K);
        IterativeSolverConfig cfg;
        cfg.tolerance = 1e-12;
        cfg.max_iterations = 800;
        auto rep = solve_cg(op, rhs, sol, cfg, status);
        if (!rep.converged && status.ok)
            status.warnings.push_back("porous: diffusion CG did not fully converge");
        if (!status.ok)
        {
            result.ok = false;
            return result;
        }
        std::copy(sol.data().begin(), sol.data().end(), result.pressure.values.begin());
        pin_faces(config, result.pressure);

        result.steps = it + 1;
        result.time += config.dt;
        result.max_change = 0.0;
        for (size_t i = 0; i < N; ++i)
            result.max_change = std::max(result.max_change,
                                         std::fabs(result.pressure.values[i] - before[i]));
        if (config.steady && result.max_change < config.steady_tolerance)
            break;
    }
    if (config.steady && result.max_change >= config.steady_tolerance)
        status.warnings.push_back("porous: steady tolerance not reached within max_steps");

    result.total_mass = 0.0;
    for (double p : result.pressure.values)
        result.total_mass += p * config.porosity * config.compressibility * cell_vol;
    result.ok = status.ok;
    return result;
}

std::unique_ptr<coupling::IScalarField3D> make_pressure_channel(const PorousResult& result)
{
    return coupling::make_scalar_grid_field(result.pressure);
}

std::unique_ptr<coupling::IVectorField3D> make_darcy_velocity_channel(
    const PorousConfig& config, const PorousResult& result)
{
    if (!result.ok) return nullptr;
    const size_t N = config.grid.node_count();
    coupling::StructuredVectorGrid v;
    v.origin = config.grid.origin;
    v.spacing = config.grid.spacing;
    v.dims = config.grid.dims;
    v.values.assign(3 * N, 0.0);

    const double coef = config.permeability / config.viscosity;
    const auto& g = config.grid;
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                auto at = [&](int32_t ii, int32_t jj, int32_t kk) {
                    const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
                    const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
                    const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
                    return result.pressure.values[static_cast<size_t>(im) +
                        static_cast<size_t>(nx) * (static_cast<size_t>(jm) +
                        static_cast<size_t>(ny) * km)];
                };
                const size_t idx = static_cast<size_t>(i) +
                    static_cast<size_t>(nx) * (static_cast<size_t>(j) +
                    static_cast<size_t>(ny) * k);
                const double dp_dx = (at(i + 1, j, k) - at(i - 1, j, k)) / (2 * g.spacing[0]);
                const double dp_dy = (at(i, j + 1, k) - at(i, j - 1, k)) / (2 * g.spacing[1]);
                const double dp_dz = (at(i, j, k + 1) - at(i, j, k - 1)) / (2 * g.spacing[2]);
                v.values[3 * idx + 0] = -coef * dp_dx;
                v.values[3 * idx + 1] = -coef * dp_dy;
                v.values[3 * idx + 2] = -coef * dp_dz;
            }
    return coupling::make_vector_grid_field(v);
}

} // namespace exd::engine::physics::porous
