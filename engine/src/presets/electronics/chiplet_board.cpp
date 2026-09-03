// Chiplet-board preset — heterogeneous steady conduction with per-node
// conductivity and sources, composed from the core operators (affine
// Dirichlet operator + CG via the shared AffineLinearPart).

#include <exd/engine/presets/electronics/chiplet_board.hpp>

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/operator.hpp>
#include <exd/engine/numerics/cg.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <algorithm>
#include <cmath>

namespace exd::engine::presets::electronics {

namespace {

using exd::engine::core::Field;
using exd::engine::core::FieldMetadata;
using exd::engine::core::FieldRank;
using exd::engine::core::FieldLocation;
using exd::engine::numerics::AffineLinearPart;
using exd::engine::numerics::IterativeSolverConfig;
using exd::engine::numerics::solve_cg;

/// The steady heterogeneous conduction operator:
///   (L·T)_i = Σ_axes (kf⁺·(T⁺ − T) − kf⁻·(T − T⁻))/h²
/// with the face conductivity kf = harmonic mean of the two node values
/// (exact for series conduction), and DIRICHLET faces (the sink).
class ConductionOperator final : public exd::engine::core::IOperator
{
public:
    ConductionOperator(mesh::StructuredGrid g, const std::vector<double>& k,
                       double sink_temp)
        : g_(std::move(g)), k_(k), sink_(sink_temp)
    {
        info_.name = "chiplet_conduction";
        info_.inputs.push_back(FieldMetadata{
            .name = "T", .rank = FieldRank::Scalar, .components = 1,
            .location = FieldLocation::Node});
        info_.outputs.push_back(info_.inputs.front());
    }
    const exd::engine::core::OperatorInfo& info() const override { return info_; }

    bool apply(const Field& in, Field& out, exd::engine::core::ModelStatus& status) const override
    {
        const int32_t nx = g_.dims[0], ny = g_.dims[1], nz = g_.dims[2];
        const size_t N = static_cast<size_t>(nx) * ny * nz;
        if (in.size() != N || out.size() != N)
        {
            status.ok = false;
            status.error = "chiplet conduction: field size mismatch";
            return false;
        }
        auto idx = [&](int i, int j, int k) {
            return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                       (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
        };
        auto T = [&](int i, int j, int k) -> double {
            // Dirichlet faces: the sink value; interior: the field
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return sink_;
            return in.data()[idx(i, j, k)];
        };
        auto kf = [&](double ka, double kb) -> double {
            // harmonic mean (series conduction); equal values → the value
            return (ka + kb > 0.0) ? 2.0 * ka * kb / (ka + kb) : 0.0;
        };
        auto kn = [&](int i, int j, int k) -> double {
            const int32_t im = (i < 0) ? 0 : (i >= nx ? nx - 1 : i);
            const int32_t jm = (j < 0) ? 0 : (j >= ny ? ny - 1 : j);
            const int32_t km = (k < 0) ? 0 : (k >= nz ? nz - 1 : k);
            return k_[idx(im, jm, km)];
        };
        const double hx = g_.spacing[0], hy = g_.spacing[1], hz = g_.spacing[2];
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const size_t id = idx(i, j, k);
                    const double txp = (i < nx - 1) ? in.data()[idx(i + 1, j, k)] : sink_;
                    const double txm = (i > 0) ? in.data()[idx(i - 1, j, k)] : sink_;
                    const double typ = (j < ny - 1) ? in.data()[idx(i, j + 1, k)] : sink_;
                    const double tym = (j > 0) ? in.data()[idx(i, j - 1, k)] : sink_;
                    const double tzp = (k < nz - 1) ? in.data()[idx(i, j, k + 1)] : sink_;
                    const double tzm = (k > 0) ? in.data()[idx(i, j, k - 1)] : sink_;
                    const double kxp = kf(kn(i, j, k), kn(i + 1, j, k));
                    const double kxm = kf(kn(i, j, k), kn(i - 1, j, k));
                    const double kyp = kf(kn(i, j, k), kn(i, j + 1, k));
                    const double kym = kf(kn(i, j, k), kn(i, j - 1, k));
                    const double kzp = kf(kn(i, j, k), kn(i, j, k + 1));
                    const double kzm = kf(kn(i, j, k), kn(i, j, k - 1));
                    // B = −∇·(k∇): the SPD operator (positive diagonal) —
                    // the steady form written as B·T = q + boundary terms
                    const double self = (kxp + kxm) / (hx * hx) +
                                        (kyp + kym) / (hy * hy) +
                                        (kzp + kzm) / (hz * hz);
                    out.data()[id] = self * in.data()[id] -
                                     (kxp * txp + kxm * txm) / (hx * hx) -
                                     (kyp * typ + kym * tym) / (hy * hy) -
                                     (kzp * tzp + kzm * tzm) / (hz * hz);
                }
        (void)T;
        return true;
    }
    bool apply_transpose(const Field& in, Field& out,
                         exd::engine::core::ModelStatus& status) const override
    {
        return apply(in, out, status);   // symmetric stencil
    }
    bool diagonal(Field& out, exd::engine::core::ModelStatus& status) const override
    {
        if (out.size() != k_.size())
        {
            status.ok = false;
            status.error = "chiplet conduction: diagonal size mismatch";
            return false;
        }
        // not needed by CG here (no preconditioner); report the true diag
        apply_diag(out);
        return true;
    }
    bool jacobian_vector_product(const Field& x, const Field& v,
                                 Field& out, exd::engine::core::ModelStatus& status) const override
    {
        return apply(v, out, status);   // linear operator
    }

    const mesh::StructuredGrid& grid() const { return g_; }
    const std::vector<double>& k() const { return k_; }

private:
    void apply_diag(Field& out) const
    {
        const int32_t nx = g_.dims[0], ny = g_.dims[1], nz = g_.dims[2];
        auto idx = [&](int i, int j, int k) {
            return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                       (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
        };
        auto kf = [](double a, double b) { return (a + b > 0.0) ? 2.0 * a * b / (a + b) : 0.0; };
        auto kn = [&](int i, int j, int k) -> double {
            const int32_t im = (i < 0) ? 0 : (i >= nx ? nx - 1 : i);
            const int32_t jm = (j < 0) ? 0 : (j >= ny ? ny - 1 : j);
            const int32_t km = (k < 0) ? 0 : (k >= nz ? nz - 1 : k);
            return k_[idx(im, jm, km)];
        };
        const double hx = g_.spacing[0], hy = g_.spacing[1], hz = g_.spacing[2];
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const double kxp = kf(kn(i, j, k), kn(i + 1, j, k));
                    const double kxm = kf(kn(i, j, k), kn(i - 1, j, k));
                    const double kyp = kf(kn(i, j, k), kn(i, j + 1, k));
                    const double kym = kf(kn(i, j, k), kn(i, j - 1, k));
                    const double kzp = kf(kn(i, j, k), kn(i, j, k + 1));
                    const double kzm = kf(kn(i, j, k), kn(i, j, k - 1));
                    out.data()[idx(i, j, k)] = (kxp + kxm) / (hx * hx) +
                                               (kyp + kym) / (hy * hy) +
                                               (kzp + kzm) / (hz * hz);
                }
    }

    mesh::StructuredGrid g_;
    std::vector<double> k_;
    double sink_;
    exd::engine::core::OperatorInfo info_;
};

} // namespace

ChipletBoardResult solve_chiplet_board(const ChipletBoardConfig& config)
{
    ChipletBoardResult result;
    core::ModelStatus& status = result.status;
    if (!config.grid.validate(status))
    {
        result.ok = false;
        return result;
    }
    const size_t N = config.grid.node_count();
    const int32_t nx = config.grid.dims[0];
    const int32_t ny = config.grid.dims[1];
    const int32_t nz = config.grid.dims[2];
    const double cell_vol = config.grid.cell_volume();

    auto idx = [&](int i, int j, int k) {
        return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                   (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
    };

    // ── per-node conductivity (base + spreader regions) ──
    std::vector<double> cond_k(N, config.base_conductivity);
    for (const auto& sp : config.spreaders)
    {
        const int i0 = static_cast<int>(std::lround(sp.x / config.grid.spacing[0]) - sp.w_cells / 2);
        const int j0 = static_cast<int>(std::lround(sp.y / config.grid.spacing[1]) - sp.h_cells / 2);
        for (int j = std::max(0, j0); j < std::min<int32_t>(ny, j0 + sp.h_cells); ++j)
            for (int i = std::max(0, i0); i < std::min<int32_t>(nx, i0 + sp.w_cells); ++i)
                for (int k = 0; k < nz; ++k)
                    cond_k[idx(i, j, k)] = sp.conductivity;
    }

    // ── per-node source (chips over the board footprint; uniform over z) ──
    std::vector<double> source(N, 0.0);
    double total_power = 0.0;
    for (const auto& chip : config.chips)
    {
        const double dx = config.grid.spacing[0];
        const double dy = config.grid.spacing[1];
        const double footprint_x = static_cast<double>(chip.w_cells) * dx;
        const double footprint_y = static_cast<double>(chip.h_cells) * dy;
        const double vol = footprint_x * footprint_y * config.board_thickness;
        if (!(vol > 0.0))
        {
            status.ok = false;
            status.error = "chiplet: chip footprint must be at least one cell";
            result.ok = false;
            return result;
        }
        const double q = chip.power_watts / vol;   // W/m³ over its footprint
        total_power += chip.power_watts;
        const int i0 = static_cast<int>(std::lround(chip.x / dx) - chip.w_cells / 2);
        const int j0 = static_cast<int>(std::lround(chip.y / dy) - chip.h_cells / 2);
        for (int j = std::max(0, j0); j < std::min<int32_t>(ny, j0 + chip.h_cells); ++j)
            for (int i = std::max(0, i0); i < std::min<int32_t>(nx, i0 + chip.w_cells); ++i)
                for (int k = 0; k < nz; ++k)
                    source[idx(i, j, k)] += q;
    }
    result.total_power = total_power;

    // ── steady solve −∇·(k∇T) = q with the sink faces: the operator is
    //    affine (L·T = M·T + c(sink)); transfer the affine part to the RHS
    //    and solve the linear part with CG (porous-media pattern). ──
    Field rhs(FieldMetadata{
        .name = "rhs", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    for (size_t i = 0; i < N; ++i) rhs.data()[i] = source[i];
    // B·T = q − B(0): B is affine in the sink values; transfer the constant
    // to the RHS and solve the linear part with CG (porous-media pattern).
    ConductionOperator op(config.grid, cond_k, config.sink_temperature);
    Field zero(FieldMetadata{
        .name = "zero", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    zero.assign(0.0);
    Field b0(FieldMetadata{
        .name = "b0", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    op.apply(zero, b0, status);              // b0 = B(0) = the affine constant
    Field target(FieldMetadata{
        .name = "target", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    for (size_t i = 0; i < N; ++i) target.data()[i] = rhs.data()[i] - b0.data()[i];
    AffineLinearPart L(op, N);

    Field sol(FieldMetadata{
        .name = "T", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    sol.assign(config.sink_temperature);       // good initial guess

    IterativeSolverConfig cfg;
    cfg.tolerance = config.tolerance;
    cfg.max_iterations = config.max_iterations;
    auto rep = solve_cg(L, target, sol, cfg, status);
    if (!rep.converged && status.ok)
        status.warnings.push_back("chiplet: CG did not fully converge");
    if (!status.ok)
    {
        result.ok = false;
        return result;
    }
    result.max_residual = rep.final_residual;

    static_cast<mesh::StructuredGrid&>(result.temperature) = config.grid;
    result.temperature.values.resize(N);
    std::copy(sol.data().begin(), sol.data().end(), result.temperature.values.begin());

    // ── peak + sink flux (energy balance) ──
    double peak = -1e300;
    int pi = 0, pj = 0, pk = 0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const double T = result.temperature.values[idx(i, j, k)];
                if (T > peak) { peak = T; pi = i; pj = j; pk = k; }
            }
    result.peak_temperature = peak;
    result.peak_x = config.grid.origin[0] + pi * config.grid.spacing[0];
    result.peak_y = config.grid.origin[1] + pj * config.grid.spacing[1];

    // sink flux: the EXACT discrete energy balance — the sum of the
    // operator rows over the boundary nodes equals the net outward flux
    // (interior fluxes telescope pairwise), so
    //   sink_flux = cell_vol · Σ_{∂Ω} (B·T)_i
    Field bt(FieldMetadata{
        .name = "bt", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    op.apply(sol, bt, status);
    double flux = 0.0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const bool boundary = (i == 0 || i == nx - 1 ||
                                       j == 0 || j == ny - 1 ||
                                       k == 0 || k == nz - 1);
                if (boundary) flux += bt.data()[idx(i, j, k)];
            }
    flux *= cell_vol;
    result.sink_flux = flux;
    result.ok = true;
    return result;
}

} // namespace exd::engine::presets::electronics
