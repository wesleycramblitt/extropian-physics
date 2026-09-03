// Heterogeneous steady conduction — the general "heat conduction" class:
// per-node k and q, per-face BC kinds, affine Dirichlet operator + CG.

#include <exd/engine/physics/thermal/heterogeneous.hpp>

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/operator.hpp>
#include <exd/engine/numerics/cg.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <algorithm>
#include <cmath>

namespace exd::engine::physics::thermal {

namespace {

using exd::engine::core::Field;
using exd::engine::core::FieldMetadata;
using exd::engine::core::FieldRank;
using exd::engine::core::FieldLocation;
using exd::engine::numerics::AffineLinearPart;
using exd::engine::numerics::IterativeSolverConfig;
using exd::engine::numerics::solve_cg;

/// The heterogeneous conduction operator:
///   (B·T)_i = Σ_axes (kf⁺·(T − T⁺) + kf⁻·(T − T⁻))/h²
/// (the SPD form, diagonal-positive).  Out-of-range neighbors are the fixed
/// face value for FixedValue faces and the MIRROR (interior neighbor) for
/// Insulated faces — the same ghost semantics as the FDM operators.
class ConductionOperator final : public exd::engine::core::IOperator
{
public:
    ConductionOperator(mesh::StructuredGrid g, const std::vector<double>& k,
                       const HeterogeneousConductionConfig& cfg)
        : g_(std::move(g)), k_(k), cfg_(cfg)
    {
        info_.name = "heterogeneous_conduction";
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
            status.error = "conduction operator: field size mismatch";
            return false;
        }
        auto idx = [&](int i, int j, int k) {
            return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                       (static_cast<size_t>(j) + static_cast<size_t>(ny) * k);
        };
        auto kn = [&](int i, int j, int k) -> double {
            const int32_t im = (i < 0) ? 0 : (i >= nx ? nx - 1 : i);
            const int32_t jm = (j < 0) ? 0 : (j >= ny ? ny - 1 : j);
            const int32_t km = (k < 0) ? 0 : (k >= nz ? nz - 1 : k);
            return k_[idx(im, jm, km)];
        };
        auto kf = [](double a, double b) { return (a + b > 0.0) ? 2.0 * a * b / (a + b) : 0.0; };
        auto is_fixed = [&](int face) {
            return cfg_.face_kind[static_cast<size_t>(face)] == ConductionFaceKind::FixedValue;
        };
        auto face_val = [&](int face) { return cfg_.face_value[static_cast<size_t>(face)]; };
        const double hx = g_.spacing[0], hy = g_.spacing[1], hz = g_.spacing[2];

        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const size_t id = idx(i, j, k);
                    // ── FixedValue face node → pure pin row (κ·(T − v)):
                    //    the interior never couples to it (eliminated
                    //    Dirichlet), so the operator stays symmetric. ──
                    double pin_value = 0.0;
                    bool pinned = false;
                    if (i == 0 && is_fixed(0)) { pin_value = face_val(0); pinned = true; }
                    if (i == nx - 1 && is_fixed(1)) { pin_value = face_val(1); pinned = true; }
                    if (j == 0 && is_fixed(2)) { pin_value = face_val(2); pinned = true; }
                    if (j == ny - 1 && is_fixed(3)) { pin_value = face_val(3); pinned = true; }
                    if (k == 0 && is_fixed(4)) { pin_value = face_val(4); pinned = true; }
                    if (k == nz - 1 && is_fixed(5)) { pin_value = face_val(5); pinned = true; }
                    if (pinned)
                    {
                        const double kap = 2.0 / (hx * hx) + 2.0 / (hy * hy) + 2.0 / (hz * hz);
                        out.data()[id] = kap * (in.data()[id] - pin_value);
                        continue;
                    }

                    // ── general row: Σ over the included sides of
                    //    kf·(T_i − T_side); the side VALUE is
                    //    * the neighbor node (in range),
                    //    * the face VALUE as a CONSTANT when the neighbor is
                    //      a FixedValue-face node (eliminated DOF),
                    //    * the mirror when the side is an Insulated face.
                    //    Each included side contributes +kf to the diagonal,
                    //    so Insulated face rows are symmetric half-rows. ──
                    auto side = [&](int axis, int sgn, int ii, int jj, int kk) {
                        int i2 = ii, j2 = jj, k2 = kk;
                        bool in_range = true;
                        if (axis == 0) { i2 += sgn; in_range = (i2 >= 0 && i2 < nx); }
                        if (axis == 1) { j2 += sgn; in_range = (j2 >= 0 && j2 < ny); }
                        if (axis == 2) { k2 += sgn; in_range = (k2 >= 0 && k2 < nz); }
                        const int face = 2 * axis + (sgn > 0 ? 1 : 0);
                        if (!in_range)
                        {
                            // out of the domain: FixedValue → value;
                            // Insulated → mirror of the interior neighbor
                            if (is_fixed(face)) return std::pair<double, bool>{face_val(face), false};
                            return std::pair<double, bool>{in.data()[idx(ii - (axis == 0 ? sgn : 0),
                                                                          jj - (axis == 1 ? sgn : 0),
                                                                          kk - (axis == 2 ? sgn : 0))], true};
                        }
                        // neighbor is a FixedValue-face node → CONSTANT
                        bool nb_fixed = false;
                        if (i2 == 0 && is_fixed(0)) { nb_fixed = true; }
                        if (i2 == nx - 1 && is_fixed(1)) { nb_fixed = true; }
                        if (j2 == 0 && is_fixed(2)) { nb_fixed = true; }
                        if (j2 == ny - 1 && is_fixed(3)) { nb_fixed = true; }
                        if (k2 == 0 && is_fixed(4)) { nb_fixed = true; }
                        if (k2 == nz - 1 && is_fixed(5)) { nb_fixed = true; }
                        if (nb_fixed)
                        {
                            if (i2 == 0) return std::pair<double, bool>{face_val(0), false};
                            if (i2 == nx - 1) return std::pair<double, bool>{face_val(1), false};
                            if (j2 == 0) return std::pair<double, bool>{face_val(2), false};
                            if (j2 == ny - 1) return std::pair<double, bool>{face_val(3), false};
                            if (k2 == 0) return std::pair<double, bool>{face_val(4), false};
                            return std::pair<double, bool>{face_val(5), false};
                        }
                        return std::pair<double, bool>{in.data()[idx(i2, j2, k2)], true};
                    };

                    double out_val = 0.0;
                    // x sides
                    auto sxm = side(0, -1, i, j, k);
                    auto sxp = side(0, +1, i, j, k);
                    const double kxm = kf(kn(i, j, k), kn(i - 1, j, k));
                    const double kxp = kf(kn(i, j, k), kn(i + 1, j, k));
                    out_val += kxm * (in.data()[id] - sxm.first) / (hx * hx);
                    out_val += kxp * (in.data()[id] - sxp.first) / (hx * hx);
                    auto sym = side(1, -1, i, j, k);
                    auto syp = side(1, +1, i, j, k);
                    const double kym = kf(kn(i, j, k), kn(i, j - 1, k));
                    const double kyp = kf(kn(i, j, k), kn(i, j + 1, k));
                    out_val += kym * (in.data()[id] - sym.first) / (hy * hy);
                    out_val += kyp * (in.data()[id] - syp.first) / (hy * hy);
                    auto szm = side(2, -1, i, j, k);
                    auto szp = side(2, +1, i, j, k);
                    const double kzm = kf(kn(i, j, k), kn(i, j, k - 1));
                    const double kzp = kf(kn(i, j, k), kn(i, j, k + 1));
                    out_val += kzm * (in.data()[id] - szm.first) / (hz * hz);
                    out_val += kzp * (in.data()[id] - szp.first) / (hz * hz);
                    out.data()[id] = out_val;
                }
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
            status.error = "conduction operator: diagonal size mismatch";
            return false;
        }
        // The solve does not precondition; report the natural row diagonal
        // (the pin rows use the κ scale, half-rows for Insulated faces).
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
        auto is_fixed = [&](int face) {
            return cfg_.face_kind[static_cast<size_t>(face)] == ConductionFaceKind::FixedValue;
        };
        const double hx = g_.spacing[0], hy = g_.spacing[1], hz = g_.spacing[2];
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const size_t id = idx(i, j, k);
                    bool pinned = false;
                    if (i == 0 && is_fixed(0)) pinned = true;
                    if (i == nx - 1 && is_fixed(1)) pinned = true;
                    if (j == 0 && is_fixed(2)) pinned = true;
                    if (j == ny - 1 && is_fixed(3)) pinned = true;
                    if (k == 0 && is_fixed(4)) pinned = true;
                    if (k == nz - 1 && is_fixed(5)) pinned = true;
                    if (pinned)
                    {
                        out.data()[id] = 2.0 / (hx * hx) + 2.0 / (hy * hy) + 2.0 / (hz * hz);
                        continue;
                    }
                    // included sides: -x is excluded only when i==0 && x- is
                    // Insulated (mirror half-row); kf based on the clamped k
                    double self = 0.0;
                    if (i > 0 || !is_fixed(0)) self += kf(kn(i, j, k), kn(i - 1, j, k)) / (hx * hx);
                    if (i < nx - 1 || !is_fixed(1)) self += kf(kn(i, j, k), kn(i + 1, j, k)) / (hx * hx);
                    if (j > 0 || !is_fixed(2)) self += kf(kn(i, j, k), kn(i, j - 1, k)) / (hy * hy);
                    if (j < ny - 1 || !is_fixed(3)) self += kf(kn(i, j, k), kn(i, j + 1, k)) / (hy * hy);
                    if (k > 0 || !is_fixed(4)) self += kf(kn(i, j, k), kn(i, j, k - 1)) / (hz * hz);
                    if (k < nz - 1 || !is_fixed(5)) self += kf(kn(i, j, k), kn(i, j, k + 1)) / (hz * hz);
                    out.data()[id] = self;
                }
        return true;
    }
    bool jacobian_vector_product(const Field& x, const Field& v,
                                 Field& out, exd::engine::core::ModelStatus& status) const override
    {
        return apply(v, out, status);   // linear operator
    }

private:
    mesh::StructuredGrid g_;
    std::vector<double> k_;
    HeterogeneousConductionConfig cfg_;
    exd::engine::core::OperatorInfo info_;
};

/// Populate the per-node k and q from the region lists (box-in-box test,
/// deterministic; the LAST matching region wins).
void build_fields(const HeterogeneousConductionConfig& c, const size_t N,
                  std::vector<double>& k, std::vector<double>& q)
{
    k.assign(N, c.base_conductivity);
    q.assign(N, 0.0);
    const int32_t nx = c.grid.dims[0];
    const int32_t ny = c.grid.dims[1];
    const int32_t nz = c.grid.dims[2];
    auto idx = [&](int i, int j, int k3) {
        return static_cast<size_t>(i) + static_cast<size_t>(nx) *
                   (static_cast<size_t>(j) + static_cast<size_t>(ny) * k3);
    };
    auto in_box = [&](int i, int j, int k3, const std::array<double, 3>& center,
                      const std::array<double, 3>& half) -> bool {
        const std::array<int32_t, 3> index = {i, j, k3};
        for (int a = 0; a < 3; ++a)
        {
            const double pos = c.grid.origin[a] + static_cast<double>(index[a]) * c.grid.spacing[a];
            if (std::fabs(pos - center[a]) > half[a]) return false;
        }
        return true;
    };
    for (int k3 = 0; k3 < nz; ++k3)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                for (const auto& m : c.materials)
                    if (in_box(i, j, k3, m.center, m.half_extents))
                        k[idx(i, j, k3)] = m.conductivity;
                for (const auto& s : c.sources)
                    if (in_box(i, j, k3, s.center, s.half_extents))
                        q[idx(i, j, k3)] = s.volumetric_source;
            }
    // data-driven overrides replace the region-built fields entirely
    if (c.conductivity_field.size() == N) k = c.conductivity_field;
    if (c.source_field.size() == N) q = c.source_field;
}

} // namespace

HeterogeneousConductionResult solve_heterogeneous_conduction(
    const HeterogeneousConductionConfig& config)
{
    HeterogeneousConductionResult result;
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

    std::vector<double> k_field, q_field;
    build_fields(config, N, k_field, q_field);

    double total_power = 0.0;
    for (double qn : q_field) total_power += qn * cell_vol;
    result.total_power = total_power;

    // ── steady solve: B·T = q − B(0), affine transfer + linear part + CG ──
    ConductionOperator op(config.grid, k_field, config);
    Field zero(FieldMetadata{
        .name = "zero", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    zero.assign(0.0);
    Field b0(FieldMetadata{
        .name = "b0", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    op.apply(zero, b0, status);
    Field target(FieldMetadata{
        .name = "target", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    for (size_t i = 0; i < N; ++i) target.data()[i] = q_field[i] - b0.data()[i];
    AffineLinearPart L(op, N);

    Field sol(FieldMetadata{
        .name = "T", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    sol.assign(config.face_value[0]);       // good initial guess

    IterativeSolverConfig cfg;
    cfg.tolerance = config.tolerance;
    cfg.max_iterations = config.max_iterations;
    auto rep = solve_cg(L, target, sol, cfg, status);
    if (!rep.converged && status.ok)
        status.warnings.push_back("conduction: CG did not fully converge");
    if (!status.ok)
    {
        result.ok = false;
        return result;
    }
    result.max_residual = rep.final_residual;

    static_cast<mesh::StructuredGrid&>(result.temperature) = config.grid;
    result.temperature.values.resize(N);
    std::copy(sol.data().begin(), sol.data().end(), result.temperature.values.begin());

    // ── peak ──
    double peak = -1e300;
    int pi = 0, pj = 0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const double T = result.temperature.values[idx(i, j, k)];
                if (T > peak) { peak = T; pi = i; pj = j; }
            }
    result.peak_temperature = peak;
    result.peak_x = config.grid.origin[0] + pi * config.grid.spacing[0];
    result.peak_y = config.grid.origin[1] + pj * config.grid.spacing[1];

    // ── sink flux: the exact discrete energy balance — the FULL sum of the
    //    operator rows telescopes: interior face fluxes cancel pairwise and
    //    only the domain-boundary fluxes survive.  (Summing boundary rows
    //    alone is NOT the flux — the boundary nodes' interior pairings do
    //    not cancel within the boundary subset.)
    Field bt(FieldMetadata{
        .name = "bt", .rank = FieldRank::Scalar, .components = 1,
        .location = FieldLocation::Node}, N);
    op.apply(sol, bt, status);
    double flux = 0.0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                flux += bt.data()[idx(i, j, k)];
    result.sink_flux = flux * cell_vol;
    result.ok = true;
    return result;
}

} // namespace exd::engine::physics::thermal
