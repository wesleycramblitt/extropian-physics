#pragma once

// ─────────────────────────────────────────────────────
// FDM operators (implementation_spec §26, §34–§35).
//
// Reusable stencil operators on structured grids:
//   gradient, divergence, curl, Laplacian, upwind advection
// with ghost-cell boundary stencils (mirror/Neumann or
// Dirichlet per face) and contiguous-memory indexing
// (index i + nx·(j + ny·k)).  Operators are physics-free
// and backend-agnostic: they run through the CPU backend
// primitives today and can be re-targeted to CUDA without
// touching physics modules (spec §3.3, §40).
//
// `FdmLaplacianOperator` is a core::IOperator (matrix-free:
// apply/apply_transpose/diagonal/jvp) so iterative linear
// solvers (numerics/cg.hpp etc.) consume it directly —
// no global matrix assembly (§35).
// ─────────────────────────────────────────────────────

#include <exd/engine/backends/cpu.hpp>
#include <exd/engine/core/field.hpp>
#include <exd/engine/core/operator.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace exd::engine::discretization::fdm {

// Structured lattice types are owned by the mesh layer (spec §9).
using mesh::StructuredGrid;

/// Ghost-cell boundary treatment per structured-box face.
struct FaceGhostSpec
{
    bool dirichlet = false;   // false → zero-gradient mirror
    double value = 0.0;       // Dirichlet value (or 0 for mirror)
};

struct GhostSpec
{
    std::array<FaceGhostSpec, 6> faces{};  // order: XNeg, XPos, YNeg, YPos, ZNeg, ZPos
};

// ── indexing helpers (contiguous, spec §3.1) ─────────
inline size_t idx(const StructuredGrid& g, int32_t i, int32_t j, int32_t k)
{
    const size_t nx = static_cast<size_t>(g.dims[0]);
    const size_t ny = static_cast<size_t>(g.dims[1]);
    return static_cast<size_t>(i) + nx * (static_cast<size_t>(j) + ny * static_cast<size_t>(k));
}

/// Gradient of a scalar node field (central differences, O(h²)).
inline bool gradient(const StructuredGrid& g, std::span<const double> u,
                     std::span<double> gx, std::span<double> gy, std::span<double> gz,
                     const GhostSpec& ghosts, ModelStatus& status)
{
    if (u.size() != g.node_count() || gx.size() < g.node_count() ||
        gy.size() < g.node_count() || gz.size() < g.node_count())
    {
        status.ok = false;
        status.error = "fdm gradient: field size mismatch";
        return false;
    }
    const auto& b = backends::cpu_backend();
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    b.parallel_for(static_cast<size_t>(nx * ny * nz), [&](size_t n) {
        const int32_t k = static_cast<int32_t>(n / (nx * ny));
        const int32_t j = static_cast<int32_t>((n / nx) % ny);
        const int32_t i = static_cast<int32_t>(n % nx);
        const double hx = g.spacing[0], hy = g.spacing[1], hz = g.spacing[2];
        auto at = [&](int32_t ii, int32_t jj, int32_t kk) -> double {
            const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
            const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
            const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
            return u[idx(g, im, jm, km)];
        };
        auto ghost = [&](int32_t face) -> const FaceGhostSpec& { return ghosts.faces[static_cast<size_t>(face)]; };
        // per-axis: interior central; boundary reads through ghost value
        const double u_ip = (i < nx - 1) ? at(i + 1, j, k) : (ghost(1).dirichlet ? 2 * ghost(1).value - at(i - 1, j, k) : at(i - 1, j, k));
        const double u_im = (i > 0) ? at(i - 1, j, k) : (ghost(0).dirichlet ? 2 * ghost(0).value - at(i + 1, j, k) : at(i + 1, j, k));
        const double u_jp = (j < ny - 1) ? at(i, j + 1, k) : (ghost(3).dirichlet ? 2 * ghost(3).value - at(i, j - 1, k) : at(i, j - 1, k));
        const double u_jm = (j > 0) ? at(i, j - 1, k) : (ghost(2).dirichlet ? 2 * ghost(2).value - at(i, j + 1, k) : at(i, j + 1, k));
        const double u_kp = (k < nz - 1) ? at(i, j, k + 1) : (ghost(5).dirichlet ? 2 * ghost(5).value - at(i, j, k - 1) : at(i, j, k - 1));
        const double u_km = (k > 0) ? at(i, j, k - 1) : (ghost(4).dirichlet ? 2 * ghost(4).value - at(i, j, k + 1) : at(i, j, k + 1));
        gx[n] = (u_ip - u_im) / (2 * hx);
        gy[n] = (u_jp - u_jm) / (2 * hy);
        gz[n] = (u_kp - u_km) / (2 * hz);
    });
    return true;
}

/// Divergence of a vector node field (central differences, O(h²)).
inline bool divergence(const StructuredGrid& g, std::span<const double> vx,
                       std::span<const double> vy, std::span<const double> vz,
                       std::span<double> out, const GhostSpec& ghosts,
                       ModelStatus& status)
{
    if (vx.size() != g.node_count() || vy.size() != g.node_count() ||
        vz.size() != g.node_count() || out.size() < g.node_count())
    {
        status.ok = false;
        status.error = "fdm divergence: field size mismatch";
        return false;
    }
    const auto& b = backends::cpu_backend();
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    b.parallel_for(static_cast<size_t>(nx * ny * nz), [&](size_t n) {
        const int32_t k = static_cast<int32_t>(n / (nx * ny));
        const int32_t j = static_cast<int32_t>((n / nx) % ny);
        const int32_t i = static_cast<int32_t>(n % nx);
        const double hx = g.spacing[0], hy = g.spacing[1], hz = g.spacing[2];
        auto at = [&](std::span<const double> f, int32_t ii, int32_t jj, int32_t kk) -> double {
            const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
            const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
            const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
            return f[idx(g, im, jm, km)];
        };
        auto ghost = [&](int32_t face) -> const FaceGhostSpec& { return ghosts.faces[static_cast<size_t>(face)]; };
        const double xp = (i < nx - 1) ? at(vx, i + 1, j, k) : (ghost(1).dirichlet ? 2 * ghost(1).value - at(vx, i - 1, j, k) : at(vx, i - 1, j, k));
        const double xm = (i > 0) ? at(vx, i - 1, j, k) : (ghost(0).dirichlet ? 2 * ghost(0).value - at(vx, i + 1, j, k) : at(vx, i + 1, j, k));
        const double yp = (j < ny - 1) ? at(vy, i, j + 1, k) : (ghost(3).dirichlet ? 2 * ghost(3).value - at(vy, i, j - 1, k) : at(vy, i, j - 1, k));
        const double ym = (j > 0) ? at(vy, i, j - 1, k) : (ghost(2).dirichlet ? 2 * ghost(2).value - at(vy, i, j + 1, k) : at(vy, i, j + 1, k));
        const double zp = (k < nz - 1) ? at(vz, i, j, k + 1) : (ghost(5).dirichlet ? 2 * ghost(5).value - at(vz, i, j, k - 1) : at(vz, i, j, k - 1));
        const double zm = (k > 0) ? at(vz, i, j, k - 1) : (ghost(4).dirichlet ? 2 * ghost(4).value - at(vz, i, j, k + 1) : at(vz, i, j, k + 1));
        out[n] = (xp - xm) / (2 * hx) + (yp - ym) / (2 * hy) + (zp - zm) / (2 * hz);
    });
    return true;
}

/// Curl of a vector node field (central differences, O(h²)).
inline bool curl(const StructuredGrid& g, std::span<const double> vx,
                 std::span<const double> vy, std::span<const double> vz,
                 std::span<double> cx, std::span<double> cy, std::span<double> cz,
                 const GhostSpec& ghosts, ModelStatus& status)
{
    if (vx.size() != g.node_count() || vy.size() != g.node_count() ||
        vz.size() != g.node_count() || cx.size() < g.node_count() ||
        cy.size() < g.node_count() || cz.size() < g.node_count())
    {
        status.ok = false;
        status.error = "fdm curl: field size mismatch";
        return false;
    }
    const auto& b = backends::cpu_backend();
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    b.parallel_for(static_cast<size_t>(nx * ny * nz), [&](size_t n) {
        const int32_t k = static_cast<int32_t>(n / (nx * ny));
        const int32_t j = static_cast<int32_t>((n / nx) % ny);
        const int32_t i = static_cast<int32_t>(n % nx);
        const double hx = g.spacing[0], hy = g.spacing[1], hz = g.spacing[2];
        auto at = [&](std::span<const double> f, int32_t ii, int32_t jj, int32_t kk) -> double {
            const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
            const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
            const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
            return f[idx(g, im, jm, km)];
        };
        auto ghost = [&](int32_t face) -> const FaceGhostSpec& { return ghosts.faces[static_cast<size_t>(face)]; };
        auto deriv = [&](std::span<const double> f, int dim, int32_t ii, int32_t jj, int32_t kk) -> double {
            const double h = g.spacing[static_cast<size_t>(dim)];
            const int32_t nx0 = g.dims[0], ny0 = g.dims[1], nz0 = g.dims[2];
            if (dim == 0)
            {
                const double fp = (ii < nx0 - 1) ? at(f, ii + 1, jj, kk)
                                 : (ghost(1).dirichlet ? 2 * ghost(1).value - at(f, ii - 1, jj, kk) : at(f, ii - 1, jj, kk));
                const double fm = (ii > 0) ? at(f, ii - 1, jj, kk)
                                 : (ghost(0).dirichlet ? 2 * ghost(0).value - at(f, ii + 1, jj, kk) : at(f, ii + 1, jj, kk));
                return (fp - fm) / (2 * h);
            }
            if (dim == 1)
            {
                const double fp = (jj < ny0 - 1) ? at(f, ii, jj + 1, kk)
                                 : (ghost(3).dirichlet ? 2 * ghost(3).value - at(f, ii, jj - 1, kk) : at(f, ii, jj - 1, kk));
                const double fm = (jj > 0) ? at(f, ii, jj - 1, kk)
                                 : (ghost(2).dirichlet ? 2 * ghost(2).value - at(f, ii, jj + 1, kk) : at(f, ii, jj + 1, kk));
                return (fp - fm) / (2 * h);
            }
            const double fp = (kk < nz0 - 1) ? at(f, ii, jj, kk + 1)
                             : (ghost(5).dirichlet ? 2 * ghost(5).value - at(f, ii, jj, kk - 1) : at(f, ii, jj, kk - 1));
            const double fm = (kk > 0) ? at(f, ii, jj, kk - 1)
                             : (ghost(4).dirichlet ? 2 * ghost(4).value - at(f, ii, jj, kk + 1) : at(f, ii, jj, kk + 1));
            return (fp - fm) / (2 * h);
        };
        cz[n] = deriv(vy, 0, i, j, k) - deriv(vx, 1, i, j, k);
        cx[n] = deriv(vz, 1, i, j, k) - deriv(vy, 2, i, j, k);
        cy[n] = deriv(vx, 2, i, j, k) - deriv(vz, 0, i, j, k);
    });
    return true;
}

/// 7-point Laplacian (central, O(h²)) with ghost-cell boundaries.
inline bool laplacian(const StructuredGrid& g, std::span<const double> u,
                      std::span<double> out, const GhostSpec& ghosts,
                      ModelStatus& status)
{
    if (u.size() != g.node_count() || out.size() < g.node_count())
    {
        status.ok = false;
        status.error = "fdm laplacian: field size mismatch";
        return false;
    }
    const auto& b = backends::cpu_backend();
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    b.parallel_for(static_cast<size_t>(nx * ny * nz), [&](size_t n) {
        const int32_t k = static_cast<int32_t>(n / (nx * ny));
        const int32_t j = static_cast<int32_t>((n / nx) % ny);
        const int32_t i = static_cast<int32_t>(n % nx);
        const double hx2 = g.spacing[0] * g.spacing[0];
        const double hy2 = g.spacing[1] * g.spacing[1];
        const double hz2 = g.spacing[2] * g.spacing[2];
        auto at = [&](int32_t ii, int32_t jj, int32_t kk) -> double {
            const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
            const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
            const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
            return u[idx(g, im, jm, km)];
        };
        auto ghost = [&](int32_t face) -> const FaceGhostSpec& { return ghosts.faces[static_cast<size_t>(face)]; };
        // Symmetric boundary stencils (spec §26 "boundary stencils"):
        //   Dirichlet face  — face node row is pure diagonal; the adjacent
        //     interior row treats the face value as a CONSTANT (the DOF is
        //     eliminated) so the operator stays symmetric.
        //   mirror face     — one-sided half-row (u_adj − u_face)/h², which
        //     matches the interior row's +1/h² coupling (symmetric Neumann).
        // Interior rows keep the central 7-point stencil (exact for
        // quadratics); face rows are first-order.
        auto second_x = [&](int32_t ii) {
            if (ii > 0 && ii < nx - 1)
                return (at(ii + 1, j, k) - 2 * at(ii, j, k) + at(ii - 1, j, k)) / hx2;
            if (ii == 0)
            {
                if (ghost(0).dirichlet) return (2 * ghost(0).value - 2 * at(0, j, k)) / hx2;
                return (at(1, j, k) - at(0, j, k)) / hx2;      // mirror half-row
            }
            // ii == nx - 1
            if (ghost(1).dirichlet) return (2 * ghost(1).value - 2 * at(nx - 1, j, k)) / hx2;
            return (at(nx - 2, j, k) - at(nx - 1, j, k)) / hx2; // mirror half-row
        };
        auto second_y = [&](int32_t jj) {
            if (jj > 0 && jj < ny - 1)
                return (at(i, jj + 1, k) - 2 * at(i, jj, k) + at(i, jj - 1, k)) / hy2;
            if (jj == 0)
            {
                if (ghost(2).dirichlet) return (2 * ghost(2).value - 2 * at(i, 0, k)) / hy2;
                return (at(i, 1, k) - at(i, 0, k)) / hy2;
            }
            if (ghost(3).dirichlet) return (2 * ghost(3).value - 2 * at(i, ny - 1, k)) / hy2;
            return (at(i, ny - 2, k) - at(i, ny - 1, k)) / hy2;
        };
        auto second_z = [&](int32_t kk) {
            if (kk > 0 && kk < nz - 1)
                return (at(i, j, kk + 1) - 2 * at(i, j, kk) + at(i, j, kk - 1)) / hz2;
            if (kk == 0)
            {
                if (ghost(4).dirichlet) return (2 * ghost(4).value - 2 * at(i, j, 0)) / hz2;
                return (at(i, j, 1) - at(i, j, 0)) / hz2;
            }
            if (ghost(5).dirichlet) return (2 * ghost(5).value - 2 * at(i, j, nz - 1)) / hz2;
            return (at(i, j, nz - 2) - at(i, j, nz - 1)) / hz2;
        };
        // Interior rows adjacent to a Dirichlet face: neighbor value is the
        // FIXED boundary value (DOF eliminated), keeping the operator
        // symmetric (§35 matrix-free operator).
        double xpart = second_x(i);
        if (i == 1 && ghost(0).dirichlet)
            xpart = (ghost(0).value - 2 * at(i, j, k) + at(i + 1, j, k)) / hx2;
        else if (i == nx - 2 && ghost(1).dirichlet)
            xpart = (ghost(1).value - 2 * at(i, j, k) + at(i - 1, j, k)) / hx2;
        double ypart = second_y(j);
        if (j == 1 && ghost(2).dirichlet)
            ypart = (ghost(2).value - 2 * at(i, j, k) + at(i, j + 1, k)) / hy2;
        else if (j == ny - 2 && ghost(3).dirichlet)
            ypart = (ghost(3).value - 2 * at(i, j, k) + at(i, j - 1, k)) / hy2;
        double zpart = second_z(k);
        if (k == 1 && ghost(4).dirichlet)
            zpart = (ghost(4).value - 2 * at(i, j, k) + at(i, j, k + 1)) / hz2;
        else if (k == nz - 2 && ghost(5).dirichlet)
            zpart = (ghost(5).value - 2 * at(i, j, k) + at(i, j, k - 1)) / hz2;
        out[n] = xpart + ypart + zpart;
    });
    return true;
}

/// First-order upwind advection of scalar φ by velocity field (u·∇φ).
/// Returns the advection term (u·∇φ) at every node.
inline bool upwind_advect(const StructuredGrid& g, std::span<const double> vx,
                          std::span<const double> vy, std::span<const double> vz,
                          std::span<const double> phi, std::span<double> out,
                          const GhostSpec& ghosts, ModelStatus& status)
{
    if (vx.size() != g.node_count() || phi.size() != g.node_count() ||
        out.size() < g.node_count())
    {
        status.ok = false;
        status.error = "fdm upwind_advect: field size mismatch";
        return false;
    }
    const auto& b = backends::cpu_backend();
    const int32_t nx = g.dims[0], ny = g.dims[1], nz = g.dims[2];
    b.parallel_for(static_cast<size_t>(nx * ny * nz), [&](size_t n) {
        const int32_t k = static_cast<int32_t>(n / (nx * ny));
        const int32_t j = static_cast<int32_t>((n / nx) % ny);
        const int32_t i = static_cast<int32_t>(n % nx);
        auto at = [&](int32_t ii, int32_t jj, int32_t kk) -> double {
            const int32_t im = (ii < 0) ? 0 : (ii >= nx ? nx - 1 : ii);
            const int32_t jm = (jj < 0) ? 0 : (jj >= ny ? ny - 1 : jj);
            const int32_t km = (kk < 0) ? 0 : (kk >= nz ? nz - 1 : kk);
            return phi[idx(g, im, jm, km)];
        };
        auto ghost = [&](int32_t face) -> const FaceGhostSpec& { return ghosts.faces[static_cast<size_t>(face)]; };
        auto upwind = [&](double v, double f_p, double f_m, double f_c, double h) -> double {
            const double f_up = (v >= 0) ? f_m : f_p;  // upwind side relative to node
            return v * (f_c - f_up) / h;
        };
        const double phi_ip = (i < nx - 1) ? at(i + 1, j, k) : (ghost(1).dirichlet ? 2 * ghost(1).value - at(i - 1, j, k) : at(i - 1, j, k));
        const double phi_im = (i > 0) ? at(i - 1, j, k) : (ghost(0).dirichlet ? 2 * ghost(0).value - at(i + 1, j, k) : at(i + 1, j, k));
        const double phi_jp = (j < ny - 1) ? at(i, j + 1, k) : (ghost(3).dirichlet ? 2 * ghost(3).value - at(i, j - 1, k) : at(i, j - 1, k));
        const double phi_jm = (j > 0) ? at(i, j - 1, k) : (ghost(2).dirichlet ? 2 * ghost(2).value - at(i, j + 1, k) : at(i, j + 1, k));
        const double phi_kp = (k < nz - 1) ? at(i, j, k + 1) : (ghost(5).dirichlet ? 2 * ghost(5).value - at(i, j, k - 1) : at(i, j, k - 1));
        const double phi_km = (k > 0) ? at(i, j, k - 1) : (ghost(4).dirichlet ? 2 * ghost(4).value - at(i, j, k + 1) : at(i, j, k + 1));
        out[n] = upwind(vx[n], phi_ip, phi_im, phi[n], g.spacing[0]) +
                 upwind(vy[n], phi_jp, phi_jm, phi[n], g.spacing[1]) +
                 upwind(vz[n], phi_kp, phi_km, phi[n], g.spacing[2]);
    });
    return true;
}

// ─────────────────────────────────────────────────────
// Matrix-free Laplacian operator (§35) — the operator
// iterative linear solvers consume (numerics/cg.hpp…).
// ─────────────────────────────────────────────────────

class FdmLaplacianOperator final : public core::IOperator
{
public:
    FdmLaplacianOperator(mesh::StructuredGrid grid, GhostSpec ghosts)
        : grid_(std::move(grid)), ghosts_(std::move(ghosts))
    {
        info_.name = "fdm_laplacian";
        info_.inputs.push_back(core::FieldMetadata{
            .name = "u", .rank = core::FieldRank::Scalar, .components = 1,
            .location = core::FieldLocation::Node,
        });
        info_.outputs.push_back(core::FieldMetadata{
            .name = "lap_u", .rank = core::FieldRank::Scalar, .components = 1,
            .location = core::FieldLocation::Node,
        });
    }

    const core::OperatorInfo& info() const override { return info_; }

    bool apply(const core::Field& in, core::Field& out, core::ModelStatus& status) const override
    {
        if (in.size() != out.size() || in.size() != grid_.node_count())
        {
            status.ok = false;
            status.error = "fdm laplacian operator: field size mismatch";
            return false;
        }
        return laplacian(grid_, in.data(), out.data(), ghosts_, status);
    }

    // Symmetric operator: apply_transpose == apply.
    bool apply_transpose(const core::Field& in, core::Field& out,
                         core::ModelStatus& status) const override
    {
        return apply(in, out, status);
    }

    bool diagonal(core::Field& out, core::ModelStatus& status) const override
    {
        if (out.size() != grid_.node_count())
        {
            status.ok = false;
            status.error = "fdm laplacian operator: diagonal size mismatch";
            return false;
        }
        const double d = -(2.0 / (grid_.spacing[0] * grid_.spacing[0]) +
                          2.0 / (grid_.spacing[1] * grid_.spacing[1]) +
                          2.0 / (grid_.spacing[2] * grid_.spacing[2]));
        for (auto& v : out.data()) v = d;  // 7-point stencil: -2/h² per axis
        return true;
    }

    bool jacobian_vector_product(const core::Field& x, const core::Field& v,
                                 core::Field& out, core::ModelStatus& status) const override
    {
        // Linear operator: J(x)·v = A·v
        return apply(v, out, status);
    }

    const mesh::StructuredGrid& grid() const { return grid_; }
    const GhostSpec& ghosts() const { return ghosts_; }

private:
    mesh::StructuredGrid grid_;
    GhostSpec ghosts_;
    core::OperatorInfo info_;
};

} // namespace exd::engine::discretization::fdm
