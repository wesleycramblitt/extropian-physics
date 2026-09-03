#pragma once

// ─────────────────────────────────────────────────────
// GMRES(m) (implementation_spec §36).
//
// Restarted GMRES for general non-symmetric operators.
// Requires apply() only (no transpose).  Arnoldi basis in
// std::vector storage; restart length from config.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <cmath>
#include <vector>

namespace exd::engine::numerics {

struct GmresConfig : IterativeSolverConfig
{
    int restart = 30;   // GMRES(m) restart length
};

inline IterativeSolveReport solve_gmres(const LinearOperator& op,
                                        const core::Field& b,
                                        core::Field& x,
                                        const GmresConfig& cfg,
                                        core::ModelStatus& status)
{
    IterativeSolveReport rep;
    if (b.size() != x.size())
    {
        status.ok = false;
        status.error = "gmres: b and x size mismatch";
        return rep;
    }
    const size_t n = b.size();
    const double b_norm = norm2(b.data());
    if (b_norm == 0.0)
    {
        x.assign(0.0);
        rep.converged = true;
        return rep;
    }

    core::Field tmp(core::FieldMetadata{
        .name = "tmp", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node}, n);
    std::vector<double> r(n);
    copy(b.data(), r);
    if (!op.apply(x, tmp, status)) return rep;
    axpy(-1.0, tmp.data(), r);

    const int m = cfg.restart > 0 ? cfg.restart : 30;
    std::vector<std::vector<double>> V(static_cast<size_t>(m + 1), std::vector<double>(n));
    std::vector<double> H((static_cast<size_t>(m) + 1) * static_cast<size_t>(m), 0.0);  // Hessenberg
    std::vector<double> cs(static_cast<size_t>(m)), sn(static_cast<size_t>(m));
    std::vector<double> g(static_cast<size_t>(m) + 1, 0.0);

    int total = 0;
    while (total < cfg.max_iterations)
    {
        const double beta = norm2(r);
        if (beta / b_norm < cfg.tolerance)
        {
            rep.converged = true;
            rep.final_residual = beta / b_norm;
            return rep;
        }
        copy(r, V[0]);
        scale(1.0 / beta, V[0]);
        g.assign(g.size(), 0.0);
        g[0] = beta;

        int k = 0;
        for (k = 0; k < m && total + k < cfg.max_iterations; ++k)
        {
            // w = A·V[k]
            core::Field wk(core::FieldMetadata{
                .name = "wk", .rank = core::FieldRank::Scalar, .components = 1,
                .location = core::FieldLocation::Node}, n);
            const core::Field vk(core::FieldMetadata{
                .name = "vk", .rank = core::FieldRank::Scalar, .components = 1,
                .location = core::FieldLocation::Node}, n);
            // copy V[k] into a Field for the operator
            core::Field vf = vk;
            copy(V[static_cast<size_t>(k)], vf.data());
            if (!op.apply(vf, wk, status)) return rep;

            // Arnoldi: orthogonalize against previous basis
            for (int j = 0; j <= k; ++j)
            {
                H[static_cast<size_t>(j) * static_cast<size_t>(m) + static_cast<size_t>(k)] =
                    dot(V[static_cast<size_t>(j)], wk.data());
                axpy(-H[static_cast<size_t>(j) * static_cast<size_t>(m) + static_cast<size_t>(k)],
                     V[static_cast<size_t>(j)], wk.data());
            }
            const double hnext = norm2(wk.data());
            H[static_cast<size_t>(k + 1) * static_cast<size_t>(m) + static_cast<size_t>(k)] = hnext;
            if (hnext > 0.0)
            {
                scale(1.0 / hnext, wk.data());
                copy(wk.data(), V[static_cast<size_t>(k + 1)]);
            }

            // Givens rotations on H column k
            for (int j = 0; j < k; ++j)
            {
                const double h1 = H[static_cast<size_t>(j) * static_cast<size_t>(m) + static_cast<size_t>(k)];
                const double h2 = H[static_cast<size_t>(j + 1) * static_cast<size_t>(m) + static_cast<size_t>(k)];
                H[static_cast<size_t>(j) * static_cast<size_t>(m) + static_cast<size_t>(k)] =
                    cs[static_cast<size_t>(j)] * h1 + sn[static_cast<size_t>(j)] * h2;
                H[static_cast<size_t>(j + 1) * static_cast<size_t>(m) + static_cast<size_t>(k)] =
                    -sn[static_cast<size_t>(j)] * h1 + cs[static_cast<size_t>(j)] * h2;
            }
            const double h_kk = H[static_cast<size_t>(k) * static_cast<size_t>(m) + static_cast<size_t>(k)];
            const double h_k1k = H[static_cast<size_t>(k + 1) * static_cast<size_t>(m) + static_cast<size_t>(k)];
            const double denom = std::sqrt(h_kk * h_kk + h_k1k * h_k1k);
            if (denom == 0.0)
            {
                status.ok = false;
                status.error = "gmres: breakdown (singular Hessenberg)";
                return rep;
            }
            cs[static_cast<size_t>(k)] = h_kk / denom;
            sn[static_cast<size_t>(k)] = h_k1k / denom;
            H[static_cast<size_t>(k) * static_cast<size_t>(m) + static_cast<size_t>(k)] = denom;
            H[static_cast<size_t>(k + 1) * static_cast<size_t>(m) + static_cast<size_t>(k)] = 0.0;

            g[static_cast<size_t>(k + 1)] = -sn[static_cast<size_t>(k)] * g[static_cast<size_t>(k)];
            g[static_cast<size_t>(k)] = cs[static_cast<size_t>(k)] * g[static_cast<size_t>(k)];

            rep.final_residual = std::abs(g[static_cast<size_t>(k + 1)]) / b_norm;
            rep.iterations = total + k + 1;
            if (rep.final_residual < cfg.tolerance || hnext < 1e-300)
            {
                ++k;
                break;
            }
        }

        // Back-substitute upper Hessenberg (k x k)
        std::vector<double> y(static_cast<size_t>(k), 0.0);
        for (int jj = k - 1; jj >= 0; --jj)
        {
            double acc = g[static_cast<size_t>(jj)];
            for (int j2 = jj + 1; j2 < k; ++j2)
                acc -= H[static_cast<size_t>(jj) * static_cast<size_t>(m) + static_cast<size_t>(j2)] *
                       y[static_cast<size_t>(j2)];
            y[static_cast<size_t>(jj)] = acc / H[static_cast<size_t>(jj) * static_cast<size_t>(m) + static_cast<size_t>(jj)];
        }
        for (int jj = 0; jj < k; ++jj)
            axpy(y[static_cast<size_t>(jj)], V[static_cast<size_t>(jj)], x.data());

        // Residual for the restart loop
        copy(b.data(), r);
        if (!op.apply(x, tmp, status)) return rep;
        axpy(-1.0, tmp.data(), r);
        total += k;
        if (rep.final_residual < cfg.tolerance)
        {
            rep.converged = true;
            return rep;
        }
    }
    status.warnings.push_back("gmres: max iterations reached without convergence");
    return rep;
}

} // namespace exd::engine::numerics
