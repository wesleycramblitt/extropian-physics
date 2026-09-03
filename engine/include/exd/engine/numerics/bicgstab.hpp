#pragma once

// ─────────────────────────────────────────────────────
// BiCGSTAB (implementation_spec §36).
//
// For non-symmetric operators: uses apply() and
// apply_transpose() (fallback: apply when transpose is
// unsupported — BiCGSTAB without transpose reduces to a
// variant of CGS; prefer GMRES for strongly non-symmetric
// problems without a transpose).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <cmath>
#include <vector>

namespace exd::engine::numerics {

inline IterativeSolveReport solve_bicgstab(const LinearOperator& op,
                                           const core::Field& b,
                                           core::Field& x,
                                           const IterativeSolverConfig& cfg,
                                           core::ModelStatus& status)
{
    IterativeSolveReport rep;
    if (b.size() != x.size())
    {
        status.ok = false;
        status.error = "bicgstab: b and x size mismatch";
        return rep;
    }
    const size_t n = b.size();
    core::FieldMetadata meta{
        .name = "bicg_work", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node};
    core::Field r(meta, n), rhat(meta, n), p(meta, n), v(meta, n), s(meta, n), t(meta, n);
    core::Field tmp(meta, n);

    copy(b.data(), r.data());
    if (!op.apply(x, tmp, status)) return rep;
    axpy(-1.0, tmp.data(), r.data());       // r0 = b - A x0
    copy(r.data(), rhat.data());

    const double b_norm = norm2(b.data());
    if (b_norm == 0.0)
    {
        x.assign(0.0);
        rep.converged = true;
        return rep;
    }

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    copy(r.data(), p.data());
    double rho_prev = 1.0;

    const auto residual = [&]() { return norm2(r.data()) / b_norm; };
    rep.final_residual = residual();
    if (rep.final_residual < cfg.tolerance)
    {
        rep.converged = true;
        return rep;
    }

    for (int it = 1; it <= cfg.max_iterations; ++it)
    {
        const double rho_new = dot(rhat.data(), r.data());
        if (rho_new == 0.0)
        {
            status.ok = false;
            status.error = "bicgstab: breakdown (rho == 0)";
            return rep;
        }
        const double beta = (rho_new / rho_prev) * (alpha / omega);
        scale(beta, p.data());
        axpy(-beta * omega, v.data(), p.data());
        axpy(1.0, r.data(), p.data());

        if (!op.apply(p, tmp, status)) return rep;
        copy(tmp.data(), v.data());
        const double rhat_v = dot(rhat.data(), v.data());
        if (rhat_v == 0.0)
        {
            status.ok = false;
            status.error = "bicgstab: breakdown (rhat·v == 0)";
            return rep;
        }
        alpha = rho_new / rhat_v;
        copy(r.data(), s.data());
        axpy(-alpha, v.data(), s.data());
        if (norm2(s.data()) / b_norm < cfg.tolerance)
        {
            axpy(alpha, p.data(), x.data());
            rep.iterations = it;
            rep.final_residual = norm2(s.data()) / b_norm;
            rep.converged = true;
            return rep;
        }
        if (!op.apply(s, tmp, status)) return rep;
        copy(tmp.data(), t.data());
        const double tt = dot(t.data(), t.data());
        if (tt == 0.0)
        {
            status.ok = false;
            status.error = "bicgstab: breakdown (t == 0)";
            return rep;
        }
        omega = dot(t.data(), s.data()) / tt;
        axpy(alpha, p.data(), x.data());
        axpy(omega, s.data(), x.data());
        copy(s.data(), r.data());
        axpy(-omega, t.data(), r.data());
        rho_prev = rho_new;
        rep.iterations = it;
        rep.final_residual = residual();
        if (rep.final_residual < cfg.tolerance ||
            norm2(r.data()) < cfg.min_absolute_residual)
        {
            rep.converged = true;
            return rep;
        }
        if (omega == 0.0)
        {
            status.ok = false;
            status.error = "bicgstab: breakdown (omega == 0)";
            return rep;
        }
    }
    status.warnings.push_back("bicgstab: max iterations reached without convergence");
    return rep;
}

} // namespace exd::engine::numerics
