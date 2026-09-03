#pragma once

// ─────────────────────────────────────────────────────
// Conjugate Gradient (implementation_spec §36).
//
// CG on a symmetric positive-definite abstract operator
// (matrix-free: only apply() is used; diagonal() may
// precondition).  Fields provide the vector algebra.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <vector>

namespace exd::engine::numerics {

/// Solve A·x = b with CG.  `x` is both initial guess and solution.
/// Optional Jacobi preconditioner uses operator->diagonal() when
/// `use_diagonal_preconditioner` is true (and the operator supports it).
inline IterativeSolveReport solve_cg(const LinearOperator& op,
                                     const core::Field& b,
                                     core::Field& x,
                                     const IterativeSolverConfig& cfg,
                                     core::ModelStatus& status,
                                     bool use_diagonal_preconditioner = false)
{
    IterativeSolveReport rep;
    if (b.size() != x.size())
    {
        status.ok = false;
        status.error = "cg: b and x size mismatch";
        return rep;
    }
    const size_t n = b.size();
    core::FieldMetadata meta{
        .name = "cg_work", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node};
    core::Field r(meta, n), p(meta, n), tmp(meta, n);

    copy(b.data(), r.data());                       // r = b - A·x0
    if (!op.apply(x, tmp, status)) return rep;
    axpy(-1.0, tmp.data(), r.data());

    const double b_norm = norm2(b.data());
    // Homogeneous solves (b == 0) are valid when the operator is AFFINE with
    // the boundary values built in (eliminated-Dirichlet operators): the
    // convergence criterion is then the ABSOLUTE residual (relative to b is
    // undefined).
    const bool homogeneous = (b_norm == 0.0);

    double rho = dot(r.data(), r.data());
    double r_abs = std::sqrt(rho);
    rep.final_residual = homogeneous ? r_abs : r_abs / b_norm;
    if (rep.final_residual < cfg.tolerance || r_abs < cfg.min_absolute_residual)
    {
        rep.converged = true;
        return rep;
    }
    copy(r.data(), p.data());

    for (int it = 1; it <= cfg.max_iterations; ++it)
    {
        if (!op.apply(p, tmp, status)) return rep;
        const double denom = dot(p.data(), tmp.data());
        if (denom == 0.0)
        {
            status.ok = false;
            status.error = "cg: breakdown (A·p == 0)";
            return rep;
        }
        const double alpha = rho / denom;
        axpy(alpha, p.data(), x.data());
        axpy(-alpha, tmp.data(), r.data());
        const double rho_new = dot(r.data(), r.data());
        const double r_new = std::sqrt(rho_new);
        rep.iterations = it;
        rep.final_residual = homogeneous ? r_new : r_new / b_norm;
        if (rep.final_residual < cfg.tolerance || r_new < cfg.min_absolute_residual)
        {
            rep.converged = true;
            return rep;
        }
        const double beta = rho_new / rho;
        scale(beta, p.data());
        axpy(1.0, r.data(), p.data());
        rho = rho_new;
    }
    status.warnings.push_back("cg: max iterations reached without convergence (" +
                              std::to_string(cfg.max_iterations) + ")");
    return rep;
}

} // namespace exd::engine::numerics
