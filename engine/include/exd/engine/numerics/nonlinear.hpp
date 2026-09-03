#pragma once

// ─────────────────────────────────────────────────────
// Nonlinear solvers (implementation_spec §37).
//
// Fixed-point iteration and Newton's method on the
// matrix-free architecture: a residual function
//   R(x) -> r
// plus either an operator-provided JVP or a finite-
// difference Jacobian-vector product.  Newton solves the
// linear system J·δ = -r with an inner iterative solver
// (default: GMRES — Newton-Krylov).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/numerics/gmres.hpp>
#include <exd/engine/numerics/linear_operator.hpp>

#include <cmath>
#include <functional>
#include <vector>

namespace exd::engine::numerics {

using ResidualFn = std::function<bool(const core::Field& x, core::Field& r, core::ModelStatus&)>;

/// Jacobian-vector product J(x)·v.  Default implementation uses forward
/// finite differences of the residual:
///   J(x)·v ≈ (R(x + εv) - R(x)) / ε
/// Override with an analytic operator when available.
using JvpFn = std::function<bool(const core::Field& x, const core::Field& v,
                                 core::Field& out, core::ModelStatus&)>;

struct NonlinearSolverConfig
{
    double tolerance = 1e-9;       // relative residual tolerance
    int max_iterations = 50;
    double relaxation = 1.0;       // under-relaxation for fixed-point/Newton step
    double fd_epsilon = 1e-7;      // finite-difference step for JVP fallback
};

struct NonlinearSolveReport
{
    bool converged = false;
    int iterations = 0;
    double final_residual = 0.0;
    double residual_norm_0 = 0.0;
};

inline JvpFn make_finite_difference_jvp(const ResidualFn& residual, double epsilon)
{
    return [residual, epsilon](const core::Field& x, const core::Field& v,
                               core::Field& out, core::ModelStatus& status) -> bool {
        core::Field xp = x;
        core::Field rp = out;
        core::Field r0 = out;
        // scale epsilon by ||v|| to keep the perturbation meaningful
        double vn = norm2(v.data());
        if (vn == 0.0)
        {
            out.assign(0.0);
            return true;
        }
        const double h = epsilon / vn;
        for (size_t i = 0; i < xp.size(); ++i) xp.data()[i] = x.at(i) + h * v.at(i);
        if (!residual(x, r0, status)) return false;
        if (!residual(xp, rp, status)) return false;
        for (size_t i = 0; i < out.size(); ++i)
            out.data()[i] = (rp.at(i) - r0.at(i)) / h;
        return true;
    };
}

/// Newton-Krylov: solve R(x) = 0 with GMRES on the JVP operator.
inline NonlinearSolveReport solve_newton(const ResidualFn& residual,
                                         core::Field& x,
                                         const NonlinearSolverConfig& cfg,
                                         core::ModelStatus& status,
                                         const JvpFn& jvp = JvpFn{})
{
    NonlinearSolveReport rep;
    const size_t n = x.size();
    core::Field r(core::FieldMetadata{
        .name = "r", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node}, n);
    core::Field d(core::FieldMetadata{
        .name = "d", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node}, n);

    if (!residual(x, r, status)) return rep;
    rep.residual_norm_0 = norm2(r.data());
    if (rep.residual_norm_0 < cfg.tolerance)
    {
        rep.converged = true;
        rep.final_residual = rep.residual_norm_0;
        return rep;
    }

    const JvpFn jvp_fn = jvp ? jvp : make_finite_difference_jvp(residual, cfg.fd_epsilon);

    class JvpOperator final : public core::IOperator
    {
    public:
        JvpOperator(core::Field x0, JvpFn fn)
            : x0_(std::move(x0)), fn_(std::move(fn))
        {
            info_.name = "newton_jvp";
            info_.inputs.push_back(core::FieldMetadata{
                .name = "v", .rank = core::FieldRank::Scalar, .components = 1});
            info_.outputs.push_back(core::FieldMetadata{
                .name = "jv", .rank = core::FieldRank::Scalar, .components = 1});
        }
        const core::OperatorInfo& info() const override { return info_; }
        bool apply(const core::Field& v, core::Field& out, core::ModelStatus& status) const override
        {
            return fn_(x0_, v, out, status);
        }
    private:
        core::Field x0_;
        JvpFn fn_;
        core::OperatorInfo info_;
    };

    for (int it = 1; it <= cfg.max_iterations; ++it)
    {
        JvpOperator jop(x, jvp_fn);
        // Solve J·d = -r
        core::Field rhs = r;
        scale(-1.0, rhs.data());
        d.assign(0.0);
        GmresConfig gcfg;
        gcfg.tolerance = cfg.tolerance * 0.1;
        gcfg.max_iterations = 200;
        auto lrep = solve_gmres(jop, rhs, d, gcfg, status);
        if (!lrep.converged && lrep.iterations == 0)
        {
            status.ok = false;
            status.error = "newton: inner linear solve failed";
            return rep;
        }
        axpy(cfg.relaxation, d.data(), x.data());
        if (!residual(x, r, status)) return rep;
        rep.iterations = it;
        rep.final_residual = norm2(r.data());
        if (rep.final_residual < cfg.tolerance * rep.residual_norm_0 ||
            rep.final_residual < 1e-14)
        {
            rep.converged = true;
            return rep;
        }
    }
    status.warnings.push_back("newton: max iterations reached without convergence");
    return rep;
}

/// Fixed-point iteration x ← g(x) until ||g(x) - x|| < tolerance.
inline NonlinearSolveReport solve_fixed_point(
    const std::function<bool(const core::Field& x, core::Field& x_new, core::ModelStatus&)>& g,
    core::Field& x, const NonlinearSolverConfig& cfg, core::ModelStatus& status)
{
    NonlinearSolveReport rep;
    const size_t n = x.size();
    core::Field x_new(core::FieldMetadata{
        .name = "x_new", .rank = core::FieldRank::Scalar, .components = 1,
        .location = core::FieldLocation::Node}, n);
    for (int it = 1; it <= cfg.max_iterations; ++it)
    {
        if (!g(x, x_new, status)) return rep;
        double diff = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            const double d = x_new.at(i) - x.at(i);
            diff += d * d;
        }
        diff = std::sqrt(diff);
        rep.iterations = it;
        rep.final_residual = diff;
        // damped update: x ← (1−r)·x + r·g(x)
        for (size_t i = 0; i < n; ++i)
            x.data()[i] = (1.0 - cfg.relaxation) * x.at(i) + cfg.relaxation * x_new.at(i);
        if (diff < cfg.tolerance)
        {
            rep.converged = true;
            return rep;
        }
    }
    status.warnings.push_back("fixed point: max iterations reached without convergence");
    return rep;
}

} // namespace exd::engine::numerics
