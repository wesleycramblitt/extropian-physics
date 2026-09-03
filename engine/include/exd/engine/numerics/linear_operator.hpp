#pragma once

// ─────────────────────────────────────────────────────
// Linear operator interface + field algebra
// (implementation_spec §35–§36).
//
// Iterative solvers operate on abstract operator
// interfaces: apply(x), apply_transpose(x), diagonal(),
// jacobian_vector_product(x, v).  `LinearOperator` is the
// aliased core::IOperator contract; the field helpers
// (dot/axpy/norm/copy) are the shared vector algebra the
// solvers are built on.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/operator.hpp>

#include <cmath>
#include <span>

namespace exd::engine::numerics {

using LinearOperator = core::IOperator;

// ── field algebra (contiguous spans; no allocation) ──

inline double dot(std::span<const double> a, std::span<const double> b)
{
    double acc = 0.0;
    for (size_t i = 0; i < a.size(); ++i) acc += a[i] * b[i];
    return acc;
}

inline double norm2(std::span<const double> a)
{
    return std::sqrt(dot(a, a));
}

inline void axpy(double alpha, std::span<const double> x, std::span<double> y)
{
    for (size_t i = 0; i < y.size(); ++i) y[i] += alpha * x[i];
}

inline void scale(double alpha, std::span<double> x)
{
    for (auto& v : x) v *= alpha;
}

inline void copy(std::span<const double> src, std::span<double> dst)
{
    std::copy(src.begin(), src.end(), dst.begin());
}

/// Operator negation: (−A)·x = −(A·x).  Composes a (negative definite)
/// Laplacian into the SPD Poisson operator −Δ that CG/BiCGSTAB require
/// (spec §35: operators compose without global assembly).
class NegatedOperator final : public core::IOperator
{
public:
    explicit NegatedOperator(const core::IOperator& inner) : inner_(inner) {}
    const core::OperatorInfo& info() const override { return inner_.info(); }
    bool apply(const core::Field& in, core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.apply(in, out, status)) return false;
        for (auto& v : out.data()) v = -v;
        return true;
    }
    bool apply_transpose(const core::Field& in, core::Field& out,
                         core::ModelStatus& status) const override
    {
        if (!inner_.apply_transpose(in, out, status)) return false;
        for (auto& v : out.data()) v = -v;
        return true;
    }
    bool diagonal(core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.diagonal(out, status)) return false;
        for (auto& v : out.data()) v = -v;
        return true;
    }
    bool jacobian_vector_product(const core::Field& x, const core::Field& v,
                                 core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.jacobian_vector_product(x, v, out, status)) return false;
        for (auto& vv : out.data()) vv = -vv;
        return true;
    }

private:
    const core::IOperator& inner_;
};

/// Scaled operator: (c·A)·x = c·(A·x) — e.g. K·(−Δ) for steady Darcy/heat
/// solves.  Composes with NegatedOperator without global assembly (§35).
class ScaledOperator final : public core::IOperator
{
public:
    ScaledOperator(const core::IOperator& inner, double scale)
        : inner_(inner), scale_(scale) {}
    const core::OperatorInfo& info() const override { return inner_.info(); }
    bool apply(const core::Field& in, core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.apply(in, out, status)) return false;
        for (auto& v : out.data()) v *= scale_;
        return true;
    }
    bool apply_transpose(const core::Field& in, core::Field& out,
                         core::ModelStatus& status) const override
    {
        if (!inner_.apply_transpose(in, out, status)) return false;
        for (auto& v : out.data()) v *= scale_;
        return true;
    }
    bool diagonal(core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.diagonal(out, status)) return false;
        for (auto& v : out.data()) v *= scale_;
        return true;
    }
    bool jacobian_vector_product(const core::Field& x, const core::Field& v,
                                 core::Field& out, core::ModelStatus& status) const override
    {
        if (!inner_.jacobian_vector_product(x, v, out, status)) return false;
        for (auto& vv : out.data()) vv *= scale_;
        return true;
    }

private:
    const core::IOperator& inner_;
    double scale_;
};

/// Absolute + relative residual criterion for iterative solvers.
struct IterativeSolverConfig
{
    double tolerance = 1e-10;   // relative residual target ||r||/||b||
    int max_iterations = 1000;
    double min_absolute_residual = 1e-14;  // early stop floor
};

struct IterativeSolveReport
{
    bool converged = false;
    int iterations = 0;
    double final_residual = 0.0;   // ||r||/||b|| (1.0 when b == 0)
};

} // namespace exd::engine::numerics
