#pragma once

// ─────────────────────────────────────────────────────
// Numerical operators (implementation_spec §34–§35).
//
// Operators are reusable mathematical transformations on
// fields/state: Gradient, Divergence, Curl, Laplacian,
// Advection, Diffusion, Projection, Residual, Jacobian…
// They expose requirements (OperatorInfo) and outputs.
//
// The matrix-free architecture (§35) prefers
//   Input Fields → Operator → Output Field
// over global matrix assembly.  Where a linear operator
// exists it supports: apply(x), apply_transpose(x),
// diagonal(), jacobian_vector_product(x, v) — the
// interfaces iterative solvers consume.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>
#include <exd/engine/core/model_status.hpp>

#include <string>
#include <vector>

namespace exd::engine::core {

/// Declared operator signature: name, input field metadata (requirements),
/// output field metadata (produces).
struct OperatorInfo
{
    std::string name;
    std::vector<FieldMetadata> inputs;
    std::vector<FieldMetadata> outputs;
};

/// Abstract linear/matrix-free operator on Fields (§35).
class IOperator
{
public:
    virtual ~IOperator() = default;

    virtual const OperatorInfo& info() const = 0;

    /// out = A·in
    virtual bool apply(const Field& in, Field& out, ModelStatus& status) const = 0;

    /// out = Aᵀ·in (default: unsupported unless overridden)
    virtual bool apply_transpose(const Field& in, Field& out, ModelStatus& status) const;

    /// out = diag(A) (default: unsupported)
    virtual bool diagonal(Field& out, ModelStatus& status) const;

    /// out = J(x)·v (default: unsupported)
    virtual bool jacobian_vector_product(const Field& x, const Field& v,
                                         Field& out, ModelStatus& status) const;
};

inline bool IOperator::apply_transpose(const Field&, Field&, ModelStatus& status) const
{
    status.ok = false;
    status.error = "operator '" + info().name + "': apply_transpose not implemented";
    return false;
}

inline bool IOperator::diagonal(Field&, ModelStatus& status) const
{
    status.ok = false;
    status.error = "operator '" + info().name + "': diagonal not implemented";
    return false;
}

inline bool IOperator::jacobian_vector_product(const Field&, const Field&,
                                               Field&, ModelStatus& status) const
{
    status.ok = false;
    status.error = "operator '" + info().name + "': jacobian_vector_product not implemented";
    return false;
}

} // namespace exd::engine::core
