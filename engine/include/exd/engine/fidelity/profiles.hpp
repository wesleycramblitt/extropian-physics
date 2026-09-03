#pragma once

// ─────────────────────────────────────────────────────
// Fidelity profiles (implementation_spec §45–§46).
//
// Fidelity is a first-class engine concept: the same
// simulation configuration runs at different fidelity
// levels.  Builtin profiles are defaults — users override
// any knob.  REALTIME → INTERACTIVE → FAST → BALANCED →
// ACCURATE → HIGH_FIDELITY.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/field.hpp>

#include <cstdint>

namespace exd::engine::fidelity {

enum class FidelityLevel : uint8_t
{
    Realtime,
    Interactive,
    Fast,
    Balanced,
    Accurate,
    HighFidelity,
};

constexpr const char* to_string(FidelityLevel l)
{
    switch (l)
    {
    case FidelityLevel::Realtime: return "REALTIME";
    case FidelityLevel::Interactive: return "INTERACTIVE";
    case FidelityLevel::Fast: return "FAST";
    case FidelityLevel::Balanced: return "BALANCED";
    case FidelityLevel::Accurate: return "ACCURATE";
    case FidelityLevel::HighFidelity: return "HIGH_FIDELITY";
    }
    return "?";
}

/// Fidelity controls (spec §46).  All knobs are defaults; presets and
/// solvers read them where the knob exists (unknown knobs are ignored).
struct FidelityProfile
{
    FidelityLevel level = FidelityLevel::Balanced;

    double solver_tolerance_scale = 1.0;   // multiply solver tolerances
    int max_solver_iterations_scale = 1;   // multiply solver iteration caps
    double dt_scale = 1.0;                 // multiply timestep targets
    double coupling_tolerance = 1e-8;      // coupling iteration tolerance
    int coupling_iterations = 1;           // max coupling sub-iterations
    core::FieldPrecision precision = core::FieldPrecision::F64;
    int output_cadence_scale = 1;          // multiply output interval
    bool simplified_physics = false;       // allow simplified sub-models
};

/// Builtin profiles (spec §45 names; defaults, not hardcoded modes).
inline const FidelityProfile& profile(FidelityLevel level)
{
    static const FidelityProfile realtime{
        .level = FidelityLevel::Realtime,
        .solver_tolerance_scale = 100.0,
        .max_solver_iterations_scale = 0,
        .dt_scale = 2.0,
        .coupling_tolerance = 1e-4,
        .coupling_iterations = 1,
        .precision = core::FieldPrecision::F32,
        .output_cadence_scale = 10,
        .simplified_physics = true,
    };
    static const FidelityProfile interactive{
        .level = FidelityLevel::Interactive,
        .solver_tolerance_scale = 30.0,
        .max_solver_iterations_scale = 0,
        .dt_scale = 1.5,
        .coupling_tolerance = 1e-5,
        .coupling_iterations = 1,
        .precision = core::FieldPrecision::F32,
        .output_cadence_scale = 5,
        .simplified_physics = true,
    };
    static const FidelityProfile fast{
        .level = FidelityLevel::Fast,
        .solver_tolerance_scale = 10.0,
        .max_solver_iterations_scale = 1,
        .dt_scale = 1.2,
        .coupling_tolerance = 1e-6,
        .coupling_iterations = 2,
        .precision = core::FieldPrecision::F64,
        .output_cadence_scale = 2,
        .simplified_physics = false,
    };
    static const FidelityProfile balanced{
        .level = FidelityLevel::Balanced,
        .solver_tolerance_scale = 1.0,
        .max_solver_iterations_scale = 1,
        .dt_scale = 1.0,
        .coupling_tolerance = 1e-8,
        .coupling_iterations = 10,
        .precision = core::FieldPrecision::F64,
        .output_cadence_scale = 1,
        .simplified_physics = false,
    };
    static const FidelityProfile accurate{
        .level = FidelityLevel::Accurate,
        .solver_tolerance_scale = 0.1,
        .max_solver_iterations_scale = 2,
        .dt_scale = 0.5,
        .coupling_tolerance = 1e-10,
        .coupling_iterations = 50,
        .precision = core::FieldPrecision::F64,
        .output_cadence_scale = 1,
        .simplified_physics = false,
    };
    static const FidelityProfile high_fidelity{
        .level = FidelityLevel::HighFidelity,
        .solver_tolerance_scale = 0.01,
        .max_solver_iterations_scale = 4,
        .dt_scale = 0.25,
        .coupling_tolerance = 1e-12,
        .coupling_iterations = 200,
        .precision = core::FieldPrecision::F64,
        .output_cadence_scale = 1,
        .simplified_physics = false,
    };
    switch (level)
    {
    case FidelityLevel::Realtime: return realtime;
    case FidelityLevel::Interactive: return interactive;
    case FidelityLevel::Fast: return fast;
    case FidelityLevel::Accurate: return accurate;
    case FidelityLevel::HighFidelity: return high_fidelity;
    case FidelityLevel::Balanced: return balanced;
    }
    return balanced;
}

} // namespace exd::engine::fidelity
