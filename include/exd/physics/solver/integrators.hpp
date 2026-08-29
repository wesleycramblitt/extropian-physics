#pragma once

#include <exd/physics/model_status.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>

namespace exd::physics::solver {

// ─────────────────────────────────────────────────────
// Shared ODE integrators. One module, many methods:
// explicit (ForwardEuler, Heun, RK4, AdaptiveRK45), implicit
// (BackwardEuler, CrankNicolson), and geometric
// (SymplecticEuler, Verlet) for Hamiltonian systems.
//
// Works on a flat state vector + derivative callback, so
// engines, circuits, control loops, rigid bodies and any
// future domain share identical stepping semantics.
// ─────────────────────────────────────────────────────

enum class IntegrationMethod : uint8_t
{
    ForwardEuler,     // y' = f(y):           y += dt·f
    BackwardEuler,    // implicit fixed point y += dt·f(y_new)  (A-stable)
    Heun,             // 2nd-order explicit RK
    RK4,              // 4th-order explicit RK
    CrankNicolson,    // trapezoidal implicit (A-stable, 2nd order)
    SymplecticEuler,  // semi-implicit Euler; state = [positions, velocities]
    Verlet,           // leapfrog;            state = [positions, velocities]
    AdaptiveRK45,     // Dormand-Prince with error control (dt may change)
};

struct IntegratorConfig
{
    IntegrationMethod method = IntegrationMethod::RK4;

    // ── implicit fixed-point (BackwardEuler, CrankNicolson) ──
    int max_fixed_point_iterations = 50;
    double fixed_point_tolerance = 1e-10;
    double relaxation = 1.0;            // under-relaxation factor in (0, 1]

    // ── geometric methods (SymplecticEuler, Verlet) ──
    // First `position_count` components of the state are positions, the
    // remaining are velocities. 0 = state.size() / 2.
    std::size_t position_count = 0;

    // ── adaptive (AdaptiveRK45) ──
    double rel_tol = 1e-6;
    double abs_tol = 1e-9;
    double dt_min = 1e-12;
    double dt_max = 1e3;
};

/// State derivative callback: fills `dstate` from `state` at time `t`.
using DerivativeFn = std::function<void(std::span<const double> state,
                                        std::span<double> dstate, double t)>;

/// Advance `state` by one step of `dt`.
///
/// - AdaptiveRK45 may shrink dt internally and reports the accepted step in
///   `dt_used` (when non-null; state is advanced only by the accepted dt).
/// - Implicit methods iterate to `fixed_point_tolerance` with under-
///   relaxation and report failure through `status` on non-convergence
///   (the best iterate is still written to `state`).
/// - Never throws. Returns false on invalid config/state.
bool integrate_step(const IntegratorConfig& config, double t, double dt,
                    std::span<double> state, const DerivativeFn& deriv,
                    exd::physics::ModelStatus& status,
                    double* dt_used = nullptr);

} // namespace exd::physics::solver
