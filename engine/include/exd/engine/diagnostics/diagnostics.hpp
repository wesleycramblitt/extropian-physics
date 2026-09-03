#pragma once

// ─────────────────────────────────────────────────────
// Diagnostics (implementation_spec §50).
//
// Basic diagnostics belong in the engine and are
// accessible programmatically: residual, CFL, mass/energy/
// momentum conservation, solver iterations, timestep,
// runtime.  Accuracy/error estimation (spec §48) is the
// roadmap layer on top of these primitives.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace exd::engine::diagnostics {

/// CFL number of an advection/diffusion step on a uniform grid (spec §50).
struct CflDiagnostics
{
    double cfl_advection = 0.0;
    double cfl_diffusion = 0.0;
    double dt = 0.0;

    [[nodiscard]] double max_cfl() const { return std::max(cfl_advection, cfl_diffusion); }

    [[nodiscard]] bool stable(double limit = 1.0) const { return max_cfl() <= limit; }
};

inline CflDiagnostics compute_cfl(const mesh::StructuredGrid& g,
                                  double u_max, double diffusivity, double dt)
{
    CflDiagnostics out;
    out.dt = dt;
    const double h_min = std::min({g.spacing[0], g.spacing[1], g.spacing[2]});
    out.cfl_advection = u_max * dt / h_min;
    out.cfl_diffusion = diffusivity * dt / (h_min * h_min);
    return out;
}

/// Conservation ledger for a scalar quantity (mass/energy) on a domain:
/// compare a tracked buffer before/after a step.
struct ConservationDiagnostics
{
    double initial_total = 0.0;
    double final_total = 0.0;

    [[nodiscard]] double absolute_error() const { return final_total - initial_total; }
    [[nodiscard]] double relative_error() const
    {
        return initial_total == 0.0 ? 0.0 : std::abs(absolute_error()) / std::abs(initial_total);
    }
    [[nodiscard]] bool conserved(double tolerance = 1e-10) const
    {
        return relative_error() <= tolerance;
    }
};

inline double field_total(std::span<const double> values)
{
    double acc = 0.0;
    for (double v : values) acc += v;
    return acc;
}

/// Residual history for convergence diagnostics.
struct ResidualHistory
{
    std::vector<double> residuals;
    size_t iterations = 0;
    double tolerance = 0.0;

    bool converged() const { return !residuals.empty() && residuals.back() < tolerance; }
    double final_residual() const { return residuals.empty() ? 0.0 : residuals.back(); }
};

/// Simple wall-clock stopwatch for performance diagnostics (spec §50:
/// runtime, timestep rate).
class Stopwatch
{
public:
    void start() { start_ = std::chrono::steady_clock::now(); }
    double seconds() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }
    double steps_per_second(size_t steps) const
    {
        const double s = seconds();
        return s > 0.0 ? static_cast<double>(steps) / s : 0.0;
    }

private:
    std::chrono::steady_clock::time_point start_{};
};

} // namespace exd::engine::diagnostics
