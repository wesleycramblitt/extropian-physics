// integrators.cpp
// Shared ODE integrators (Phase A). Implements integrate_step() for the
// exd::engine::numerics module: explicit (ForwardEuler, Heun, RK4, AdaptiveRK45),
// implicit fixed-point (BackwardEuler, CrankNicolson) and geometric
// (SymplecticEuler, Verlet) stepping over a flat state vector driven by a
// single derivative callback.

#include <exd/engine/numerics/integrators.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace exd::engine::numerics {

namespace {
// ─────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────

constexpr int kMaxAdaptiveAttempts = 20;
constexpr double kAdaptiveSafety = 0.9;
constexpr double kAdaptiveMaxGrowth = 5.0;

/// Max |a_i - b_i| over the common length.
double max_abs_error(std::span<const double> a, std::span<const double> b)
{
    double max_err = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const double err = std::fabs(a[i] - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

/// Max |a_i|.
double max_abs(std::span<const double> a)
{
    double max_v = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const double v = std::fabs(a[i]);
        if (v > max_v) max_v = v;
    }
    return max_v;
}

/// Scale-aware fixed-point convergence check:
/// max|y_new - y| < tol * (1 + max|y_new|).
bool fixed_point_converged(std::span<const double> y_new,
                           std::span<const double> y, double tol)
{
    return max_abs_error(y_new, y) < tol * (1.0 + max_abs(y_new));
}

} // namespace

bool integrate_step(const IntegratorConfig& config, double t, double dt,
                    std::span<double> state, const DerivativeFn& deriv,
                    exd::engine::core::ModelStatus& status, double* dt_used)
{
    status = exd::engine::core::ModelStatus{};
    if (dt_used != nullptr) *dt_used = 0.0;

    // ── Validation (state untouched on failure) ───────────────────────────
    if (state.empty())
    {
        status.ok = false;
        status.error = "integrator: state must not be empty";
        return false;
    }
    if (!deriv)
    {
        status.ok = false;
        status.error = "integrator: derivative callback is not callable";
        return false;
    }
    if (!(dt > 0.0))
    {
        status.ok = false;
        status.error = "integrator: dt must be positive";
        return false;
    }
    if (config.position_count > state.size())
    {
        status.ok = false;
        status.error = "integrator: position_count exceeds state size";
        return false;
    }
    if (!(config.relaxation > 0.0) || config.relaxation > 1.0)
    {
        status.ok = false;
        status.error = "integrator: relaxation must be in (0, 1]";
        return false;
    }
    if (config.max_fixed_point_iterations <= 0)
    {
        status.ok = false;
        status.error = "integrator: max_fixed_point_iterations must be positive";
        return false;
    }
    if (config.method == IntegrationMethod::AdaptiveRK45)
    {
        if (!(config.dt_min > 0.0))
        {
            status.ok = false;
            status.error = "integrator: dt_min must be positive";
            return false;
        }
        if (config.dt_max < config.dt_min)
        {
            status.ok = false;
            status.error = "integrator: dt_max must be >= dt_min";
            return false;
        }
    }

    const std::size_t n = state.size();

    // Reusable scratch buffers, allocated once per call (no heap churn inside
    // the stage loops). Each buffer is size n.
    std::vector<double> buf0(n);
    std::vector<double> buf1(n);
    std::vector<double> buf2(n);
    std::vector<double> buf3(n);
    std::vector<double> buf4(n);
    std::vector<double> buf5(n);
    std::vector<double> buf6(n);
    std::vector<double> buf7(n);

    switch (config.method)
    {
        case IntegrationMethod::ForwardEuler:
        {
            // k = deriv(state, t); state += dt * k
            deriv(state, buf0, t);
            for (std::size_t i = 0; i < n; ++i) state[i] += dt * buf0[i];
        }
        break;

        case IntegrationMethod::Heun:
        {
            // k1 = deriv(state, t); pred = state + dt*k1; k2 = deriv(pred, t+dt)
            deriv(state, buf0, t);                          // buf0 = k1
            for (std::size_t i = 0; i < n; ++i) buf1[i] = state[i] + dt * buf0[i];
            deriv(buf1, buf2, t + dt);                      // buf2 = k2
            for (std::size_t i = 0; i < n; ++i) state[i] += 0.5 * dt * (buf0[i] + buf2[i]);
        }
        break;

        case IntegrationMethod::RK4:
        {
            // Classic 4-stage Runge-Kutta with half-step intermediate states.
            deriv(state, buf0, t);                          // buf0 = k1
            for (std::size_t i = 0; i < n; ++i) buf1[i] = state[i] + 0.5 * dt * buf0[i];
            deriv(buf1, buf2, t + 0.5 * dt);                // buf2 = k2
            for (std::size_t i = 0; i < n; ++i) buf1[i] = state[i] + 0.5 * dt * buf2[i];
            deriv(buf1, buf3, t + 0.5 * dt);                // buf3 = k3
            for (std::size_t i = 0; i < n; ++i) buf1[i] = state[i] + dt * buf3[i];
            deriv(buf1, buf4, t + dt);                      // buf4 = k4
            for (std::size_t i = 0; i < n; ++i)
                state[i] += dt / 6.0 * (buf0[i] + 2.0 * buf2[i] + 2.0 * buf3[i] + buf4[i]);
        }
        break;

        case IntegrationMethod::BackwardEuler:
        {
            // Fixed-point iteration: y* = state + dt * f(t+dt, y*).
            // Relaxation damps the iteration when |lambda*dt| > 1.
            for (std::size_t i = 0; i < n; ++i) buf0[i] = state[i]; // buf0 = y
            bool converged = false;
            const double relax = config.relaxation;
            for (int iter = 0; iter < config.max_fixed_point_iterations; ++iter)
            {
                deriv(buf0, buf1, t + dt);                  // buf1 = f(y, t+dt)
                for (std::size_t i = 0; i < n; ++i)
                    buf2[i] = state[i] + dt * buf1[i];      // buf2 = y_new
                const bool done =
                    fixed_point_converged(buf2, buf0, config.fixed_point_tolerance);
                for (std::size_t i = 0; i < n; ++i)
                    buf0[i] = relax * buf2[i] + (1.0 - relax) * buf0[i];
                if (done)
                {
                    converged = true;
                    break;
                }
            }
            // Best iterate (converged or not) is written to state.
            for (std::size_t i = 0; i < n; ++i) state[i] = buf0[i];
            if (!converged)
            {
                status.ok = false;
                status.error = "BackwardEuler: fixed point not converged";
                return false;
            }
        }
        break;

        case IntegrationMethod::CrankNicolson:
        {
            // Fixed-point iteration:
            //   y* = state + dt/2 * (f(state, t) + f(y*, t+dt)).
            for (std::size_t i = 0; i < n; ++i) buf1[i] = state[i]; // buf1 = y
            deriv(state, buf0, t);                          // buf0 = fn (once)
            bool converged = false;
            const double relax = config.relaxation;
            const double half_dt = 0.5 * dt;
            for (int iter = 0; iter < config.max_fixed_point_iterations; ++iter)
            {
                deriv(buf1, buf2, t + dt);                  // buf2 = f(y, t+dt)
                for (std::size_t i = 0; i < n; ++i)
                    buf3[i] = state[i] + half_dt * (buf0[i] + buf2[i]); // buf3 = y_new
                const bool done =
                    fixed_point_converged(buf3, buf1, config.fixed_point_tolerance);
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = relax * buf3[i] + (1.0 - relax) * buf1[i];
                if (done)
                {
                    converged = true;
                    break;
                }
            }
            for (std::size_t i = 0; i < n; ++i) state[i] = buf1[i];
            if (!converged)
            {
                status.ok = false;
                status.error = "CrankNicolson: fixed point not converged";
                return false;
            }
        }
        break;

        case IntegrationMethod::SymplecticEuler:
        {
            // State = [positions (pc), velocities (n-pc)].
            const std::size_t pc =
                config.position_count != 0 ? config.position_count : n / 2;
            if (pc == 0 || pc > n - pc)
            {
                status.ok = false;
                status.error = "SymplecticEuler: state must hold position and velocity blocks";
                return false;
            }
            deriv(state, buf0, t);                          // buf0 = k
            // Positions derivative unused: only the acceleration (velocity
            // block) drives the update.
            for (std::size_t i = pc; i < n; ++i) state[i] += dt * buf0[i]; // v_new
            for (std::size_t i = 0; i < pc; ++i) state[i] += dt * state[pc + i]; // q_new
        }
        break;

        case IntegrationMethod::Verlet:
        {
            // Leapfrog: v_half -> q_new at t, then k2 at t + dt/2 -> v_new.
            const std::size_t pc =
                config.position_count != 0 ? config.position_count : n / 2;
            if (pc == 0 || pc > n - pc)
            {
                status.ok = false;
                status.error = "Verlet: state must hold position and velocity blocks";
                return false;
            }
            const std::size_t vc = n - pc;
            deriv(state, buf0, t);                          // buf0 = k1
            for (std::size_t i = 0; i < vc; ++i)
                buf1[pc + i] = state[pc + i] + 0.5 * dt * buf0[pc + i]; // v_half
            for (std::size_t i = 0; i < pc; ++i)
                buf1[i] = state[i] + dt * buf1[pc + i];     // q_new
            deriv(buf1, buf2, t + 0.5 * dt);                // buf2 = k2
            for (std::size_t i = 0; i < vc; ++i)
                buf1[pc + i] += 0.5 * dt * buf2[pc + i];    // v_new
            for (std::size_t i = 0; i < n; ++i) state[i] = buf1[i];
        }
        break;

        case IntegrationMethod::AdaptiveRK45:
        {
            // ── Dormand-Prince 5(4) tableau (RK5 with embedded RK4) ──
            // Node times c2..c7.
            constexpr double c2 = 1.0 / 5.0;
            constexpr double c3 = 3.0 / 10.0;
            constexpr double c4 = 4.0 / 5.0;
            constexpr double c5 = 8.0 / 9.0;
            constexpr double c6 = 1.0;
            constexpr double c7 = 1.0;
            // Stage matrix rows a_s,j (a72 = 0 omitted below).
            constexpr double a21 = 1.0 / 5.0;
            constexpr double a31 = 3.0 / 40.0;
            constexpr double a32 = 9.0 / 40.0;
            constexpr double a41 = 44.0 / 45.0;
            constexpr double a42 = -56.0 / 15.0;
            constexpr double a43 = 32.0 / 9.0;
            constexpr double a51 = 19372.0 / 6561.0;
            constexpr double a52 = -25360.0 / 2187.0;
            constexpr double a53 = 64448.0 / 6561.0;
            constexpr double a54 = -212.0 / 729.0;
            constexpr double a61 = 9017.0 / 3168.0;
            constexpr double a62 = -355.0 / 33.0;
            constexpr double a63 = 46732.0 / 5247.0;
            constexpr double a64 = 49.0 / 176.0;
            constexpr double a65 = -5103.0 / 18656.0;
            constexpr double a71 = 35.0 / 384.0;
            constexpr double a73 = 500.0 / 1113.0;
            constexpr double a74 = 125.0 / 192.0;
            constexpr double a75 = -2187.0 / 6784.0;
            constexpr double a76 = 11.0 / 84.0;
            // 5th-order weights (k2 and k7 are not used: b5[1] = b5[6] = 0).
            constexpr double b5[7] = {35.0 / 384.0, 0.0, 500.0 / 1113.0,
                                      125.0 / 192.0, -2187.0 / 6784.0,
                                      11.0 / 84.0, 0.0};
            // Embedded 4th-order weights for the error estimate.
            constexpr double b4[7] = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0,
                                      393.0 / 640.0, -92097.0 / 339200.0,
                                      187.0 / 2100.0, 1.0 / 40.0};

            const double rel_tol = config.rel_tol;
            const double abs_tol = config.abs_tol;
            const double dt_min = config.dt_min;
            const double dt_max = config.dt_max;

            double h = dt;          // current attempt size (may shrink)
            bool accepted = false;

            for (int attempt = 0; attempt < kMaxAdaptiveAttempts; ++attempt)
            {
                // Stage derivatives k1..k7 live in buf0/buf2..buf7; buf1 is the
                // per-stage scratch state.
                deriv(state, buf0, t);                      // k1
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a21 * buf0[i]);
                deriv(buf1, buf2, t + c2 * h);              // k2
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a31 * buf0[i] + a32 * buf2[i]);
                deriv(buf1, buf3, t + c3 * h);              // k3
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a41 * buf0[i] + a42 * buf2[i] + a43 * buf3[i]);
                deriv(buf1, buf4, t + c4 * h);              // k4
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a51 * buf0[i] + a52 * buf2[i] + a53 * buf3[i]
                                              + a54 * buf4[i]);
                deriv(buf1, buf5, t + c5 * h);              // k5
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a61 * buf0[i] + a62 * buf2[i] + a63 * buf3[i]
                                              + a64 * buf4[i] + a65 * buf5[i]);
                deriv(buf1, buf6, t + c6 * h);              // k6
                for (std::size_t i = 0; i < n; ++i)
                    buf1[i] = state[i] + h * (a71 * buf0[i] + a73 * buf3[i] + a74 * buf4[i]
                                              + a75 * buf5[i] + a76 * buf6[i]);
                deriv(buf1, buf7, t + c7 * h);              // k7

                // 5th-order solution y5 (kept in buf1) plus embedded 4th-order
                // y4 used for the relative error metric.
                double err_max = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    const double y5_i = state[i] + h * (b5[0] * buf0[i] + b5[2] * buf3[i]
                                                        + b5[3] * buf4[i] + b5[4] * buf5[i]
                                                        + b5[5] * buf6[i]);
                    const double y4_i = state[i] + h * (b4[0] * buf0[i] + b4[2] * buf3[i]
                                                        + b4[3] * buf4[i] + b4[4] * buf5[i]
                                                        + b4[5] * buf6[i] + b4[6] * buf7[i]);
                    buf1[i] = y5_i;
                    const double scale =
                        abs_tol + rel_tol * std::max(std::fabs(y5_i), std::fabs(y4_i));
                    const double err_i = std::fabs(y5_i - y4_i) / scale;
                    if (err_i > err_max) err_max = err_i;
                }

                if (err_max <= 1.0)
                {
                    // Accepted: advance by the step actually used. err_max == 0
                    // (exact step) also lands here.
                    for (std::size_t i = 0; i < n; ++i) state[i] = buf1[i];
                    if (dt_used != nullptr) *dt_used = h;
                    // Suggested next step, kept internal: dt_used reports the
                    // accepted step; the next call re-proposes dt from the caller.
                    const double err_exp = err_max > 0.0 ? std::pow(err_max, -0.2)
                                                         : kAdaptiveMaxGrowth;
                    double next_dt = kAdaptiveSafety * h * err_exp;
                    next_dt = std::clamp(next_dt, dt_min, dt_max);
                    if (next_dt > kAdaptiveMaxGrowth * h) next_dt = kAdaptiveMaxGrowth * h;
                    (void)next_dt;
                    accepted = true;
                    break;
                }

                // Rejected: shrink the step and retry.
                h = std::clamp(kAdaptiveSafety * h * std::pow(err_max, -0.2),
                               dt_min, dt_max);
            }

            if (!accepted)
            {
                status.ok = false;
                status.error = "AdaptiveRK45: step rejected 20 times";
                return false;
            }
        }
        break;
    }

    return true;
}

} // namespace exd::engine::numerics