#pragma once

// ─────────────────────────────────────────────────────
// Output cadence policy.
//
// Two independent triggers, OR-ed:
//   every_n_steps      — deterministic, offline replay
//   wall_clock_interval_s — "real-time if specified";
//                          throttles to a stable animation
//                          cadence regardless of step size
//
// The policy struct is PURE (no mutable state) so it stays
// optimizer-batchable; drivers keep the last-emit book-
// keeping in an OutputScheduler.
// ─────────────────────────────────────────────────────

#include <cstdint>

namespace exd::engine::output {

/// Pure trigger rules.
struct OutputPolicy
{
    uint32_t every_n_steps = 0;          // 0 = disabled
    double wall_clock_interval_s = 0.0;  // 0 = disabled

    [[nodiscard]] bool enabled() const { return every_n_steps > 0 || wall_clock_interval_s > 0; }

    /// Pure decision for step-based triggering.
    [[nodiscard]] bool step_triggered(uint64_t step) const
    {
        return every_n_steps > 0 && (step % every_n_steps) == 0;
    }
};

/// Mutable throttle state around an OutputPolicy. `now_s`
/// comes from an injected clock so tests are deterministic.
class OutputScheduler
{
public:
    explicit OutputScheduler(OutputPolicy policy,
                             double initial_now_s = 0.0,
                             double initial_last_emit_s = -1e300)
        : policy_(policy), now_s_(initial_now_s), last_emit_s_(initial_last_emit_s)
    {
    }

    /// Advance the clock (wall time for real-time mode).
    void set_now(double now_s) { now_s_ = now_s; }

    /// Should a stamp be emitted for `step` at the current
    /// clock time? Updates last-emit bookkeeping when yes.
    bool should_emit(uint64_t step)
    {
        if (policy_.step_triggered(step))
        {
            last_emit_s_ = now_s_;
            return true;
        }
        if (policy_.wall_clock_interval_s > 0.0
            && (now_s_ - last_emit_s_) >= policy_.wall_clock_interval_s)
        {
            last_emit_s_ = now_s_;
            return true;
        }
        return false;
    }

    [[nodiscard]] const OutputPolicy& policy() const { return policy_; }

private:
    OutputPolicy policy_;
    double now_s_;
    double last_emit_s_;
};

} // namespace exd::engine::output
