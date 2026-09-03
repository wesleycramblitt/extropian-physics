#pragma once

#include <array>

namespace exd::engine::physics::rigid_body {

// ─────────────────────────────────────────────────────
// Quaternion utilities for 3D orientation (Hamilton
// convention, w-first storage: {w, x, y, z}).
// ─────────────────────────────────────────────────────

/// Hamilton product a ⊗ b.
std::array<double, 4> quat_multiply(const std::array<double, 4>& a,
                                    const std::array<double, 4>& b);

/// Rotate vector `v` by quaternion `q`: world = q ⊗ (0, v) ⊗ q*.
/// Requires |q| ≈ 1 (call quat_normalize first if unsure).
std::array<double, 3> quat_rotate(const std::array<double, 4>& q,
                                  const std::array<double, 3>& v);

/// Normalize `q` in place (no-op for the zero quaternion).
void quat_normalize(std::array<double, 4>& q);

/// Integrate orientation under a constant body-frame angular velocity:
/// q_new = normalize(q + 0.5·dt·[0, ω_body] ⊗ q).
std::array<double, 4> quat_integrate_body(const std::array<double, 4>& q,
                                          const std::array<double, 3>& omega_body,
                                          double dt);

} // namespace exd::engine::physics::rigid_body
