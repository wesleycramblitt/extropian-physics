// 3D orientation utilities using the Hamilton quaternion convention with
// w-first storage {w, x, y, z}. Single-pose helpers used by the rigid-body
// integrators; explicit array math only, no external math library.

#include <exd/engine/physics/rigid_body/quaternion.hpp>

#include <array>
#include <cmath>

namespace exd::engine::physics::rigid_body
{

std::array<double, 4> quat_multiply(const std::array<double, 4>& a,
                                    const std::array<double, 4>& b)
{
    // Hamilton product a ⊗ b (w-first components).
    return {
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    };
}

std::array<double, 3> quat_rotate(const std::array<double, 4>& q,
                                  const std::array<double, 3>& v)
{
    // World vector = q ⊗ (0, v) ⊗ q*, with q* = (w, −x, −y, −z).
    // Composed directly from the Hamilton product as specified.
    const std::array<double, 4> vector_quat = {0.0, v[0], v[1], v[2]};
    const std::array<double, 4> product = quat_multiply(q, vector_quat);
    const std::array<double, 4> conjugate = {q[0], -q[1], -q[2], -q[3]};
    const std::array<double, 4> rotated = quat_multiply(product, conjugate);
    return {rotated[1], rotated[2], rotated[3]};
}

void quat_normalize(std::array<double, 4>& q)
{
    const double length =
        std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);

    // The zero quaternion has no direction; leave it unchanged (documented no-op).
    if (length < 1e-12)
        return;

    q[0] /= length;
    q[1] /= length;
    q[2] /= length;
    q[3] /= length;
}

std::array<double, 4> quat_integrate_body(const std::array<double, 4>& q,
                                          const std::array<double, 3>& omega_body,
                                          double dt)
{
    // dq = 0.5·dt·[0, ωx, ωy, ωz] ⊗ q; result = normalize(q + dq).
    const std::array<double, 4> omega = {0.0, omega_body[0], omega_body[1], omega_body[2]};
    const std::array<double, 4> dq = quat_multiply(omega, q);

    std::array<double, 4> next = {
        q[0] + 0.5 * dt * dq[0],
        q[1] + 0.5 * dt * dq[1],
        q[2] + 0.5 * dt * dq[2],
        q[3] + 0.5 * dt * dq[3],
    };
    quat_normalize(next);
    return next;
}

} // namespace exd::engine::physics::rigid_body