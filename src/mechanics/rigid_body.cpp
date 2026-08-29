// 6-DOF rigid body dynamics.
// Translation is integrated in the world frame (F = m·dv/dt); rotation in the
// body frame (I·dω/dt = τ − ω × I·ω) with quaternion orientation
// (dq/dt = 0.5·[0, ω_body] ⊗ q). Two integrators are provided:
//   - SymplecticEuler: semi-implicit, bounded energy drift on Hamiltonian systems.
//   - RK4: classic fourth-order explicit, renormalizes the orientation each step.
// Loads are held constant for the duration of a step (advance() contract).

#include <exd/physics/mechanics/rigid_body.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string_view>

namespace exd::physics::mechanics
{

namespace
{

/// 13-component integration state: [p(3), q(4), v(3), ω(3)].
using RigidState = std::array<double, 13>;

std::array<double, 3> cross(const std::array<double, 3>& a,
                            const std::array<double, 3>& b)
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

RigidState to_state(const RigidBodyState& s)
{
    return {
        s.position[0], s.position[1], s.position[2],
        s.orientation[0], s.orientation[1], s.orientation[2], s.orientation[3],
        s.linear_velocity[0], s.linear_velocity[1], s.linear_velocity[2],
        s.angular_velocity[0], s.angular_velocity[1], s.angular_velocity[2],
    };
}

RigidBodyState from_state(const RigidState& s)
{
    RigidBodyState out;
    out.position = {s[0], s[1], s[2]};
    out.orientation = {s[3], s[4], s[5], s[6]};
    out.linear_velocity = {s[7], s[8], s[9]};
    out.angular_velocity = {s[10], s[11], s[12]};
    return out;
}

void normalize_orientation(RigidState& s)
{
    const double length =
        std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5] + s[6] * s[6]);
    if (length < 1e-12)
        return;
    s[3] /= length;
    s[4] /= length;
    s[5] /= length;
    s[6] /= length;
}

/// State derivative of the 13-component dynamics under constant loads:
///   d/dt[p] = v                     (world frame)
///   d/dt[q] = 0.5·[0, ω_body] ⊗ q
///   d/dt[v] = F/m                   (world frame)
///   d/dt[ω] = I⁻¹·(τ − ω × (I·ω))  (body frame)
RigidState derivative(const RigidState& s, const RigidBodyConfig& config,
                      const RigidBodyForces& loads)
{
    const std::array<double, 3> omega = {s[10], s[11], s[12]};

    // I·ω then ω × (I·ω)
    const std::array<double, 3> iw = {
        config.inertia_principal[0] * omega[0],
        config.inertia_principal[1] * omega[1],
        config.inertia_principal[2] * omega[2],
    };
    const std::array<double, 3> w_cross_iw = cross(omega, iw);

    // [0, ω] ⊗ q
    const std::array<double, 4> omega_quat = {0.0, omega[0], omega[1], omega[2]};
    const std::array<double, 4> q = {s[3], s[4], s[5], s[6]};
    const std::array<double, 4> dq = quat_multiply(omega_quat, q);

    return {
        // d/dt[p] = v
        s[7], s[8], s[9],
        // d/dt[q] = 0.5·[0, ω] ⊗ q
        0.5 * dq[0], 0.5 * dq[1], 0.5 * dq[2], 0.5 * dq[3],
        // d/dt[v] = F/m
        loads.force[0] / config.mass,
        loads.force[1] / config.mass,
        loads.force[2] / config.mass,
        // d/dt[ω] = I⁻¹·(τ − ω × (I·ω))
        (loads.torque[0] - w_cross_iw[0]) / config.inertia_principal[0],
        (loads.torque[1] - w_cross_iw[1]) / config.inertia_principal[1],
        (loads.torque[2] - w_cross_iw[2]) / config.inertia_principal[2],
    };
}

class RigidBodyDynamicsImpl final : public IRigidBodyDynamics
{
public:
    RigidBodyDynamicsImpl(RigidBodyIntegration method, RigidBodyConfig config)
        : method_(method), config_(config)
    {
    }

    std::string_view name() const override { return "rigid_body"; }

    RigidBodyState advance(double dt, const RigidBodyForces& loads,
                           const RigidBodyState& state,
                           ModelStatus& status) const override
    {
        // Validate dt, mass and each principal inertia. On failure the input
        // state is returned unchanged.
        if (dt <= 0.0)
        {
            status.ok = false;
            status.error = "rigid body: dt must be positive";
            return state;
        }
        if (config_.mass <= 0.0)
        {
            status.ok = false;
            status.error = "rigid body: mass must be positive";
            return state;
        }
        if (config_.inertia_principal[0] <= 0.0 ||
            config_.inertia_principal[1] <= 0.0 ||
            config_.inertia_principal[2] <= 0.0)
        {
            status.ok = false;
            status.error = "rigid body: principal inertias must be positive";
            return state;
        }

        if (method_ == RigidBodyIntegration::SymplecticEuler)
            return advance_symplectic_euler(dt, loads, state);

        return advance_rk4(dt, loads, state);
    }

private:
    // ── Symplectic (semi-implicit) Euler ─────────────────────────────
    //   1. v_new = v + dt·F/m            (world; F constant over the step)
    //   2. ω_new = ω + dt·I⁻¹·(τ − ω×Iω) (explicit in ω, body frame)
    //   3. q_new = quat_integrate_body(q, ω_new, dt)  (normalized)
    //   4. p_new = p + dt·v_new
    RigidBodyState advance_symplectic_euler(double dt, const RigidBodyForces& loads,
                                            const RigidBodyState& state) const
    {
        const std::array<double, 3> inv_inertia = {
            1.0 / config_.inertia_principal[0],
            1.0 / config_.inertia_principal[1],
            1.0 / config_.inertia_principal[2],
        };

        const std::array<double, 3> iw = {
            config_.inertia_principal[0] * state.angular_velocity[0],
            config_.inertia_principal[1] * state.angular_velocity[1],
            config_.inertia_principal[2] * state.angular_velocity[2],
        };
        const std::array<double, 3> w_cross_iw = cross(state.angular_velocity, iw);

        std::array<double, 3> v_new = {
            state.linear_velocity[0] + dt * loads.force[0] / config_.mass,
            state.linear_velocity[1] + dt * loads.force[1] / config_.mass,
            state.linear_velocity[2] + dt * loads.force[2] / config_.mass,
        };

        std::array<double, 3> w_new = {
            state.angular_velocity[0] + dt * (loads.torque[0] - w_cross_iw[0]) * inv_inertia[0],
            state.angular_velocity[1] + dt * (loads.torque[1] - w_cross_iw[1]) * inv_inertia[1],
            state.angular_velocity[2] + dt * (loads.torque[2] - w_cross_iw[2]) * inv_inertia[2],
        };

        std::array<double, 4> q_new = quat_integrate_body(state.orientation, w_new, dt);
        quat_normalize(q_new);

        std::array<double, 3> p_new = {
            state.position[0] + dt * v_new[0],
            state.position[1] + dt * v_new[1],
            state.position[2] + dt * v_new[2],
        };

        return {p_new, q_new, v_new, w_new};
    }

    // ── Classic RK4 over the 13-component state ──────────────────────
    //   Stage orientations are normalized before each derivative
    //   evaluation (quaternions drift); the resulting orientation is
    //   normalized once more after the step.
    RigidBodyState advance_rk4(double dt, const RigidBodyForces& loads,
                               const RigidBodyState& state) const
    {
        const RigidState s0 = to_state(state);

        const RigidState k1 = derivative(s0, config_, loads);

        RigidState s1 = s0;
        for (std::size_t i = 0; i < s1.size(); ++i)
            s1[i] += 0.5 * dt * k1[i];
        normalize_orientation(s1);
        const RigidState k2 = derivative(s1, config_, loads);

        RigidState s2 = s0;
        for (std::size_t i = 0; i < s2.size(); ++i)
            s2[i] += 0.5 * dt * k2[i];
        normalize_orientation(s2);
        const RigidState k3 = derivative(s2, config_, loads);

        RigidState s3 = s0;
        for (std::size_t i = 0; i < s3.size(); ++i)
            s3[i] += dt * k3[i];
        normalize_orientation(s3);
        const RigidState k4 = derivative(s3, config_, loads);

        RigidState next = s0;
        for (std::size_t i = 0; i < next.size(); ++i)
            next[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        normalize_orientation(next);

        return from_state(next);
    }

    RigidBodyIntegration method_;
    RigidBodyConfig config_;
};

} // anonymous namespace

// ── Factory functions ──────────────────────────────────────────────

std::unique_ptr<IRigidBodyDynamics> make_rigid_body_dynamics(
    RigidBodyIntegration method, const RigidBodyConfig& config)
{
    return std::make_unique<RigidBodyDynamicsImpl>(method, config);
}

// ── Convenience load helpers (world-frame force/torque sources) ────

std::array<double, 3> force_gravity(double mass, const std::array<double, 3>& g)
{
    return {
        mass * g[0],
        mass * g[1],
        mass * g[2],
    };
}

std::array<double, 3> force_linear_spring(const std::array<double, 3>& position,
                                          const std::array<double, 3>& anchor,
                                          double stiffness)
{
    return {
        -stiffness * (position[0] - anchor[0]),
        -stiffness * (position[1] - anchor[1]),
        -stiffness * (position[2] - anchor[2]),
    };
}

} // namespace exd::physics::mechanics