// Two-link planar arm — minimal-coordinate rigid-body chains with revolute
// joints, torque limits, joint stops and PD control (robotic-arm use case).

#include <exd/engine/physics/robotics/arm.hpp>

#include <cmath>

namespace exd::engine::physics::robotics {

void arm_mass(const TwoLinkArmConfig& c, double q2,
              double& m11, double& m12, double& m22)
{
    const double c2 = std::cos(q2);
    m11 = c.m1 * c.l1 * c.l1 +
          c.m2 * (c.l1 * c.l1 + c.l2 * c.l2 + 2.0 * c.l1 * c.l2 * c2);
    m12 = c.m2 * (c.l2 * c.l2 + c.l1 * c.l2 * c2);
    m22 = c.m2 * c.l2 * c.l2;
}

void arm_coriolis(const TwoLinkArmConfig& c, double q2, double dq1, double dq2,
                  double& c1, double& c2)
{
    const double s2 = std::sin(q2);
    const double h = c.m2 * c.l1 * c.l2 * s2;
    c1 = -h * (2.0 * dq1 * dq2 + dq2 * dq2);
    c2 = h * dq1 * dq1;
}

void arm_gravity(const TwoLinkArmConfig& c, double q1, double q2,
                 double& g1, double& g2)
{
    const double c1q = std::cos(q1);
    const double c12 = std::cos(q1 + q2);
    g1 = (c.m1 + c.m2) * c.l1 * c.gravity * c1q + c.m2 * c.l2 * c.gravity * c12;
    g2 = c.m2 * c.l2 * c.gravity * c12;
}

void arm_forward_kinematics(const TwoLinkArmConfig& c, double q1, double q2,
                            double& x, double& y)
{
    const double c1 = std::cos(q1), s1 = std::sin(q1);
    const double c12 = std::cos(q1 + q2), s12 = std::sin(q1 + q2);
    x = c.l1 * c1 + c.l2 * c12;
    y = c.l1 * s1 + c.l2 * s12;
}

double pd_torque(double q, double dq, double q_ref, double dq_ref,
                 const PidGains& gains, double dt, double& integral)
{
    const double e = q_ref - q;
    const double de = dq_ref - dq;
    integral += gains.ki * e * dt;
    const double clamp = gains.integral_clamp;
    integral = std::max(-clamp, std::min(clamp, integral));
    return gains.kp * e + gains.kd * de + integral;
}

double arm_kinetic_energy(const TwoLinkArmConfig& c, const ArmState& s)
{
    double m11, m12, m22;
    arm_mass(c, s.q2, m11, m12, m22);
    return 0.5 * (m11 * s.dq1 * s.dq1 + 2.0 * m12 * s.dq1 * s.dq2 +
                  m22 * s.dq2 * s.dq2);
}

namespace {

/// Effective torque at a joint: actuator torque clamped to ±t_max plus the
/// joint-stop penalty when a limit is crossed and the motion pushes further.
double effective_torque(double tau, double t_max,
                        double q, double dq, double q_min, double q_max,
                        double stiffness, double damping)
{
    double out = std::max(-t_max, std::min(t_max, tau));
    if (q_max > q_min)   // stops enabled
    {
        if (q < q_min && dq < 0.0)
            out += stiffness * (q_min - q) - damping * dq;
        else if (q > q_max && dq > 0.0)
            out += stiffness * (q_max - q) - damping * dq;
    }
    return out;
}

} // namespace

bool step_arm(const TwoLinkArmConfig& config, ArmState& state,
              double tau1, double tau2, double dt, core::ModelStatus& status)
{
    if (!(dt > 0.0))
    {
        status.ok = false;
        status.error = "arm: dt must be positive";
        return false;
    }
    // accel: q̈ = M⁻¹·(τ − C − G); RK4 on the 4-dim state (q1,q2,dq1,dq2).
    // The torque (limits + joint-stop penalty) is evaluated at the CURRENT
    // state — piecewise-constant within the step.
    const double t1 = effective_torque(tau1, config.t1_max, state.q1, state.dq1,
                                       config.q1_min, config.q1_max,
                                       config.stop_stiffness, config.stop_damping);
    const double t2 = effective_torque(tau2, config.t2_max, state.q2, state.dq2,
                                       config.q2_min, config.q2_max,
                                       config.stop_stiffness, config.stop_damping);

    auto deriv = [&](double q1, double q2, double dq1, double dq2,
                     double out[4]) -> bool {
        double m11, m12, m22;
        arm_mass(config, q2, m11, m12, m22);
        double c1, c2;
        arm_coriolis(config, q2, dq1, dq2, c1, c2);
        double g1, g2;
        arm_gravity(config, q1, q2, g1, g2);
        const double det = m11 * m22 - m12 * m12;
        if (!(det > 0.0)) return false;
        const double r1 = t1 - c1 - g1;
        const double r2 = t2 - c2 - g2;
        out[0] = dq1;
        out[1] = dq2;
        out[2] = (m22 * r1 - m12 * r2) / det;
        out[3] = (-m12 * r1 + m11 * r2) / det;
        return true;
    };

    const double y0[4] = {state.q1, state.q2, state.dq1, state.dq2};
    double k1[4], k2[4], k3[4], k4[4];
    if (!deriv(y0[0], y0[1], y0[2], y0[3], k1)) { status.ok = false; status.error = "arm: singular mass matrix"; return false; }
    if (!deriv(y0[0] + dt / 2 * k1[0], y0[1] + dt / 2 * k1[1],
               y0[2] + dt / 2 * k1[2], y0[3] + dt / 2 * k1[3], k2)) { status.ok = false; status.error = "arm: singular mass matrix"; return false; }
    if (!deriv(y0[0] + dt / 2 * k2[0], y0[1] + dt / 2 * k2[1],
               y0[2] + dt / 2 * k2[2], y0[3] + dt / 2 * k2[3], k3)) { status.ok = false; status.error = "arm: singular mass matrix"; return false; }
    if (!deriv(y0[0] + dt * k3[0], y0[1] + dt * k3[1],
               y0[2] + dt * k3[2], y0[3] + dt * k3[3], k4)) { status.ok = false; status.error = "arm: singular mass matrix"; return false; }

    const double w = dt / 6.0;
    state.q1  = y0[0] + w * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    state.q2  = y0[1] + w * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    state.dq1 = y0[2] + w * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    state.dq2 = y0[3] + w * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
    return true;
}

} // namespace exd::engine::physics::robotics
