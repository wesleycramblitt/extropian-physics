// Serial manipulator — N-link planar chain in minimal coordinates:
// Jacobian-based mass matrix, Christoffel Coriolis, gravity, torque
// limits, joint stops, RK4 integration, PD control.

#include <exd/engine/physics/robotics/manipulator.hpp>

#include <cmath>

namespace exd::engine::physics::robotics {

namespace {

/// Angle of link m (cumulative) and its (cos, sin).
inline void link_angle(const std::vector<double>& q, int m, double& th, double& ct, double& st)
{
    th = 0.0;
    for (int j = 0; j <= m; ++j) th += q[static_cast<size_t>(j)];
    ct = std::cos(th);
    st = std::sin(th);
}

/// Jacobian of the position of the distal point of link i (2×n), flattened
/// row-major: J[i][2*n + k] etc.
void point_jacobian(const SerialManipulatorConfig& c, const std::vector<double>& q,
                    int i, std::vector<double>& jac)
{
    const int n = static_cast<int>(c.links.size());
    jac.assign(static_cast<size_t>(2 * n), 0.0);
    // p_i = Σ_{m ≤ i} l_m·(cos θ_m, sin θ_m); ∂p_i/∂q_k = Σ_{m = k..i} l_m·(−sin θ_m, cos θ_m)
    for (int k = 0; k <= i; ++k)
    {
        double th, ct, st;
        link_angle(q, k, th, ct, st);
        double sx = 0.0, sy = 0.0;
        for (int m = k; m <= i; ++m)
        {
            double thm, ctm, stm;
            link_angle(q, m, thm, ctm, stm);
            sx += -c.links[static_cast<size_t>(m)].length * stm;
            sy += c.links[static_cast<size_t>(m)].length * ctm;
        }
        jac[static_cast<size_t>(0 * n + k)] = sx;
        jac[static_cast<size_t>(1 * n + k)] = sy;
    }
}

} // namespace

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

bool validate_manipulator(const SerialManipulatorConfig& config,
                          exd::engine::core::ModelStatus& status)
{
    const size_t n = config.links.size();
    if (n == 0)
    {
        status.ok = false;
        status.error = "manipulator: at least one link required";
        return false;
    }
    auto check = [&](const std::vector<double>& v, const char* what) {
        if (!v.empty() && v.size() != n)
        {
            status.ok = false;
            status.error = std::string("manipulator: ") + what + " must be empty or per-joint";
            return false;
        }
        return true;
    };
    if (!check(config.q_min, "q_min") || !check(config.q_max, "q_max") ||
        !check(config.torque_max, "torque_max"))
        return false;
    for (auto& l : config.links)
        if (!(l.length > 0.0) || !(l.mass >= 0.0))
        {
            status.ok = false;
            status.error = "manipulator: link length must be > 0, mass >= 0";
            return false;
        }
    return true;
}

void mass_matrix(const SerialManipulatorConfig& c, const std::vector<double>& q,
                 std::vector<double>& m)
{
    const int n = static_cast<int>(c.links.size());
    m.assign(static_cast<size_t>(n * n), 0.0);
    std::vector<double> jac;
    for (int i = 0; i < n; ++i)
    {
        point_jacobian(c, q, i, jac);
        const double mi = c.links[static_cast<size_t>(i)].mass;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                m[static_cast<size_t>(a * n + b)] +=
                    mi * (jac[static_cast<size_t>(0 * n + a)] * jac[static_cast<size_t>(0 * n + b)] +
                          jac[static_cast<size_t>(1 * n + a)] * jac[static_cast<size_t>(1 * n + b)]);
    }
}

void gravity(const SerialManipulatorConfig& c, const std::vector<double>& q,
             std::vector<double>& g)
{
    const int n = static_cast<int>(c.links.size());
    g.assign(static_cast<size_t>(n), 0.0);
    if (c.gravity == 0.0) return;
    std::vector<double> jac;
    for (int i = 0; i < n; ++i)
    {
        point_jacobian(c, q, i, jac);
        const double mi = c.links[static_cast<size_t>(i)].mass;
        for (int k = 0; k < n; ++k)
            g[static_cast<size_t>(k)] += mi * c.gravity *
                jac[static_cast<size_t>(1 * n + k)];
    }
}

void coriolis(const SerialManipulatorConfig& c, const std::vector<double>& q,
              const std::vector<double>& dq, std::vector<double>& cc)
{
    // Christoffel symbols from finite-difference ∂M/∂q:
    //   c_ijk = ½·(∂M_ij/∂q_k + ∂M_ik/∂q_j − ∂M_jk/∂q_i)
    //   C_j = Σ_{i,k} c_ijk · q̇_i · q̇_k
    const int n = static_cast<int>(c.links.size());
    cc.assign(static_cast<size_t>(n), 0.0);
    const double eps = 1e-6;
    std::vector<double> mp, mm;
    std::vector<std::vector<double>> dMdq(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k)
    {
        std::vector<double> qp = q, qm = q;
        qp[static_cast<size_t>(k)] += eps;
        qm[static_cast<size_t>(k)] -= eps;
        mass_matrix(c, qp, mp);
        mass_matrix(c, qm, mm);
        dMdq[static_cast<size_t>(k)].assign(static_cast<size_t>(n * n), 0.0);
        for (int idx = 0; idx < n * n; ++idx)
            dMdq[static_cast<size_t>(k)][static_cast<size_t>(idx)] =
                (mp[static_cast<size_t>(idx)] - mm[static_cast<size_t>(idx)]) / (2.0 * eps);
    }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k)
            {
                const double c_ijk = 0.5 *
                    (dMdq[static_cast<size_t>(k)][static_cast<size_t>(j * n + i)] +
                     dMdq[static_cast<size_t>(j)][static_cast<size_t>(i * n + k)] -
                     dMdq[static_cast<size_t>(i)][static_cast<size_t>(j * n + k)]);
                cc[static_cast<size_t>(j)] += c_ijk * dq[static_cast<size_t>(i)] *
                                              dq[static_cast<size_t>(k)];
            }
}

void forward_kinematics(const SerialManipulatorConfig& c, const std::vector<double>& q,
                        double& x, double& y)
{
    x = 0.0;
    y = 0.0;
    for (int m = 0; m < static_cast<int>(c.links.size()); ++m)
    {
        double th, ct, st;
        link_angle(q, m, th, ct, st);
        x += c.links[static_cast<size_t>(m)].length * ct;
        y += c.links[static_cast<size_t>(m)].length * st;
    }
}

double kinetic_energy(const SerialManipulatorConfig& c, const ManipulatorState& s)
{
    std::vector<double> m;
    mass_matrix(c, s.q, m);
    const int n = static_cast<int>(c.links.size());
    double ke = 0.0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            ke += 0.5 * m[static_cast<size_t>(a * n + b)] *
                  s.dq[static_cast<size_t>(a)] * s.dq[static_cast<size_t>(b)];
    return ke;
}

namespace {

/// Effective torque at a joint: actuator torque clamped to ±torque_max plus
/// the joint-stop penalty when a limit is crossed and the motion pushes
/// further outward.
double effective_torque(double tau, double t_max, bool limited,
                        double q, double dq, double q_min, double q_max,
                        bool has_stop, double stiffness, double damping)
{
    double out = limited ? std::max(-t_max, std::min(t_max, tau)) : tau;
    if (has_stop)
    {
        if (q < q_min && dq < 0.0)
            out += stiffness * (q_min - q) - damping * dq;
        else if (q > q_max && dq > 0.0)
            out += stiffness * (q_max - q) - damping * dq;
    }
    return out;
}

} // namespace

bool step_manipulator(const SerialManipulatorConfig& config, ManipulatorState& state,
                      const std::vector<double>& tau, double dt,
                      exd::engine::core::ModelStatus& status)
{
    if (!validate_manipulator(config, status)) return false;
    const size_t n = config.links.size();
    if (state.q.size() != n || state.dq.size() != n || tau.size() != n)
    {
        status.ok = false;
        status.error = "manipulator: state/torque size mismatch";
        return false;
    }
    if (!(dt > 0.0))
    {
        status.ok = false;
        status.error = "manipulator: dt must be positive";
        return false;
    }

    // torque limits + joint stops at the CURRENT state (piecewise-constant
    // within the step)
    std::vector<double> t_eff(n);
    for (size_t j = 0; j < n; ++j)
    {
        const bool limited = !config.torque_max.empty();
        const double tmax = limited ? config.torque_max[j] : 0.0;
        const bool has_stop = !config.q_min.empty() && !config.q_max.empty() &&
                              config.q_max[j] > config.q_min[j];
        const double qmin = has_stop ? config.q_min[j] : 0.0;
        const double qmax = has_stop ? config.q_max[j] : 0.0;
        t_eff[j] = effective_torque(tau[j], tmax, limited,
                                    state.q[j], state.dq[j], qmin, qmax,
                                    has_stop, config.stop_stiffness, config.stop_damping);
    }

    auto deriv = [&](const std::vector<double>& q, const std::vector<double>& dq,
                     std::vector<double>& out) -> bool {
        out.resize(2 * n);
        std::vector<double> m, cc, g;
        mass_matrix(config, q, m);
        coriolis(config, q, dq, cc);
        gravity(config, q, g);
        // solve M·q̈ = τ − C − G (dense Gaussian elimination, n small)
        std::vector<double> a = m;
        std::vector<double> r(n);
        for (size_t j = 0; j < n; ++j)
            r[j] = t_eff[j] - cc[j] - g[j];
        for (size_t col = 0; col < n; ++col)
        {
            size_t piv = col;
            for (size_t row = col + 1; row < n; ++row)
                if (std::fabs(a[row * n + col]) > std::fabs(a[piv * n + col]))
                    piv = row;
            if (std::fabs(a[piv * n + col]) < 1e-14) return false;   // singular
            if (piv != col)
            {
                for (size_t c2 = 0; c2 < n; ++c2) std::swap(a[col * n + c2], a[piv * n + c2]);
                std::swap(r[col], r[piv]);
            }
            const double piv_val = a[col * n + col];
            for (size_t c2 = col; c2 < n; ++c2) a[col * n + c2] /= piv_val;
            r[col] /= piv_val;
            for (size_t row = 0; row < n; ++row)
            {
                if (row == col) continue;
                const double f = a[row * n + col];
                for (size_t c2 = col; c2 < n; ++c2) a[row * n + c2] -= f * a[col * n + c2];
                r[row] -= f * r[col];
            }
        }
        for (size_t j = 0; j < n; ++j) { out[j] = dq[j]; out[n + j] = r[j]; }
        return true;
    };

    std::vector<double> y(2 * n), k1(2 * n), k2(2 * n), k3(2 * n), k4(2 * n), tmp;
    for (size_t j = 0; j < n; ++j) { y[j] = state.q[j]; y[n + j] = state.dq[j]; }
    if (!deriv(state.q, state.dq, k1)) { status.ok = false; status.error = "manipulator: singular mass matrix"; return false; }
    auto y_plus = [&](const std::vector<double>& kk, std::vector<double>& q, std::vector<double>& dq, double h) {
        q.resize(n); dq.resize(n);
        for (size_t j = 0; j < n; ++j) { q[j] = y[j] + h * kk[j]; dq[j] = y[n + j] + h * kk[n + j]; }
    };
    {
        std::vector<double> q2, dq2;
        y_plus(k1, q2, dq2, dt / 2);
        if (!deriv(q2, dq2, k2)) { status.ok = false; status.error = "manipulator: singular mass matrix"; return false; }
        y_plus(k2, q2, dq2, dt / 2);
        if (!deriv(q2, dq2, k3)) { status.ok = false; status.error = "manipulator: singular mass matrix"; return false; }
        y_plus(k3, q2, dq2, dt);
        if (!deriv(q2, dq2, k4)) { status.ok = false; status.error = "manipulator: singular mass matrix"; return false; }
    }
    const double w = dt / 6.0;
    for (size_t j = 0; j < n; ++j)
    {
        state.q[j] = y[j] + w * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        state.dq[j] = y[n + j] + w * (k1[n + j] + 2 * k2[n + j] + 2 * k3[n + j] + k4[n + j]);
    }
    return true;
}

} // namespace exd::engine::physics::robotics
