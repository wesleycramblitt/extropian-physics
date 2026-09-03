// elasticity.cpp
// Static linear elasticity (Navier-Cauchy) in displacement form with central
// finite differences, ghost-node Robin free-surface treatment and component-wise
// SOR relaxation.  See the header for the model and numerics write-up.

#include <exd/engine/physics/structural/elasticity.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace exd::engine::physics::structural {

namespace {

constexpr double kSorOmega = 1.5;

std::size_t node_idx(int i, int j, int k, int nx, int ny)
{
    return static_cast<std::size_t>(i)
           + static_cast<std::size_t>(nx)
               * (static_cast<std::size_t>(j)
                  + static_cast<std::size_t>(ny) * static_cast<std::size_t>(k));
}

int clamp_int(int p, int n)
{
    return p < 0 ? 0 : (p >= n ? n - 1 : p);
}

/// Discretized strain tensor components.
struct Strain
{
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
};

// ─────────────────────────────────────────────────────────────────────
// Shared solve context.  All stencil lookups are ghost-aware: out-of-bounds
// points first mirror the displacement across the boundary node and then add
// the traction-Robin face corrections (one per out-of-bounds face; edge and
// corner ghosts accumulate the participating faces).
// ─────────────────────────────────────────────────────────────────────

class Solver
{
public:
    Solver(const ElasticityConfig& cfg, std::vector<double>& u, const std::vector<double>& dT,
           bool mirror_only = false)
        : cfg_(cfg), u_(u), dT_(dT), mirror_only_(mirror_only)
    {
        nx_ = cfg.grid.dims[0];
        ny_ = cfg.grid.dims[1];
        nz_ = cfg.grid.dims[2];
        dx_ = cfg.grid.spacing[0];
        dy_ = cfg.grid.spacing[1];
        dz_ = cfg.grid.spacing[2];
        const double E = cfg.material.elastic_modulus;
        const double nu = cfg.material.poisson_ratio;
        lam_ = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu_ = E / (2.0 * (1.0 + nu));
        kappa_ = E / (1.0 - 2.0 * nu); // 3*lam + 2*mu
        alpha_ = cfg.thermal_expansion_coefficient;
        have_temp_ = (alpha_ > 0.0) && !dT_.empty();
        N_ = static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_)
             * static_cast<std::size_t>(nz_);
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    std::size_t count() const { return N_; }

    double lam() const { return lam_; }
    double mu() const { return mu_; }
    double kappa() const { return kappa_; }
    bool mirror_only() const { return mirror_only_; }

    double value(int comp, int i, int j, int k) const
    {
        if (in(i, j, k)) return raw(comp, i, j, k);
        const int ir = reflect(i, nx_);
        const int jr = reflect(j, ny_);
        const int kr = reflect(k, nz_);
        double v = raw(comp, ir, jr, kr);
        if (mirror_only_) return v;   // dynamic free surface: pure mirror ghost
        // Count out-of-bounds faces: tangential (shear) corrections are only
        // applied for single-face ghosts; edge/corner ghosts mirror with only
        // their normal-component corrections, which keeps the coupled shear
        // extrapolation out of the corners (patch tests are unaffected).
        const bool tang = (i >= nx_) + (i < 0) + (j >= ny_) + (j < 0)
                              + (k >= nz_) + (k < 0) == 1;
        if (i < 0)      v += face_correction(comp, 0, -1, i, j, k, tang);
        if (i >= nx_)   v += face_correction(comp, 0, +1, i, j, k, tang);
        if (j < 0)      v += face_correction(comp, 1, -1, i, j, k, tang);
        if (j >= ny_)   v += face_correction(comp, 1, +1, i, j, k, tang);
        if (k < 0)      v += face_correction(comp, 2, -1, i, j, k, tang);
        if (k >= nz_)   v += face_correction(comp, 2, +1, i, j, k, tang);
        return v;
    }

    /// Discrete equilibrium LHS at node (i,j,k) for component `comp`; the
    /// equation solved (SOR) is L == 0.  Central second derivatives, 4-point
    /// mixed derivatives, ghost-aware neighbor reads, plus body force and the
    /// thermal-gradient pseudo-force.
    double equilibrium(int comp, int i, int j, int k) const
    {
        double out = 0.0;
        for (int t = 0; t < 3; ++t)
        {
            const double coef = (t == comp) ? (lam_ + 2.0 * mu_) : mu_;
            out += coef * d2(comp, t, i, j, k);
        }
        for (int t = 0; t < 3; ++t)
        {
            if (t == comp) continue;
            out += (lam_ + mu_) * mixed(t, comp, t, i, j, k);
        }
        out += cfg_.body_force[static_cast<std::size_t>(comp)];
        out -= kappa_ * alpha_ * dT_grad(comp, i, j, k);
        return out;
    }

    /// Magnitude of the diagonal (own-node) coefficient of component `comp`.
    double A_diag(int comp) const
    {
        const double hx = 2.0 / (dx_ * dx_);
        const double hy = 2.0 / (dy_ * dy_);
        const double hz = 2.0 / (dz_ * dz_);
        const double self = lam_ + 2.0 * mu_;
        switch (comp)
        {
            case 0: return self * hx + mu_ * (hy + hz);
            case 1: return self * hy + mu_ * (hx + hz);
            default: return self * hz + mu_ * (hx + hy);
        }
    }

    /// Nodal strain at (i,j,k), ghost-aware central differences.
    Strain nodal_strain(int i, int j, int k) const
    {
        Strain e;
        e.xx = (value(0, i + 1, j, k) - value(0, i - 1, j, k)) / (2.0 * dx_);
        e.yy = (value(1, i, j + 1, k) - value(1, i, j - 1, k)) / (2.0 * dy_);
        e.zz = (value(2, i, j, k + 1) - value(2, i, j, k - 1)) / (2.0 * dz_);
        e.xy = 0.5 * ((value(0, i, j + 1, k) - value(0, i, j - 1, k)) / (2.0 * dy_)
                      + (value(1, i + 1, j, k) - value(1, i - 1, j, k)) / (2.0 * dx_));
        e.xz = 0.5 * ((value(0, i, j, k + 1) - value(0, i, j, k - 1)) / (2.0 * dz_)
                      + (value(2, i + 1, j, k) - value(2, i - 1, j, k)) / (2.0 * dx_));
        e.yz = 0.5 * ((value(1, i, j, k + 1) - value(1, i, j, k - 1)) / (2.0 * dz_)
                      + (value(2, i, j + 1, k) - value(2, i, j - 1, k)) / (2.0 * dy_));
        return e;
    }

    /// Delta T at (i,j,k), clamped to the grid (zero when no temperature data).
    double dT_at(int i, int j, int k) const
    {
        if (!have_temp_) return 0.0;
        const int ci = clamp_int(i, nx_);
        const int cj = clamp_int(j, ny_);
        const int ck = clamp_int(k, nz_);
        return dT_[node_idx(ci, cj, ck, nx_, ny_)];
    }

private:
    std::size_t idx(int i, int j, int k) const
    {
        return node_idx(i, j, k, nx_, ny_);
    }

    bool in(int i, int j, int k) const
    {
        return i >= 0 && i < nx_ && j >= 0 && j < ny_ && k >= 0 && k < nz_;
    }

    int reflect(int p, int n) const
    {
        if (p < 0) return -p;
        if (p >= n) return 2 * (n - 1) - p;
        return p;
    }

    double axis_h(int ax) const
    {
        return ax == 0 ? dx_ : (ax == 1 ? dy_ : dz_);
    }

    double raw(int comp, int i, int j, int k) const
    {
        return u_[static_cast<std::size_t>(comp) * N_ + idx(i, j, k)];
    }

    /// Traction vector on the face (ax, sgn) at the boundary node (bi,bj,bk).
    /// Only the +z face carries config.surface_traction (optionally masked).
    std::array<double, 3> face_traction(int ax, int sgn, int bi, int bj, int bk) const
    {
        std::array<double, 3> t = {0.0, 0.0, 0.0};
        if (ax == 2 && sgn > 0)
        {
            const std::size_t ni = idx(bi, bj, bk);
            const bool masked = cfg_.traction_mask.empty() || cfg_.traction_mask[ni];
            if (masked) t = cfg_.surface_traction;
        }
        return t;
    }

    /// Derivative of component `comp` along axis `ax` at (i,j,k), using only
    /// in-domain neighbors (one-sided at the grid boundary; both central and
    /// one-sided are exact for linear fields, so patch tests are unaffected).
    double deriv_in_domain(int comp, int ax, int i, int j, int k) const
    {
        const double h = axis_h(ax);
        int p1[3] = {i, j, k};
        int p0[3] = {i, j, k};
        p1[ax] += 1;
        p0[ax] -= 1;
        const bool b1 = in(p1[0], p1[1], p1[2]);
        const bool b0 = in(p0[0], p0[1], p0[2]);
        if (b0 && b1)
            return (raw(comp, p1[0], p1[1], p1[2]) - raw(comp, p0[0], p0[1], p0[2]))
                   / (2.0 * h);
        if (b0)
            return (raw(comp, i, j, k) - raw(comp, p0[0], p0[1], p0[2])) / h;
        if (b1)
            return (raw(comp, p1[0], p1[1], p1[2]) - raw(comp, i, j, k)) / h;
        return 0.0;
    }

    /// Second derivative of component `comp` along axis `ax` (ghost-aware).
    double d2(int comp, int ax, int i, int j, int k) const
    {
        const double h = axis_h(ax);
        int p1[3] = {i, j, k};
        int p0[3] = {i, j, k};
        p1[ax] += 1;
        p0[ax] -= 1;
        return (value(comp, p1[0], p1[1], p1[2])
                - 2.0 * value(comp, i, j, k)
                + value(comp, p0[0], p0[1], p0[2])) / (h * h);
    }

    /// 4-point mixed derivative d2/d(ax_a)d(ax_b) of component `comp`.
    double mixed(int comp, int ax_a, int ax_b, int i, int j, int k) const
    {
        const double h_a = axis_h(ax_a);
        const double h_b = axis_h(ax_b);
        int p[3] = {i, j, k};
        p[ax_a] += 1; p[ax_b] += 1;
        const double vpp = value(comp, p[0], p[1], p[2]);
        p[ax_a] -= 2;
        const double vmp = value(comp, p[0], p[1], p[2]);
        p[ax_a] += 2; p[ax_b] -= 2;
        const double vpm = value(comp, p[0], p[1], p[2]);
        p[ax_a] -= 2;
        const double vmm = value(comp, p[0], p[1], p[2]);
        return ((vpp - vmp) - (vpm - vmm)) / (4.0 * h_a * h_b);
    }

    /// Derivative of the temperature field along axis `ax` (zero without data).
    double dT_grad(int ax, int i, int j, int k) const
    {
        if (!have_temp_) return 0.0;
        const double h = axis_h(ax);
        int p1[3] = {i, j, k};
        int p0[3] = {i, j, k};
        p1[ax] += 1;
        p0[ax] -= 1;
        const bool b1 = in(p1[0], p1[1], p1[2]);
        const bool b0 = in(p0[0], p0[1], p0[2]);
        if (b0 && b1) return (dT_at(p1[0], p1[1], p1[2]) - dT_at(p0[0], p0[1], p0[2])) / (2.0 * h);
        if (b0) return (dT_at(i, j, k) - dT_at(p0[0], p0[1], p0[2])) / h;
        if (b1) return (dT_at(p1[0], p1[1], p1[2]) - dT_at(i, j, k)) / h;
        return 0.0;
    }

    /// Traction-Robin correction for a ghost value of component `comp` caused
    /// by the single face (ax, sgn) adjacent to the ghost point (gi,gj,gk).
    ///   normal  comp: u[g] = mirror + 2*h/(lam+2mu)*(t_n + kappa*alpha*dT
    ///                                                    - lam*(eps_bb+eps_cc))
    ///   tangent comp: u[g] = mirror + 2*h/mu*(t_t - d(u_normal)/d(tangent))
    /// Derivative samples use only in-domain neighbors, keeping edge/corner
    /// ghost evaluation non-recursive.
    double face_correction(int comp, int ax, int sgn, int gi, int gj, int gk,
                               bool tangential_ok) const
    {
        const int bi = clamp_int(gi, nx_);
        const int bj = clamp_int(gj, ny_);
        const int bk = clamp_int(gk, nz_);
        const double h = axis_h(ax);
        const std::array<double, 3> t = face_traction(ax, sgn, bi, bj, bk);

        if (comp == ax)
        {
            const int t1 = (ax + 1) % 3;
            const int t2 = (ax + 2) % 3;
            const double e1 = deriv_in_domain(t1, t1, bi, bj, bk);
            const double e2 = deriv_in_domain(t2, t2, bi, bj, bk);
            const double th = kappa_ * alpha_ * dT_at(bi, bj, bk);
            return (2.0 * h / (lam_ + 2.0 * mu_))
                   * (t[static_cast<std::size_t>(ax)] + th - lam_ * (e1 + e2));
        }
        // tangential component: sigma_{normal,tangent} = t_tangent:
        //   (u_c[g] - u_c[mirror])/(2h) = t_tangent/mu - d(u_normal)/d(x_tangent)
        if (!tangential_ok) return 0.0; // edge/corner ghost: mirror only
        const double du_normal_d_t = deriv_in_domain(ax, comp, bi, bj, bk);
        return (2.0 * h / mu_) * t[static_cast<std::size_t>(comp)]
               - 2.0 * h * du_normal_d_t;
    }

    const ElasticityConfig& cfg_;
    std::vector<double>& u_;
    const std::vector<double>& dT_;
    int32_t nx_ = 0;
    int32_t ny_ = 0;
    int32_t nz_ = 0;
    double dx_ = 0.0;
    double dy_ = 0.0;
    double dz_ = 0.0;
    double lam_ = 0.0;
    double mu_ = 0.0;
    double kappa_ = 0.0;
    double alpha_ = 0.0;
    bool have_temp_ = false;
    bool mirror_only_ = false;
    std::size_t N_ = 0;
};


/// Mechanical energy of a transient state: KE (lumped nodal) + PE (cell
/// strain energy, same convention as the static result).
struct MechanicalEnergy
{
    double kinetic = 0.0;
    double potential = 0.0;
    double total() const { return kinetic + potential; }
};

MechanicalEnergy mechanical_energy(const Solver& solver, const std::vector<double>& u,
                                   const std::vector<double>& v, double rho,
                                   std::size_t N, const ElasticityConfig& config)
{
    const int nx = config.grid.dims[0], ny = config.grid.dims[1], nz = config.grid.dims[2];
    const double cell_vol = config.grid.spacing[0] * config.grid.spacing[1]
                            * config.grid.spacing[2];
    double ke = 0.0;
    for (std::size_t n = 0; n < N; ++n)
    {
        const double vx = v[0 * N + n], vy = v[1 * N + n], vz = v[2 * N + n];
        ke += 0.5 * rho * (vx * vx + vy * vy + vz * vz) * cell_vol;
    }
    std::vector<Strain> nodal(N);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                nodal[node_idx(i, j, k, nx, ny)] = solver.nodal_strain(i, j, k);
    const double lam = solver.lam(), mu = solver.mu();
    double pe = 0.0;
    for (int ck = 0; ck < nz - 1; ++ck)
        for (int cj = 0; cj < ny - 1; ++cj)
            for (int ci = 0; ci < nx - 1; ++ci)
            {
                Strain e = {};
                for (int dk = 0; dk <= 1; ++dk)
                    for (int dj = 0; dj <= 1; ++dj)
                        for (int di = 0; di <= 1; ++di)
                        {
                            const Strain& c = nodal[node_idx(ci + di, cj + dj, ck + dk, nx, ny)];
                            e.xx += c.xx; e.yy += c.yy; e.zz += c.zz;
                            e.xy += c.xy; e.xz += c.xz; e.yz += c.yz;
                        }
                e.xx /= 8.0; e.yy /= 8.0; e.zz /= 8.0;
                e.xy /= 8.0; e.xz /= 8.0; e.yz /= 8.0;
                const double tr = e.xx + e.yy + e.zz;
                const double sxxe = lam * tr + 2.0 * mu * e.xx;
                const double syye = lam * tr + 2.0 * mu * e.yy;
                const double szze = lam * tr + 2.0 * mu * e.zz;
                const double sxy = 2.0 * mu * e.xy;
                const double sxz = 2.0 * mu * e.xz;
                const double syz = 2.0 * mu * e.yz;
                pe += 0.5 * (sxxe * e.xx + syye * e.yy + szze * e.zz
                             + 2.0 * (sxy * e.xy + sxz * e.xz + syz * e.yz)) * cell_vol;
            }
    (void)config.thermal_expansion_coefficient;
    return MechanicalEnergy{ke, pe};
}

} // namespace

bool validate_elasticity_config(const ElasticityConfig& config,
                                std::string& error,
                                std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    const auto& dims = config.grid.dims;
    if (dims[0] < 2 || dims[1] < 2 || dims[2] < 2)
    {
        error = "elasticity: grid dims must be >= 2 per axis";
        return false;
    }
    const auto& s = config.grid.spacing;
    if (!(s[0] > 0.0) || !(s[1] > 0.0) || !(s[2] > 0.0))
    {
        error = "elasticity: grid spacing must be > 0";
        return false;
    }
    if (!(config.material.elastic_modulus > 0.0))
    {
        error = "elasticity: elastic_modulus must be > 0";
        return false;
    }
    const double nu = config.material.poisson_ratio;
    if (!(nu > -1.0 && nu < 0.5))
    {
        error = "elasticity: poisson_ratio must be in (-1, 0.5)";
        return false;
    }
    if (!(config.tolerance > 0.0))
    {
        error = "elasticity: tolerance must be > 0";
        return false;
    }
    if (config.max_iterations == 0)
    {
        error = "elasticity: max_iterations must be > 0";
        return false;
    }
    if (config.thermal_expansion_coefficient < 0.0)
    {
        error = "elasticity: thermal_expansion_coefficient must be >= 0";
        return false;
    }
    for (int c = 0; c < 3; ++c)
    {
        if (!std::isfinite(config.body_force[static_cast<std::size_t>(c)]))
        {
            error = "elasticity: body_force must be finite";
            return false;
        }
        if (!std::isfinite(config.surface_traction[static_cast<std::size_t>(c)]))
        {
            error = "elasticity: surface_traction must be finite";
            return false;
        }
    }

    const std::size_t n = static_cast<std::size_t>(dims[0])
                          * static_cast<std::size_t>(dims[1])
                          * static_cast<std::size_t>(dims[2]);
    if (!config.fixed_mask.empty() && config.fixed_mask.size() != n)
    {
        error = "elasticity: fixed_mask size must match the node count";
        return false;
    }
    if (!config.fixed_displacement.empty() && config.fixed_displacement.size() != n)
    {
        error = "elasticity: fixed_displacement size must match the node count";
        return false;
    }
    if (!config.traction_mask.empty() && config.traction_mask.size() != n)
    {
        error = "elasticity: traction_mask size must match the node count";
        return false;
    }

    if (config.fixed_mask.empty() && !config.fixed_displacement.empty())
        warnings.push_back("elasticity: fixed_displacement supplied but fixed_mask is empty; ignored");
    if (!config.fixed_mask.empty() && config.fixed_displacement.empty())
        warnings.push_back("elasticity: no fixed_displacement supplied; pinned nodes use zero displacement");
    return true;
}

ElasticityResult solve_elasticity(const ElasticityConfig& config,
                                  ModelStatus& status,
                                  const exd::engine::coupling::IScalarField3D* temperature_channel)
{
    ElasticityResult res;
    std::string error;
    std::vector<std::string> warnings;
    if (!validate_elasticity_config(config, error, warnings))
    {
        res.status = ModelStatus{false, error, warnings};
        status = res.status;
        return res;
    }

    const int32_t nx = config.grid.dims[0];
    const int32_t ny = config.grid.dims[1];
    const int32_t nz = config.grid.dims[2];
    const std::size_t N = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny)
                          * static_cast<std::size_t>(nz);
    const std::size_t stride = 3 * N;

    // ── optional temperature channel: delta T vs the grid-origin node ──
    const double alpha = config.thermal_expansion_coefficient;
    std::vector<double> dT(N, 0.0);
    bool temp_oob = false;
    if (temperature_channel != nullptr && alpha > 0.0)
    {
        double t_ref = 0.0;
        temperature_channel->sample(config.grid.origin, t_ref);
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const std::array<double, 3> p = {
                        config.grid.origin[0] + static_cast<double>(i) * config.grid.spacing[0],
                        config.grid.origin[1] + static_cast<double>(j) * config.grid.spacing[1],
                        config.grid.origin[2] + static_cast<double>(k) * config.grid.spacing[2]};
                    const std::size_t n = node_idx(i, j, k, nx, ny);
                    double t = 0.0;
                    if (temperature_channel->sample(p, t)) dT[n] = t - t_ref;
                    else temp_oob = true;
                }
        if (temp_oob)
            warnings.push_back("elasticity: temperature channel went out of bounds; thermal strain set to 0 there");
    }
    else if (temperature_channel != nullptr && alpha == 0.0)
    {
        warnings.push_back("elasticity: temperature channel supplied but thermal_expansion_coefficient == 0; ignored");
    }
    else if (alpha > 0.0)
    {
        warnings.push_back("elasticity: thermal_expansion_coefficient > 0 with no temperature channel; thermal strain ignored");
    }

    // ── state + Dirichlet initialization ──
    std::vector<double> u(stride, 0.0);
    std::vector<std::array<bool, 3>> pinned(N, {false, false, false});
    if (!config.fixed_mask.empty()) pinned = config.fixed_mask;
    if (!config.fixed_mask.empty())
    {
        for (std::size_t n = 0; n < N; ++n)
            for (int c = 0; c < 3; ++c)
            {
                if (!pinned[n][static_cast<std::size_t>(c)]) continue;
                const double v = config.fixed_displacement.empty()
                                     ? 0.0
                                     : config.fixed_displacement[n][static_cast<std::size_t>(c)];
                u[static_cast<std::size_t>(c) * N + n] = v;
            }
    }

    Solver solver(config, u, dT);

    // ── transient path: rho·d²u/dt² = div σ + rho·f (velocity-Verlet) ──
    if (config.transient)
    {
        const double rho = config.density;
        const double h_min = std::min({config.grid.spacing[0],
                                       config.grid.spacing[1],
                                       config.grid.spacing[2]});
        const double cp = p_wave_speed(config.material, rho);
        double dt = config.dt;
        if (dt <= 0.0) dt = 0.3 * h_min / cp;
        if (dt > 0.3 * h_min / cp)
        {
            dt = 0.3 * h_min / cp;
            warnings.push_back("elasticity: dt clamped to the CFL limit 0.3·h/c_p");
        }
        if (!config.initial_displacement.empty())
        {
            if (config.initial_displacement.size() != stride)
            {
                res.status = ModelStatus{false, "elasticity: initial_displacement must be 3·N", warnings};
                status = res.status;
                return res;
            }
            u = config.initial_displacement;
        }
        std::vector<double> v(stride, 0.0);
        if (!config.initial_velocity.empty())
        {
            if (config.initial_velocity.size() != stride)
            {
                res.status = ModelStatus{false, "elasticity: initial_velocity must be 3·N", warnings};
                status = res.status;
                return res;
            }
            v = config.initial_velocity;
        }
        for (std::size_t n = 0; n < N; ++n)
            for (int c = 0; c < 3; ++c)
                if (pinned[n][static_cast<std::size_t>(c)])
                {
                    const double pv = config.fixed_displacement.empty()
                                          ? 0.0 : config.fixed_displacement[n][static_cast<std::size_t>(c)];
                    u[static_cast<std::size_t>(c) * N + n] = pv;
                    v[static_cast<std::size_t>(c) * N + n] = 0.0;
                }

        // finalize the solver for the transient state (u/dT already set)
        Solver tsolver(config, u, dT, /*mirror_only=*/true);
        std::vector<double> a(stride, 0.0);
        double max_disp = 0.0;
        std::vector<double> probe;
        std::vector<double> energy_history;   // total mechanical energy per step
        if (config.probe_index >= 0)
            probe.push_back(u[0 * N + static_cast<std::size_t>(config.probe_index)]);
        for (uint64_t step = 0; step < config.max_steps; ++step)
        {
            // energy each step (symplectic integrators conserve a MODIFIED
            // energy; the flight total is the conserved quantity)
            energy_history.push_back(
                mechanical_energy(tsolver, u, v, rho, N, config).total());

            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        const std::size_t n = node_idx(i, j, k, nx, ny);
                        for (int comp = 0; comp < 3; ++comp)
                        {
                            if (pinned[n][static_cast<std::size_t>(comp)])
                            {
                                a[static_cast<std::size_t>(comp) * N + n] = 0.0;
                                continue;
                            }
                            const double L = tsolver.equilibrium(comp, i, j, k);
                            // L = div σ + f (specific-force convention); the
                            // equation of motion is rho·a = L0 + rho·f.
                            const double L0 = L - config.body_force[static_cast<std::size_t>(comp)];
                            a[static_cast<std::size_t>(comp) * N + n]
                                = L0 / rho + config.body_force[static_cast<std::size_t>(comp)];
                            if (!std::isfinite(a[static_cast<std::size_t>(comp) * N + n]))
                            {
                                res.status = ModelStatus{false,
                                    "elasticity: non-finite acceleration during transient", warnings};
                                status = res.status;
                                return res;
                            }
                        }
                    }
            // velocity-Verlet (symplectic): v += dt·a(u); u += dt·v
            for (std::size_t n = 0; n < N; ++n)
            {
                for (int c = 0; c < 3; ++c)
                {
                    if (pinned[n][static_cast<std::size_t>(c)]) continue;
                    const std::size_t o = static_cast<std::size_t>(c) * N + n;
                    v[o] += dt * a[o];
                    u[o] += dt * v[o];
                    max_disp = std::max(max_disp, std::fabs(u[o]));
                }
            }
            if (config.probe_index >= 0)
                probe.push_back(u[0 * N + static_cast<std::size_t>(config.probe_index)]);
        }
        {
            // final energy + flight drift (after the loop): the symplectic
            // modified energy is conserved during flight; the drift is the
            // relative spread of the total over the run (step 0 excluded —
            // the raw initial PE is not the conserved quantity).
            Solver fsolver(config, u, dT, /*mirror_only=*/true);
            const auto e1 = mechanical_energy(fsolver, u, v, rho, N, config);
            res.kinetic_energy = e1.kinetic;
            res.strain_energy = e1.potential;
            // drift over the settled flight: the first ~10% of steps cover
            // the D'Alembert split transient, whose measured totals swing
            // while the sharp IC devolves; the flight total is conserved.
            double emin = e1.total(), emax = e1.total(), esum = 0.0;
            const size_t n_e = energy_history.size();
            const size_t begin = std::max<size_t>(1, n_e / 10);
            size_t count = 0;
            for (size_t i = begin; i < n_e; ++i)
            {
                emin = std::min(emin, energy_history[i]);
                emax = std::max(emax, energy_history[i]);
                esum += energy_history[i];
                ++count;
            }
            const double emean = count > 0 ? esum / static_cast<double>(count) : e1.total();
            res.energy_drift = (emax - emin) / std::max(emean, 1e-30);
        }
        res.steps = config.max_steps;
        res.dt_used = dt;
        res.time_elapsed = dt * static_cast<double>(config.max_steps);
        res.max_displacement = max_disp;
        res.probe_history = std::move(probe);
        // final state output (same layout as the static path)
        exd::engine::coupling::StructuredVectorGrid disp;
        disp.origin = config.grid.origin;
        disp.spacing = config.grid.spacing;
        disp.dims = config.grid.dims;
        disp.values.assign(stride, 0.0);
        for (std::size_t n = 0; n < N; ++n)
            for (int c = 0; c < 3; ++c)
                disp.values[3 * n + static_cast<std::size_t>(c)]
                    = u[static_cast<std::size_t>(c) * N + n];
        res.displacement = std::move(disp);
        exd::engine::coupling::StructuredVectorGrid vel;
        vel.origin = config.grid.origin;
        vel.spacing = config.grid.spacing;
        vel.dims = config.grid.dims;
        vel.values.assign(stride, 0.0);
        for (std::size_t n = 0; n < N; ++n)
            for (int c = 0; c < 3; ++c)
                vel.values[3 * n + static_cast<std::size_t>(c)]
                    = v[static_cast<std::size_t>(c) * N + n];
        res.velocity = std::move(vel);
        res.ok = true;
        res.status = ModelStatus{true, "", warnings};
        status = res.status;
        return res;
    }

    // The convergence metric is the residual of the discrete equilibrium
    // relative to its initial (load-driven) scale.  The scale is measured on
    // the initial state so the slow global (bending) modes are not masked by
    // tiny displacement magnitudes.
    double eq_scale = 0.0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int comp = 0; comp < 3; ++comp)
                {
                    const std::size_t n = node_idx(i, j, k, nx, ny);
                    if (pinned[n][static_cast<std::size_t>(comp)]) continue;
                    const double e = solver.equilibrium(comp, i, j, k);
                    if (std::fabs(e) > eq_scale) eq_scale = std::fabs(e);
                }
    if (eq_scale == 0.0) eq_scale = 1.0;

    // ── component-wise SOR sweep ──
    double rel_res = 1.0;
    double best_res = 1.0;
    uint64_t best_iter = 0;
    uint64_t iter = 0;
    bool converged = false;
    for (uint64_t it = 0; it < config.max_iterations; ++it)
    {
        double max_corr = 0.0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    for (int comp = 0; comp < 3; ++comp)
                    {
                        const std::size_t n = node_idx(i, j, k, nx, ny);
                        if (pinned[n][static_cast<std::size_t>(comp)]) continue;
                        const double e = solver.equilibrium(comp, i, j, k);
                        if (!std::isfinite(e))
                        {
                            res.status = ModelStatus{
                                false,
                                "elasticity: non-finite residual during SOR",
                                warnings};
                            status = res.status;
                            return res;
                        }
                        // SOR update: eq = -A_diag*u_p + (neighbor terms + f),
                        // so the value that zeroes eq is u_p + eq/A_diag.
                        const double corr = e / solver.A_diag(comp);
                        const double nv = u[static_cast<std::size_t>(comp) * N + n]
                                          + kSorOmega * corr;
                        u[static_cast<std::size_t>(comp) * N + n] = nv;
                        if (std::fabs(e) > max_corr) max_corr = std::fabs(e);
                    }
                }
        iter = it + 1;
        // Residual of the current iterate relative to the initial load scale.
        const double res_now = max_corr / eq_scale;
        // Report the best (smallest) relative residual achieved.
        if (res_now < rel_res) rel_res = res_now;
        if (res_now < best_res)
        {
            best_res = res_now;
            best_iter = it;
        }
        if (res_now < config.tolerance)
        {
            converged = true;
            break;
        }
        // Stagnation guard: if the relative residual has not improved by
        // more than 1e-6 over the last 5000 sweeps, further SOR sweeps will
        // not reach the tolerance (poor conditioning); report and stop.
        if (it > 5000 && best_res > 0.0 && best_res / rel_res > 1.0 - 1e-6
            && it - best_iter > 5000)
        {
            warnings.push_back(
                "elasticity: SOR stagnation (best relative residual "
                + std::to_string(best_res) + "); stopping early");
            break;
        }
    }
    if (!converged)
    {
        warnings.push_back("elasticity: SOR reached max_iterations without meeting tolerance "
                           "(relative residual " + std::to_string(rel_res) + ")");
    }

    // ── output: displacement grid, strain energy, effective strain ──
    exd::engine::coupling::StructuredVectorGrid disp;
    disp.origin = config.grid.origin;
    disp.spacing = config.grid.spacing;
    disp.dims = config.grid.dims;
    disp.values.assign(stride, 0.0);
    for (std::size_t n = 0; n < N; ++n)
        for (int c = 0; c < 3; ++c)
            disp.values[3 * n + static_cast<std::size_t>(c)]
                = u[static_cast<std::size_t>(c) * N + n];

    std::vector<Strain> nodal(N);
    double max_eff = 0.0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t n = node_idx(i, j, k, nx, ny);
                const Strain e = solver.nodal_strain(i, j, k);
                nodal[n] = e;
                const double tr = e.xx + e.yy + e.zz;
                const double exd = e.xx - tr / 3.0;
                const double eyd = e.yy - tr / 3.0;
                const double ezd = e.zz - tr / 3.0;
                const double eq = std::sqrt((2.0 / 3.0)
                                            * (exd * exd + eyd * eyd + ezd * ezd
                                               + 2.0 * (e.xy * e.xy + e.xz * e.xz + e.yz * e.yz)));
                if (eq > max_eff) max_eff = eq;
            }

    const double lam = solver.lam();
    const double mu = solver.mu();
    const double kappa = solver.kappa();
    const double cell_vol = config.grid.spacing[0] * config.grid.spacing[1]
                            * config.grid.spacing[2];
    double strain_energy = 0.0;
    for (int ck = 0; ck < nz - 1; ++ck)
        for (int cj = 0; cj < ny - 1; ++cj)
            for (int ci = 0; ci < nx - 1; ++ci)
            {
                // Cell strain = average of the 8 corner nodal strains.
                Strain e = {};
                double dT_cell = 0.0;
                double w = 0.0;
                for (int dk = 0; dk <= 1; ++dk)
                    for (int dj = 0; dj <= 1; ++dj)
                        for (int di = 0; di <= 1; ++di)
                        {
                            const std::size_t n = node_idx(ci + di, cj + dj, ck + dk, nx, ny);
                            const Strain& c = nodal[n];
                            e.xx += c.xx; e.yy += c.yy; e.zz += c.zz;
                            e.xy += c.xy; e.xz += c.xz; e.yz += c.yz;
                            dT_cell += solver.dT_at(ci + di, cj + dj, ck + dk);
                            w += 1.0;
                        }
                e.xx /= w; e.yy /= w; e.zz /= w;
                e.xy /= w; e.xz /= w; e.yz /= w;
                dT_cell /= w;
                const double tr = e.xx + e.yy + e.zz;
                const double th_diag = -kappa * alpha * dT_cell;
                const double sxx = lam * tr + 2.0 * mu * e.xx + th_diag;
                const double syy = lam * tr + 2.0 * mu * e.yy + th_diag;
                const double szz = lam * tr + 2.0 * mu * e.zz + th_diag;
                const double sxy = 2.0 * mu * e.xy;
                const double sxz = 2.0 * mu * e.xz;
                const double syz = 2.0 * mu * e.yz;
                const double dot = sxx * e.xx + syy * e.yy + szz * e.zz
                                   + 2.0 * (sxy * e.xy + sxz * e.xz + syz * e.yz);
                strain_energy += 0.5 * dot * cell_vol;
            }

    res.ok = true;
    res.status = ModelStatus{true, "", warnings};
    status = res.status;
    res.displacement = std::move(disp);
    res.max_residual = rel_res;
    res.iterations = iter;
    res.strain_energy = strain_energy;
    res.max_effective_strain = max_eff;
    return res;
}

} // namespace exd::engine::physics::structural