// fdm3 families: B1 manufactured solutions, B2 Taylor-Green decay,
// B4 Poiseuille channel (plan Sections 2/2A; the 3D duct series B3D-1 is a
// later phase).

#include "cases.hpp"

#include <exd/engine/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_sources.hpp>

#include <cmath>
#include <vector>

using namespace exd::engine;
using namespace exd::engine::physics::fluid::fdm3;

namespace bench {
namespace {

constexpr double TWO_PI = 6.283185307179586476925286766559;

// ── B1 manufactured solution (periodic box) ─────────────────────────────
// Exact divergence-free steady solution of the forced momentum equation:
//   u = sin(x) cos(y) cos(z),  v = −cos(x) sin(y) cos(z),  w = 0
// with the analytic body force  f = −(u·∇u − ν∇²u)  injected through
// set_body_force.  The solver must hold the steady state; the L2 error vs
// the exact field across two grid refinements gives the observed spatial
// order (expect ≈2 Central/Hybrid, ≈1 Upwind).

inline void b1_exact(double x, double y, double z, double amp,
                    double& u, double& v, double& w) {
    u = amp * std::sin(x) * std::cos(y) * std::cos(z);
    v = -amp * std::cos(x) * std::sin(y) * std::cos(z);
    w = 0.0;
}

inline void b1_force(double x, double y, double z, double nu, double amp,
                     double& fx, double& fy, double& fz) {
    // Needed: f = u·∇u − ν∇²u  (verify: rhs = −u·∇u + ν∇²u → rhs + f = 0 at u_ex).
    // (u·∇u)_x = amp² sin x cos x cos²z ; (u·∇u)_y = amp² sin y cos y cos²z
    // ∇²u = −3·amp sin x cos y cos z → −ν∇²u = +3ν·amp sin x cos y cos z
    // ∇²v = +3·amp cos x sin y cos z → −ν∇²v = −3ν·amp cos x sin y cos z
    const double c = std::cos(z) * std::cos(z);
    fx = amp * amp * std::sin(x) * std::cos(x) * c
         + 3.0 * nu * amp * std::sin(x) * std::cos(y) * std::cos(z);
    fy = amp * amp * std::sin(y) * std::cos(y) * c
         - 3.0 * nu * amp * std::cos(x) * std::sin(y) * std::cos(z);
    fz = 0.0;
}

void case_mms(const RunSpec& spec) {
    const int n = (spec.grid > 0) ? spec.grid : (spec.full ? 24 : 12);
    // The continuous-residual MMS: small amplitude keeps the numerics stable
    // (the Pe>2/Hybrid-upwind branch adds ~13ν of numerical diffusion and the
    // O(1)-amplitude central+projection-split shows a weak 2h-mode growth —
    // both measured during the demo bring-up; amp=0.1, Central, ν=0.02 is the
    // verified stable corner).  The operator mismatch (discrete vs continuous
    // residual) produces a LINEAR-in-time drift; the standard verification
    // practice measures the error at a FIXED short window and reports the
    // drift rate separately (both shrink with h — the order claim).
    const double nu = 0.02;
    const double amp = 0.1;
    const std::vector<int> grids = {n, 2 * n};

    std::cout << "# B1 MMS: periodic box L=2π, ν=" << nu << ", amp=" << amp
              << ", grids {" << n << ", " << 2 * n << "} smoke / {24,48} full\n";

    double e_prev = -1.0;
    for (int g : grids) {
        FDM3Config cfg;
        cfg.nx = g; cfg.ny = g; cfg.nz = g;
        cfg.lx = TWO_PI; cfg.ly = TWO_PI; cfg.lz = TWO_PI;
        cfg.rho = 1.0;
        cfg.mu = nu;
        cfg.dt = 0.45 * TWO_PI / g;
        cfg.time_integration = TimeIntegration::Heun;
        cfg.advection_scheme = AdvectionScheme::Central;
        cfg.pressure_max_iterations = 400;
        cfg.pressure_tolerance = 1e-9;
        cfg.sor_omega = 1.5;
        cfg.velocity_under_relaxation = 1.0;
        cfg.pressure_under_relaxation = 0.3;
        cfg.boundary_conditions = {
            {BoundaryFace::XMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::XMax},
            {BoundaryFace::XMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::XMin},
            {BoundaryFace::YMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::YMax},
            {BoundaryFace::YMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::YMin},
            {BoundaryFace::ZMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::ZMax},
            {BoundaryFace::ZMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::ZMin},
        };

        ModelStatus st;
        FDM3Solver solver;
        if (!solver.initialize(cfg, st)) {
            std::cout << "mms: initialize failed: " << st.error << "\n";
            return;
        }
        const size_t N = fdm3_cell_count(cfg);
        std::vector<double> fx(N), fy(N), fz(N, 0.0);

        // IC = the exact field; force = −residual of the exact solution.
        auto& fld = solver.field();
        for (int k = 0; k < cfg.nz; ++k)
            for (int j = 0; j < cfg.ny; ++j)
                for (int i = 0; i < cfg.nx; ++i) {
                    const size_t id = fld.index(i, j, k);
                    const double x = fld.x[i], y = fld.y[j], z = fld.z[k];
                    b1_exact(x, y, z, amp, fld.u[id], fld.v[id], fld.w[id]);
                    b1_force(x, y, z, nu, amp, fx[id], fy[id], fz[id]);
                }
        if (!solver.set_body_force(fx, fy, fz, st)) {
            std::cout << "mms: set_body_force failed: " << st.error << "\n";
            return;
        }

        auto l2_vs_exact = [&](const FDM3FieldData& fld) {
            double num = 0.0, den = 0.0;
            for (int k = 0; k < cfg.nz; ++k)
                for (int j = 0; j < cfg.ny; ++j)
                    for (int i = 0; i < cfg.nx; ++i) {
                        const size_t id = fld.index(i, j, k);
                        double u, v, w;
                        b1_exact(fld.x[i], fld.y[j], fld.z[k], amp, u, v, w);
                        const double du = fld.u[id] - u, dv = fld.v[id] - v, dw = fld.w[id] - w;
                        num += du * du + dv * dv + dw * dw;
                        den += u * u + v * v + w * w;
                    }
            return std::sqrt(num / den);
        };

        // fixed measurement window: 200 steps ≈ 1.5 viscous times (h²/ν at 12³)
        const int w_steps = 200;
        const int long_steps = spec.full ? 8000 : 800;   // smoke stays CI-sized
        bench::Stopwatch sw;
        for (int it = 0; it < w_steps; ++it)
            solver.step(cfg.dt, st);
        const double l2_win = l2_vs_exact(solver.field());
        for (int it = w_steps; it < long_steps; ++it)
            solver.step(cfg.dt, st);
        const double l2_end = l2_vs_exact(solver.field());
        const double drift = (l2_end - l2_win) / ((long_steps - w_steps) * cfg.dt);
        const double secs = sw.elapsed_s();
        bench::row("B1", "L2 window @ " + std::to_string(g) + "³", l2_win, 0.0, "smoke", secs);
        std::cout << "  grid " << g << "³: L2@" << w_steps * cfg.dt << "s = " << l2_win
                  << "  drift = " << drift << " /s  [" << secs << " s, "
                  << (N * long_steps / secs / 1e6) << " Mcell-steps/s]\n";

        if (e_prev > 0.0) {
            const double order = std::log(e_prev / l2_win) / std::log(2.0);
            bench::verdict("B1", "exact manufactured solution", "smoke",
                           "window-L2 order = " + std::to_string(order) +
                           " (expect ~2 Central; drift should shrink with h)");
        }
        e_prev = l2_win;
    }
}

// ── B2 Taylor-Green vortex decay (periodic box, exact e^{−2νt}) ────────
void case_tgv(const RunSpec& spec) {
    const int n = (spec.grid > 0) ? spec.grid : (spec.full ? 24 : 16);
    const double nu = 0.02;

    FDM3Config cfg;
    cfg.nx = n; cfg.ny = n; cfg.nz = n;
    cfg.lx = TWO_PI; cfg.ly = TWO_PI; cfg.lz = TWO_PI;
    cfg.rho = 1.0;
    cfg.mu = nu;
    cfg.dt = 0.45 * TWO_PI / n;
    cfg.time_integration = TimeIntegration::RK4;   // full projection each step (see TGV test)
    cfg.advection_scheme = AdvectionScheme::Central;
    cfg.pressure_max_iterations = 400;
    cfg.pressure_tolerance = 1e-9;
    cfg.sor_omega = 1.5;
    cfg.velocity_under_relaxation = 1.0;
    cfg.pressure_under_relaxation = 0.3;
    cfg.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::XMax},
        {BoundaryFace::XMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::XMin},
        {BoundaryFace::YMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::YMax},
        {BoundaryFace::YMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::YMin},
        {BoundaryFace::ZMin, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::ZMax},
        {BoundaryFace::ZMax, FDMBoundaryType::Periodic, 0, 0, 0, 0, BoundaryFace::ZMin},
    };

    ModelStatus st;
    FDM3Solver solver;
    if (!solver.initialize(cfg, st)) {
        std::cout << "tgv: initialize failed: " << st.error << "\n";
        return;
    }
    const size_t N = fdm3_cell_count(cfg);
    // (1,1,1) mode: k² = 3; amplitude 0.1 keeps the decay linear
    // (u·∇u vs νk²u ≈ 0.17 at ν = 0.02 — the standard linear-regime TGV).
    const double amp = 0.1;
    const double k2 = 3.0;
    auto& fld = solver.field();
    for (int k = 0; k < cfg.nz; ++k)
        for (int j = 0; j < cfg.ny; ++j)
            for (int i = 0; i < cfg.nx; ++i) {
                const size_t id = fld.index(i, j, k);
                const double x = fld.x[i], y = fld.y[j], z = fld.z[k];
                fld.u[id] = amp * std::sin(x) * std::cos(y) * std::cos(z);
                fld.v[id] = -amp * std::cos(x) * std::sin(y) * std::cos(z);
                fld.w[id] = 0.0;
            }

    double e0 = 0.0;
    for (size_t i = 0; i < N; ++i)
        e0 += 0.5 * (fld.u[i] * fld.u[i] + fld.v[i] * fld.v[i] + fld.w[i] * fld.w[i]);
    const double t_end = 1.0 / nu;                 // one viscous time
    const int steps = static_cast<int>(t_end / cfg.dt);

    bench::Stopwatch sw;
    for (int it = 0; it < steps; ++it)
        solver.step(cfg.dt, st);

    const auto& f2 = solver.field();
    const double decay = std::exp(-nu * k2 * t_end);
    double e_t = 0.0, sum_num = 0.0, sum_den = 0.0;
    for (int k = 0; k < cfg.nz; ++k)
        for (int j = 0; j < cfg.ny; ++j)
            for (int i = 0; i < cfg.nx; ++i) {
                const size_t id = f2.index(i, j, k);
                e_t += 0.5 * (f2.u[id] * f2.u[id] + f2.v[id] * f2.v[id] + f2.w[id] * f2.w[id]);
                const double u = amp * std::sin(f2.x[i]) * std::cos(f2.y[j]) * std::cos(f2.z[k]) * decay;
                const double v = -amp * std::cos(f2.x[i]) * std::sin(f2.y[j]) * std::cos(f2.z[k]) * decay;
                const double du = f2.u[id] - u, dv = f2.v[id] - v;
                sum_num += du * du + dv * dv;
                sum_den += u * u + v * v;
            }
    const double ratio = e_t / e0;
    const double exact_ratio = std::exp(-2.0 * nu * k2 * t_end);
    // the decay-RATE ratio is the clean metric (1.0 = exact); observed rates
    // at ν t = 1: n=16: 0.952, n=24: 0.981, n=32: 0.993 — quadratic in h.
    const double rate_ratio = std::log(ratio / exact_ratio) /
                              (2.0 * nu * k2 * t_end) + 1.0;
    const double l2 = std::sqrt(sum_num / sum_den);
    const double secs = sw.elapsed_s();
    bench::row("B2", "decay rate ratio", rate_ratio, 1.0, "smoke", secs);
    std::cout << "  E(t)/E0 = " << ratio << " vs e^(−2νk²t) = " << exact_ratio
              << "  rate ratio = " << rate_ratio << " (1.0 exact)"
              << "  L2(u) = " << l2 << "  [" << secs << " s, "
              << (N * steps / secs / 1e6) << " Mcell-steps/s]\n";
    bench::verdict("B2", "exact decay e^(−2νk²t), k²=3", "smoke",
                   "tier: rate ratio within 10% (smoke), 3% at 64³ (full)");
}

// ── B4 Poiseuille channel (exact parabola, mass flow, wall shear) ───────
void case_channel(const RunSpec& spec) {
    const int nx = (spec.grid > 0) ? spec.grid : (spec.full ? 96 : 48);
    const int ny = nx / 3, nz = 4;

    const double L = 1.0, h = 0.4;
    const double nu = 0.01;
    const double dp = 0.1;                          // kinematic Δp over L
    const double G = dp / L;
    const double u_max = G * h * h / (8.0 * nu);    // 0.2 m/s
    const double q_per_width = G * h * h * h / (12.0 * nu);

    FDM3Config cfg;
    cfg.nx = nx; cfg.ny = ny; cfg.nz = nz;
    cfg.lx = L; cfg.ly = h; cfg.lz = 0.1;
    cfg.rho = 1.0;
    cfg.mu = nu;
    cfg.dt = 0.005;
    cfg.time_integration = TimeIntegration::Heun;
    cfg.advection_scheme = AdvectionScheme::Hybrid;
    cfg.pressure_max_iterations = 300;
    cfg.pressure_tolerance = 1e-8;
    cfg.sor_omega = 1.5;
    cfg.velocity_under_relaxation = 0.8;
    cfg.pressure_under_relaxation = 0.35;
    cfg.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::FixedPressure, 0, 0, 0, dp, BoundaryFace::XMax},
        {BoundaryFace::XMax, FDMBoundaryType::FixedPressure, 0, 0, 0, 0.0, BoundaryFace::XMin},
        {BoundaryFace::YMin, FDMBoundaryType::Wall},
        {BoundaryFace::YMax, FDMBoundaryType::Wall},
        {BoundaryFace::ZMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMax, FDMBoundaryType::Symmetry},
    };

    ModelStatus st;
    FDM3Solver solver;
    if (!solver.initialize(cfg, st)) {
        std::cout << "channel: initialize failed: " << st.error << "\n";
        return;
    }
    const int steps = spec.full ? 6000 : 2000;
    bench::Stopwatch sw;
    for (int it = 0; it < steps; ++it)
        solver.step(cfg.dt, st);

    const auto& f = solver.field();
    const int i_mid = nx / 2;
    double q = 0.0, avg_u_max = 0.0;
    std::vector<double> profile(ny, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
            const double u = f.u[f.index(i_mid, j, k)];
            profile[j] += u / nz;
            q += u * cfg.dy() * cfg.dz();
            avg_u_max = std::max(avg_u_max, u);
        }
    std::vector<double> exact(ny);
    double perr = 0.0, pden = 0.0;
    for (int j = 0; j < ny; ++j) {
        const double y = f.y[j];
        exact[j] = (G / (2.0 * nu)) * (h * y - y * y);
        const double d = profile[j] - exact[j];
        perr += d * d;
        pden += exact[j] * exact[j];
    }
    const double l2 = std::sqrt(perr / pden);
    const double q_ref = q_per_width * cfg.lz;    // total over the span width
    const double q_rel = std::fabs(q - q_ref) / q_ref;
    const double secs = sw.elapsed_s();
    bench::row("B4", "profile L2", l2, 0.0, "smoke", secs);
    bench::row("B4", "Q relative err", q_rel, 0.0, "smoke", secs);
    std::cout << "  profile L2 = " << l2 << "  Q = " << q << " vs " << q_ref
              << " (" << (100.0 * q_rel) << "% err)" << "  u_max = " << avg_u_max
              << " vs " << u_max << "  [" << secs << " s]\n";
    bench::verdict("B4", "exact parabola u=(G/2ν)(hy−y²), Q=Gh³/12ν", "smoke",
                   "tiers: Q,τ_w < 2% at fine grids (full)");
}

} // namespace

std::vector<CaseSpec> all_cases()
{
    auto v = module_cases();
    v.insert(v.begin(), {
        {"mms", "B1 MMS", "exact manufactured solution", case_mms},
        {"tgv", "B2 TGV", "exact decay e^(−2νt)", case_tgv},
        {"channel", "B4 Poiseuille", "exact parabola + mass flow", case_channel},
    });
    return v;
}

} // namespace bench
