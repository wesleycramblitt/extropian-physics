// plenum_test.cpp
// Lumped Greitzer surge cell (isentropic plenum mapping):
//   - derivative signs at representative states
//   - numerically-fixed operating point is a fixed point (drift-free march)
//   - linearized stability matches time-march on BOTH sides of neutral
//   - sustained surge limit cycle with flow collapse
//   - domain validation (never throws)

#include <exd/physics/fluid/lumped/plenum.hpp>
#include <exd/physics/thermo/eos.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

using namespace exd::physics::fluid::lumped;
using exd::physics::ModelStatus;
using exd::physics::solver::IntegratorConfig;
using exd::physics::solver::IntegrationMethod;

namespace
{

std::unique_ptr<exd::physics::thermo::IEos> make_air()
{
    return exd::physics::thermo::make_ideal_gas(exd::physics::thermo::IdealGasConfig{});
}

// Greitzer cubic compressor characteristic (pressure RISE over ambient, Pa).
double cubic_dp(double mdot, double A, double B, double C)
{
    return A + B * mdot - C * mdot * mdot * mdot;
}

// Sqrt throttle: gain * sqrt(max(0, p - p_ambient)) (kg/s).
double sqrt_throttle(double p, double p_amb, double gain)
{
    const double dp = p - p_amb;
    return gain * std::sqrt(dp > 0.0 ? dp : 0.0);
}

// Bisection on f over [lo, hi]; requires a sign change.
double bisect_positive(const std::function<double(double)>& f, double lo, double hi)
{
    const double f_lo = f(lo);
    const double f_hi = f(hi);
    REQUIRE((f_lo <= 0.0) != (f_hi <= 0.0));
    for (int i = 0; i < 120; ++i)
    {
        const double mid = 0.5 * (lo + hi);
        const double f_mid = f(mid);
        if ((f_lo <= 0.0) == (f_mid <= 0.0)) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// 2x2 Jacobian of the plenum ODE by central finite differences.
struct PlenumJacobian
{
    double j00 = 0.0; // d(dmdot)/d(mdot)
    double j01 = 0.0; // d(dmdot)/d(p)
    double j10 = 0.0; // d(dp)/d(mdot)
    double j11 = 0.0; // d(dp)/d(p)
    double trace() const { return j00 + j11; }
    double det() const { return j00 * j11 - j01 * j10; }
};

PlenumJacobian plenum_jacobian(const PlenumModelConfig& cfg,
                               const std::function<double(double)>& comp,
                               const std::function<double(double)>& thr,
                               double p, double mdot,
                               const exd::physics::thermo::IEos& eos)
{
    ModelStatus st;
    const auto deriv = [&](double pv, double mv)
    {
        PlenumDerivative d = plenum_derivative(cfg, PlenumState{pv, mv}, comp, thr,
                                               eos, st);
        return std::make_pair(d.dmdot_dt, d.dp_dt);
    };
    const double h = 1e-6;
    const auto dm_p = deriv(p, mdot + h);
    const auto dm_m = deriv(p, mdot - h);
    const auto dp_p = deriv(p + h, mdot);
    const auto dp_m = deriv(p - h, mdot);
    PlenumJacobian j;
    j.j00 = (dm_p.first - dm_m.first) / (2.0 * h);
    j.j01 = (dp_p.first - dp_m.first) / (2.0 * h);
    j.j10 = (dm_p.second - dm_m.second) / (2.0 * h);
    j.j11 = (dp_p.second - dp_m.second) / (2.0 * h);
    return j;
}

} // anonymous namespace

TEST_CASE("plenum: derivative signs at representative states")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);
    const PlenumModelConfig cfg; // volume 0.5, duct L/A = 40, ambient air

    const auto comp = [](double m) { return 30000.0 - 4000.0 * m; }; // rising->falling
    const auto thr = [](double p) { return sqrt_throttle(p, 101325.0, 0.02); };

    ModelStatus st;
    // Low duct flow from ambient: the compressor rise exceeds the plenum back
    // pressure, so the duct flow accelerates.
    PlenumDerivative d_flow_up = plenum_derivative(cfg, PlenumState{101325.0, 0.1},
                                                   comp, thr, *eos, st);
    REQUIRE(d_flow_up.ok);
    CHECK(st.ok);
    CHECK(d_flow_up.dmdot_dt > 0.0);

    // Same point: the plenum is below throttle demand, so it fills up.
    PlenumDerivative d_press_up = plenum_derivative(cfg, PlenumState{101325.0, 0.1},
                                                    comp, thr, *eos, st);
    CHECK(d_press_up.dp_dt > 0.0);

    // High back pressure: the plenum back-pressures the duct, flow decelerates.
    PlenumDerivative d_flow_down = plenum_derivative(cfg, PlenumState{151325.0, 0.1},
                                                     comp, thr, *eos, st);
    CHECK(d_flow_down.dmdot_dt < 0.0);

    // Reverse duct flow from ambient: mass leaves the plenum.
    PlenumDerivative d_press_down = plenum_derivative(cfg, PlenumState{101325.0, -0.05},
                                                      comp, thr, *eos, st);
    CHECK(d_press_down.dp_dt < 0.0);

    // The daughter struct reports the acoustic scale (isentropic mapping).
    CHECK(d_flow_up.speed_of_sound_sq > 0.0);
    CHECK(d_flow_up.inertia == doctest::Approx(2.0 / 0.05).epsilon(1e-12));
}

TEST_CASE("plenum: steady operating point is a fixed point of the ODE")
{
    auto eos = make_air();
    const PlenumModelConfig cfg;

    // Linear compressor dp_c = A - B*mdot with the sqrt throttle:
    // at the operating point A - B*mdot = (mdot/gain)^2.
    const double A = 30000.0;
    const double B = 4000.0;
    const double gain = 0.02;
    const auto comp = [&](double m) { return A - B * m; };
    const auto thr = [&](double p) { return sqrt_throttle(p, cfg.p_ambient, gain); };

    // Fixed point: f(mdot) = dp_c(mdot) - (mdot/gain)^2 = 0.
    auto f = [&](double m) { return (A - B * m) - (m / gain) * (m / gain); };
    const double ms = bisect_positive(f, 0.0, A / B);
    const double ps = cfg.p_ambient + (ms / gain) * (ms / gain);
    CHECK(ms > 0.0);
    CHECK(ps > cfg.p_ambient);

    ModelStatus st;
    PlenumDerivative d = plenum_derivative(cfg, PlenumState{ps, ms}, comp, thr,
                                           *eos, st);
    REQUIRE(st.ok);
    REQUIRE(d.ok);

    // |derivatives| < 1e-6 relative to the model's natural rate scales.
    const double inertia = cfg.duct_length / cfg.duct_area;
    const double gamma = eos->gamma();
    const double R = eos->gas_constant();
    const double dmdot_scale = A / inertia;
    const double dp_scale = (gamma * R * cfg.T_ambient / cfg.volume) * ms;
    CHECK(std::fabs(d.dmdot_dt) < 1e-6 * dmdot_scale);
    CHECK(std::fabs(d.dp_dt) < 1e-6 * dp_scale);

    // Time-march a long spin from the exact fixed point with RK4: the drift
    // must be integrator-level (< 0.1% of the operating-point scale).
    IntegratorConfig integ;
    integ.method = IntegrationMethod::RK4;
    PlenumState s{ps, ms};
    for (int i = 0; i < 2000; ++i)
    {
        s = step_plenum(1.0e-4, cfg, s, comp, thr, *eos, integ, st);
        REQUIRE(st.ok);
    }
    CHECK(std::fabs(s.p_plenum - ps) / (ps - cfg.p_ambient) < 1.0e-3);
    CHECK(std::fabs(s.mdot_duct - ms) / ms < 1.0e-3);
}

TEST_CASE("plenum: stable operating point -- Jacobian and time-march agree")
{
    auto eos = make_air();
    const PlenumModelConfig cfg;

    // Negative compressor slope (stable): dp_c = A - B*mdot, B > 0.
    const double A = 30000.0;
    const double B = 4000.0;
    const double gain = 0.02;
    const auto comp = [&](double m) { return A - B * m; };
    const auto thr = [&](double p) { return sqrt_throttle(p, cfg.p_ambient, gain); };
    auto f = [&](double m) { return (A - B * m) - (m / gain) * (m / gain); };
    const double ms = bisect_positive(f, 0.0, A / B);
    const double ps = cfg.p_ambient + (ms / gain) * (ms / gain);

    // Linearization: tr < 0 and det > 0 (both eigenvalues in the left half-plane).
    PlenumJacobian j = plenum_jacobian(cfg, comp, thr, ps, ms, *eos);
    CHECK(j.trace() < 0.0);
    CHECK(j.det() > 0.0);

    // Time-march a small perturbation: it decays.
    ModelStatus st;
    IntegratorConfig integ;
    integ.method = IntegrationMethod::RK4;
    PlenumState s{ps + 500.0, ms + 0.01};
    const double n0 = (500.0 / (ps - cfg.p_ambient)) + (0.01 / ms);
    for (int i = 0; i < 5000; ++i)
    {
        s = step_plenum(1.0e-4, cfg, s, comp, thr, *eos, integ, st);
        REQUIRE(st.ok);
    }
    const double n1 = std::fabs(s.p_plenum - ps) / (ps - cfg.p_ambient)
                      + std::fabs(s.mdot_duct - ms) / ms;
    CHECK(n1 < 0.2 * n0);
}

TEST_CASE("plenum: positive-slope instability -- weighted criterion and growth")
{
    auto eos = make_air();
    const PlenumModelConfig cfg;

    // Classical Greitzer cubic with a positive-slope middle segment.
    //   dp_c(m) = A + B*m - C*m^3  (A,B,C all > 0)
    // The middle fixed point (rising side) is selected with the sqrt throttle.
    const double A = 40000.0;
    const double B = 12000.0;
    const double C = 2500.0;
    const double gain = 0.005;
    const auto comp = [&](double m) { return cubic_dp(m, A, B, C); };
    const auto thr = [&](double p) { return sqrt_throttle(p, cfg.p_ambient, gain); };
    auto f = [&](double m) { return m - sqrt_throttle(cfg.p_ambient + cubic_dp(m, A, B, C),
                                                      cfg.p_ambient, gain); };
    const double ms = bisect_positive(f, 0.5, 1.2649); // rising side (peak at sqrt(B/3C))
    const double ps = cfg.p_ambient + (ms / gain) * (ms / gain);
    const double slope = B - 3.0 * C * ms * ms;
    REQUIRE(slope > 0.0); // the operating point sits on the positive-slope segment

    // NOTE on the criterion: with a positive-pressure-rise compressor and the
    // sqrt throttle the product slope*d(throttle)/dp is ALWAYS <= 1/2 (fixed
    // point algebra: m = g*sqrt(dp_c), so s*C_th = s*m/(2*dp_c) < 1 since
    // 2*dp_c - s*m = 2A + B*m + ... > 0). A saddle (det < 0) is therefore
    // structurally impossible; the repulsion here is the positive TRACE
    // (unstable spiral) -- exactly the weighted criterion the naive slope
    // comparison miss-estimates (see the stable-positive-slope case below).
    PlenumJacobian j = plenum_jacobian(cfg, comp, thr, ps, ms, *eos);
    CHECK(j.trace() > 0.0);

    // Time-march a tiny perturbation: it grows by > 10x before saturation.
    ModelStatus st;
    IntegratorConfig integ;
    integ.method = IntegrationMethod::RK4;
    PlenumState s{ps, ms + 1.0e-4};
    double max_growth = 0.0;
    for (int i = 0; i < 15000; ++i)
    {
        s = step_plenum(1.0e-4, cfg, s, comp, thr, *eos, integ, st);
        REQUIRE(st.ok);
        const double np = std::fabs(s.p_plenum - ps) / (ps - cfg.p_ambient)
                          + std::fabs(s.mdot_duct - ms) / ms;
        max_growth = std::max(max_growth, np / (1.0e-4 / ms));
    }
    CHECK(max_growth > 10.0);

    // CONTRAST (naive-slope check): move the throttle open so the same
    // compressor's operating point sits just below the peak where the slope is
    // still POSITIVE but small. The trace now goes negative -> STABLE. A naive
    // "positive slope is always unstable" argument would get this wrong; the
    // weighted criterion (s/I vs (a^2/V)*slope_throttle) gets it right.
    const double gain2 = 0.005628;
    const auto comp2 = comp;
    const auto thr2 = [&](double p) { return sqrt_throttle(p, cfg.p_ambient, gain2); };
    auto f2 = [&](double m) { return m - sqrt_throttle(cfg.p_ambient + cubic_dp(m, A, B, C),
                                                       cfg.p_ambient, gain2); };
    const double ms2 = bisect_positive(f2, 1.0, 1.2649);
    const double ps2 = cfg.p_ambient + (ms2 / gain2) * (ms2 / gain2);
    const double slope2 = B - 3.0 * C * ms2 * ms2;
    REQUIRE(slope2 > 0.0); // still a positive compressor slope
    PlenumJacobian j2 = plenum_jacobian(cfg, comp2, thr2, ps2, ms2, *eos);
    CHECK(j2.trace() < 0.0);   // positive slope yet linearly STABLE
    CHECK(j2.det() > 0.0);
}

TEST_CASE("plenum: sustained surge limit cycle (deep surge signature)")
{
    auto eos = make_air();

    // Big volume, long thin duct: the Greitzer B-parameter puts the system
    // into the surge regime with the cubic characteristic.
    PlenumModelConfig big;
    big.volume = 2.0;
    big.duct_area = 0.02;
    big.duct_length = 5.0;

    const double A = 40000.0;
    const double B = 12000.0;
    const double C = 2500.0;
    const double gain = 0.005;
    const auto comp = [&](double m) { return cubic_dp(m, A, B, C); };
    const auto thr = [&](double p) { return sqrt_throttle(p, big.p_ambient, gain); };
    auto f = [&](double m) { return m - sqrt_throttle(big.p_ambient + cubic_dp(m, A, B, C),
                                                      big.p_ambient, gain); };
    const double ms = bisect_positive(f, 0.5, 1.2649);
    const double ps = big.p_ambient + (ms / gain) * (ms / gain);

    // Start ON the unstable fixed point plus a tiny nudge.
    ModelStatus st;
    IntegratorConfig integ;
    integ.method = IntegrationMethod::RK4;
    constexpr int N = 400000;    // 40 s at dt = 1e-4; the nonlinear surge cycle
                                  // runs at ~0.74 s, so this is ~54 cycles.
    constexpr int late_step = 74000; // last ~10 cycles (~7.4 s)
    PlenumState s{ps, ms + 1.0e-4};

    double pmin = 1e300, pmax = -1e300;
    double lmin = 1e300, lmax = -1e300;
    double early_max = 0.0;
    double msum = 0.0;
    int mcount = 0;
    int peaks = 0;
    bool saw_peak = false;
    double prev_p = s.p_plenum;
    for (int i = 0; i < N; ++i)
    {
        s = step_plenum(1.0e-4, big, s, comp, thr, *eos, integ, st);
        REQUIRE(st.ok);
        REQUIRE(std::isfinite(s.p_plenum));
        REQUIRE(std::isfinite(s.mdot_duct));
        REQUIRE(s.p_plenum > 0.0);

        pmin = std::min(pmin, s.p_plenum);
        pmax = std::max(pmax, s.p_plenum);
        if (i < N / 10) early_max = std::max(early_max, s.p_plenum);
        if (i >= N - late_step)
        {
            lmin = std::min(lmin, s.p_plenum);
            lmax = std::max(lmax, s.p_plenum);
            msum += s.mdot_duct;
            ++mcount;
        }
        // Peak counting on the full-resolution pressure (limit cycle is clean).
        if (i > 0 && i < N - 1)
        {
            if (s.p_plenum > prev_p) saw_peak = true;
            else if (saw_peak && s.p_plenum < prev_p) { ++peaks; saw_peak = false; }
        }
        prev_p = s.p_plenum;
    }

    // 50+ cycles of sustained oscillation.
    CHECK(peaks > 50);

    // The oscillation does NOT decay: the late peak-to-peak pressure amplitude
    // stays above 50% of the ambient-referenced transient peak pressure.
    const double late_pp = lmax - lmin;
    CHECK(late_pp > 0.5 * (early_max - big.p_ambient));

    // Classic surge signature: the time-mean through-flow over the last cycles
    // collapses below the steady (unstable) operating-point flow.
    CHECK(msum / mcount < ms);

    // The cycle keeps the cell sane (positive pressure, finite states).
    CHECK(pmin > 0.0);
    CHECK(pmax > big.p_ambient);
}

TEST_CASE("plenum: domain validation reports errors, never throws")
{
    auto eos = make_air();

    auto comp = [](double m) { return 30000.0 - 4000.0 * m; };
    auto thr = [](double p) { return sqrt_throttle(p, 101325.0, 0.02); };

    std::string error;
    std::vector<std::string> warnings;
    PlenumModelConfig bad;
    bad.volume = 0.0;
    CHECK_FALSE(validate_plenum_config(bad, error, warnings));
    CHECK_FALSE(error.empty());

    ModelStatus st;
    PlenumDerivative d = plenum_derivative(bad, PlenumState{101325.0, 0.1}, comp,
                                           thr, *eos, st);
    CHECK_FALSE(st.ok);
    CHECK_FALSE(d.ok);
    CHECK_FALSE(d.status.error.empty());

    IntegratorConfig integ;
    PlenumState s = step_plenum(1.0e-4, bad, PlenumState{101325.0, 0.1}, comp,
                                thr, *eos, integ, st);
    CHECK_FALSE(st.ok);
    CHECK(s.p_plenum == 101325.0); // state untouched on failure

    // dt <= 0 is rejected by step_plenum.
    st.ok = true;
    st.error.clear();
    PlenumState s2 = step_plenum(0.0, PlenumModelConfig{}, PlenumState{101325.0, 0.1},
                                 comp, thr, *eos, integ, st);
    CHECK_FALSE(st.ok);
    CHECK(st.error.find("dt") != std::string::npos);
}