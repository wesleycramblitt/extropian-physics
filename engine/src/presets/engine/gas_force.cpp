// gas_force.cpp
// Cycle pressure/temperature model — PURE function of the
// crank angle (no persisted cycle state → restartable,
// optimizer-batchable, deterministic under any integrator).
//
// Otto (4-stroke; theta_c = mod(theta, 4π)):
//   compression  540°–720°: p = p_intake·(V_bdc/V)^γc
//   power          0°–180°: polytropic expansion from
//                           p_TDC = p_intake·(V_bdc/V_tdc)^γc
//                           + Wiebe heat-release pressure
//                           (γe−1)·Q_in·x_b(θ)/V
//   exhaust      180°–360°: sin² blend from p_exp_end → p_exhaust
//   intake       360°–540°: sin² blend → p_intake
//
// Steam (single-acting; theta_c = mod(theta, 2π)):
//   admission 0°–cutoff: p = p_boiler (sin² blend in)
//   expansion cutoff°–180°: p = p_boiler·(V_cut/V)^n, n = steam_gamma
//   exhaust   180°–360°: sin² blend → p_condenser
//
// Temperatures: trapped mass m = p_intake·V_bdc/(R·T_intake)
// (Otto) / p_boiler·V_cut/(R·T_intake) (steam); T = pV/(mR).
// Valve-open phases: temperature is nominal (open cylinder).

#include "engine_internal.hpp"

#include <exd/engine/physics/thermo/steam.hpp>

#include <cmath>

namespace exd::engine::presets::engine {

namespace
{
constexpr double PI = 3.14159265358979323846;
constexpr double DEG = PI / 180.0;

/// sin² ramp 0→1 over fraction f ∈ [0,1] (smooth, C¹).
double smooth_blend(double f)
{
    if (f <= 0.0) return 0.0;
    if (f >= 1.0) return 1.0;
    const double s = std::sin(0.5 * PI * f);
    return s * s;
}

/// Wiebe cumulative heat-release fraction.
double wiebe_fraction(double theta_deg, const EngineThermoConfig& t)
{
    const double start = t.wiebe_ignition_deg;
    const double span = t.wiebe_burn_duration_deg;
    if (span <= 0.0) return theta_deg >= start ? 1.0 : 0.0;
    if (theta_deg <= start) return 0.0;
    if (theta_deg >= start + span) return 1.0;
    const double x = (theta_deg - start) / span;
    return 1.0 - std::exp(-t.wiebe_a * std::pow(x, t.wiebe_m + 1.0));
}

double blend_window_deg()
{
    return 15.0; // valve ramp width
}
} // anonymous namespace

void cylinder_state(double theta,
                    const EngineGeometryConfig& g,
                    const EngineThermoConfig& t,
                    double& p_cyl,
                    double& T_cyl)
{
    const double A = piston_area(g);
    const double r = g.crank_radius;
    const double l = g.rod_length;
    const double v_tdc = g.clearance_volume;

    if (t.cycle == EngineCycleType::Steam)
    {
        const double period = 2.0 * PI;
        double tc = std::fmod(theta, period);
        if (tc < 0.0) tc += period;
        const double deg = tc / DEG;

        const double mc = t.steam_cutoff_deg * DEG;
        const double v_cut = cylinder_volume(mc, g);
        const double p_180 = t.p_boiler * std::pow(v_cut / cylinder_volume(PI, g), t.steam_gamma);

        double p = t.p_condenser;
        if (deg >= 0.0 && deg < t.steam_cutoff_deg)
        {
            const double f = smooth_blend(deg / 3.0); // 3° admission ramp
            p = t.p_condenser + f * (t.p_boiler - t.p_condenser);
        }
        else if (deg >= t.steam_cutoff_deg && deg < 180.0)
        {
            const double v = cylinder_volume(theta, g);
            p = t.p_boiler * std::pow(v_cut / v, t.steam_gamma);
        }
        else if (deg >= 180.0 && deg < 180.0 + blend_window_deg())
        {
            const double f = smooth_blend((deg - 180.0) / blend_window_deg());
            p = p_180 + f * (t.p_condenser - p_180);
        }
        else
        {
            p = t.p_condenser;
        }
        p_cyl = p;

        // ── Steam temperature: saturation-based when wet ──
        // Dryness fraction during expansion follows the standard
        // wet-polytrope approximation x = x_cut·(V_cut/V)^(n−1); while
        // x < 1 the cylinder gas sits on the saturation line (T = T_sat(p)),
        // above it (superheated) the ideal-gas law applies with the trapped
        // mass m = ρ_g(p_b, T_sat(p_b))·V_cut·x_cut.
        const double x = deg >= t.steam_cutoff_deg && deg < 180.0
            ? t.steam_quality_cutoff
                  * std::pow(v_cut / cylinder_volume(theta, g), t.steam_gamma - 1.0)
            : deg < t.steam_cutoff_deg ? t.steam_quality_cutoff
                                       : t.steam_quality_cutoff
                                             * std::pow(v_cut / cylinder_volume(PI, g),
                                                        t.steam_gamma - 1.0);
        if (x < 1.0)
        {
            T_cyl = exd::engine::physics::thermo::saturation_temperature(p);
        }
        else
        {
            const double t_sat_b = exd::engine::physics::thermo::saturation_temperature(t.p_boiler);
            const double m = exd::engine::physics::thermo::rho_g(t.p_boiler, t_sat_b) * v_cut
                             * t.steam_quality_cutoff;
            T_cyl = p * cylinder_volume(theta, g) / (m * t.r_gas);
        }
        return;
    }

    // ── Otto, 4-stroke ────────────────────────────────
    const double period = 4.0 * PI;
    double tc = std::fmod(theta, period);
    if (tc < 0.0) tc += period;
    const double deg = tc / DEG;

    const double v_bdc = g.clearance_volume + A * 2.0 * r; // V(180°)
    const double p_tdc = t.p_intake * std::pow(v_bdc / v_tdc, t.gamma_compression);

    const double v = cylinder_volume(theta, g);
    double p = t.p_exhaust;

    if (deg >= 540.0 && deg < 720.0)
    {
        // Compression
        p = t.p_intake * std::pow(v_bdc / v, t.gamma_compression);
    }
    else if (deg >= 0.0 && deg < 180.0)
    {
        // Power: polytropic expansion + Wiebe heat-release pressure.
        // The heat is anchored at TDC volume and expands polytropically
        // ((γ−1)·q·x_b/V_tdc · (V_tdc/V)^γe): a raw (γ−1)·q·x_b/V term
        // decays ∝ 1/V and over-produces work beyond the Otto bound
        // (measured η ≈ 0.55 vs the 0.489 limit at r_c = 6.8, γ = 1.35).
        p = p_tdc * std::pow(v_tdc / v, t.gamma_expansion);
        const double xb = wiebe_fraction(deg, t);
        p += (t.gamma_expansion - 1.0) * t.q_in_cycle * xb / v_tdc
             * std::pow(v_tdc / v, t.gamma_expansion);
    }
    else if (deg >= 180.0 && deg < 180.0 + blend_window_deg())
    {
        // Exhaust opening ramp: p_exp_end → p_exhaust (xb ≡ 1 here)
        const double p_end = p_tdc * std::pow(v_tdc / v_bdc, t.gamma_expansion)
                             + (t.gamma_expansion - 1.0) * t.q_in_cycle / v_tdc
                                   * std::pow(v_tdc / v_bdc, t.gamma_expansion);
        const double f = smooth_blend((deg - 180.0) / blend_window_deg());
        p = p_end + f * (t.p_exhaust - p_end);
    }
    else if (deg >= 360.0 && deg < 360.0 + blend_window_deg())
    {
        // Intake opening ramp: p_exhaust → p_intake
        const double f = smooth_blend((deg - 360.0) / blend_window_deg());
        p = t.p_exhaust + f * (t.p_intake - t.p_exhaust);
    }
    else
    {
        p = t.p_exhaust; // exhaust hold / intake hold
    }
    p_cyl = p;

    const double m = t.p_intake * v_bdc / (t.r_gas * t.T_intake);
    T_cyl = p * v / (m * t.r_gas);
}

double load_moment(const EngineLoadConfig& load, double omega)
{
    double m = load.friction_constant + load.friction_viscous * omega;
    if (load.generator_enabled)
    {
        const auto& wp = load.generator_omega_pts;
        const auto& tp = load.generator_torque_pts;
        if (wp.size() == tp.size() && wp.size() >= 2)
        {
            if (omega <= wp.front()) m += tp.front();
            else if (omega >= wp.back()) m += tp.back();
            else
            {
                // piecewise-linear
                for (size_t i = 1; i < wp.size(); ++i)
                {
                    if (omega <= wp[i])
                    {
                        const double f = (omega - wp[i - 1]) / (wp[i] - wp[i - 1]);
                        m += tp[i - 1] + f * (tp[i] - tp[i - 1]);
                        break;
                    }
                }
            }
        }
    }
    return m;
}

} // namespace exd::engine::presets::engine
