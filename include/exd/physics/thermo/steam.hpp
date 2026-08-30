#pragma once

// ─────────────────────────────────────────────────────
// Steam properties (engineering-grade saturation model).
//
// Saturation line via the Clausius–Clapeyron relation with
// a single latent-heat anchor — NOT IAPWS-IF97: the
// two-point approximation tracks the true saturation
// pressure to ~5–10% across 0–200 °C (documented
// engineering model; swap for IF97 tables when needed).
//
// Anchor: p0 = 101.325 kPa at T0 = 373.15 K (100 °C),
// h_fg = 2.257 MJ/kg, R = 461.5 J/(kg·K).
// ─────────────────────────────────────────────────────

#include <cmath>

namespace exd::physics::thermo {

struct SteamConstants
{
    double p0 = 101325.0;     // Pa, anchor saturation pressure
    double T0 = 373.15;       // K, anchor saturation temperature (100 °C)
    double h_fg = 2.257e6;    // J/kg latent heat at the anchor
    double r_steam = 461.5;   // J/(kg·K) steam gas constant
    double c_pw = 4216.0;     // J/(kg·K) liquid water specific heat
    double h_f0 = 419.1e3;    // J/kg saturated liquid enthalpy at T0
};

/// Saturation pressure at temperature T (Pa), from the
/// Clausius–Clapeyron anchor. Valid ~[273, 500] K.
inline double satured_pressure(double T, const SteamConstants& c = {})
{
    return c.p0 * std::exp(-(c.h_fg / c.r_steam) * (1.0 / T - 1.0 / c.T0));
}

/// Saturation temperature at pressure p (K); inverse of
/// satured_pressure. p must be > 0.
inline double saturation_temperature(double p, const SteamConstants& c = {})
{
    return 1.0 / (1.0 / c.T0 - (c.r_steam / c.h_fg) * std::log(p / c.p0));
}

/// Saturated liquid enthalpy at T (J/kg).
inline double h_f(double T, const SteamConstants& c = {})
{
    return c.h_f0 + c.c_pw * (T - c.T0);
}

/// Latent heat at T (J/kg) via the anchor (Clausius–Clapeyron consistency).
inline double h_fg(double T, const SteamConstants& c = {})
{
    return c.h_fg; // single-anchor engineering model
}

/// Saturated vapor enthalpy at T (J/kg).
inline double h_g(double T, const SteamConstants& c = {})
{
    return h_f(T, c) + h_fg(T, c);
}

/// Wet-steam enthalpy at T with dryness fraction x ∈ [0, 1] (J/kg).
inline double h_wet(double T, double x, const SteamConstants& c = {})
{
    return h_f(T, c) + x * h_fg(T, c);
}

/// Saturated vapor density (kg/m³) at saturation pressure p and
/// temperature T — ideal-gas vapor approximation (engineering model;
/// the true saturated vapor density is lower near the critical point).
inline double rho_g(double p, double T, const SteamConstants& c = {})
{
    return p / (c.r_steam * T);
}

} // namespace exd::physics::thermo
