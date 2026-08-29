#pragma once

namespace exd::physics::thermo {

// ─────────────────────────────────────────────────────
// Thermophysical property helpers.
// ─────────────────────────────────────────────────────

/// Sutherland's law: mu = mu_ref · (T/T_ref)^1.5 · (T_ref + S)/(T + S).
/// Defaults are air: T_ref = 273.15 K, mu_ref = 1.716e-5 Pa·s, S = 110.4 K.
double sutherland_viscosity(double T,
                            double T_ref = 273.15,
                            double mu_ref = 1.716e-5,
                            double S = 110.4);

} // namespace exd::physics::thermo
