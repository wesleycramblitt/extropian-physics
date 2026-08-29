// Thermophysical property helpers.
//
// Sutherland's law returns 0.0 for T <= 0 — a degenerate input the formula
// cannot represent (it is undefined at T = -S and meaningless below). There
// is no status channel in this signature, so the zero return value signals
// the degenerate case (documented).

#include <exd/physics/thermo/properties.hpp>

#include <cmath>

namespace exd::physics::thermo
{

double sutherland_viscosity(double T, double T_ref, double mu_ref, double S)
{
    if (T <= 0.0)
        return 0.0;
    return mu_ref * std::pow(T / T_ref, 1.5) * (T_ref + S) / (T + S);
}

} // namespace exd::physics::thermo