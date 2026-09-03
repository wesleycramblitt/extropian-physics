#pragma once

// ─────────────────────────────────────────────────────
// Dimensional units (implementation_spec §25).
//
// A Unit is a vector of SI base-dimension exponents
// (m, kg, s, A, K, mol, cd).  Quantities carry units and
// interfaces validate dimensional consistency at build
// time where possible and at configuration/coupling time
// otherwise.  Unit conversions happen at interfaces, not
// inside performance-critical kernels.
// ─────────────────────────────────────────────────────

#include <cstdint>
#include <string>

namespace exd::engine::core {

struct Unit
{
    int8_t m = 0;    // length (meter)
    int8_t kg = 0;   // mass (kilogram)
    int8_t s = 0;    // time (second)
    int8_t A = 0;    // electric current (ampere)
    int8_t K = 0;    // temperature (kelvin)
    int8_t mol = 0;  // amount of substance (mole)
    int8_t cd = 0;   // luminous intensity (candela)

    constexpr bool is_dimensionless() const
    {
        return m == 0 && kg == 0 && s == 0 && A == 0 && K == 0 && mol == 0 && cd == 0;
    }

    constexpr bool operator==(const Unit&) const = default;

    constexpr Unit operator*(const Unit& o) const
    {
        return {static_cast<int8_t>(m + o.m),     static_cast<int8_t>(kg + o.kg),
                static_cast<int8_t>(s + o.s),     static_cast<int8_t>(A + o.A),
                static_cast<int8_t>(K + o.K),     static_cast<int8_t>(mol + o.mol),
                static_cast<int8_t>(cd + o.cd)};
    }

    constexpr Unit operator/(const Unit& o) const
    {
        return {static_cast<int8_t>(m - o.m),     static_cast<int8_t>(kg - o.kg),
                static_cast<int8_t>(s - o.s),     static_cast<int8_t>(A - o.A),
                static_cast<int8_t>(K - o.K),     static_cast<int8_t>(mol - o.mol),
                static_cast<int8_t>(cd - o.cd)};
    }

    constexpr Unit pow(int n) const
    {
        return {static_cast<int8_t>(m * n),     static_cast<int8_t>(kg * n),
                static_cast<int8_t>(s * n),     static_cast<int8_t>(A * n),
                static_cast<int8_t>(K * n),     static_cast<int8_t>(mol * n),
                static_cast<int8_t>(cd * n)};
    }
};

/// SI base and derived unit constants.
namespace units {
inline constexpr Unit dimensionless{};
inline constexpr Unit meter{1, 0, 0, 0, 0, 0, 0};
inline constexpr Unit kilogram{0, 1, 0, 0, 0, 0, 0};
inline constexpr Unit second{0, 0, 1, 0, 0, 0, 0};
inline constexpr Unit ampere{0, 0, 0, 1, 0, 0, 0};
inline constexpr Unit kelvin{0, 0, 0, 0, 1, 0, 0};
inline constexpr Unit mole{0, 0, 0, 0, 0, 1, 0};
inline constexpr Unit candela{0, 0, 0, 0, 0, 0, 1};

inline constexpr Unit hertz = second.pow(-1);
inline constexpr Unit newton = kilogram * meter / second.pow(2);
inline constexpr Unit pascal = newton / meter.pow(2);
inline constexpr Unit joule = newton * meter;
inline constexpr Unit watt = joule / second;
inline constexpr Unit volt = watt / ampere;
inline constexpr Unit tesla = kilogram / (ampere * second.pow(2));
inline constexpr Unit coulomb = ampere * second;
inline constexpr Unit ohm = volt / ampere;
inline constexpr Unit farad = coulomb / volt;
inline constexpr Unit velocity = meter / second;
inline constexpr Unit acceleration = meter / second.pow(2);
inline constexpr Unit density = kilogram / meter.pow(3);
inline constexpr Unit force = newton;
inline constexpr Unit pressure = pascal;
inline constexpr Unit energy = joule;
inline constexpr Unit power = watt;
inline constexpr Unit torque = newton * meter;
inline constexpr Unit heat_flux = watt / meter.pow(2);
inline constexpr Unit thermal_conductivity = watt / (meter * kelvin);
inline constexpr Unit specific_heat = joule / (kilogram * kelvin);
inline constexpr Unit dynamic_viscosity = pascal * second;
inline constexpr Unit kinematic_viscosity = meter.pow(2) / second;
inline constexpr Unit angular_velocity = second.pow(-1);
} // namespace units

/// True when `a` and `b` have identical dimensions (spec §25: dimensional
/// validation of couplings and interfaces).
constexpr bool units_compatible(const Unit& a, const Unit& b) { return a == b; }

/// Human-readable SI representation, e.g. "kg·m⁻¹·s⁻²" for pressure.
inline std::string to_string(const Unit& u)
{
    const char* syms[] = {"m", "kg", "s", "A", "K", "mol", "cd"};
    const int8_t exps[] = {u.m, u.kg, u.s, u.A, u.K, u.mol, u.cd};
    std::string out;
    for (int i = 0; i < 7; ++i)
    {
        if (exps[i] == 0) continue;
        if (!out.empty()) out += "·";
        out += syms[i];
        const int e = exps[i];
        if (e != 1)
        {
            // UTF-8 superscripts — indexed as STRINGS, not bytes
            static const char* const sup[] = {"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
            if (e < 0) out += "⁻";
            const int mag = e < 0 ? -e : e;
            if (mag < 10) out += sup[mag];
            else out += std::to_string(e);   // fallback for large exponents
        }
    }
    return out.empty() ? "—" : out;
}

/// Quantity = value + unit metadata (configuration/interfaces; hot paths use
/// raw doubles — spec §25: conversions at interfaces).
struct Quantity
{
    double value = 0.0;
    Unit unit{};
};

constexpr Quantity operator*(double v, const Unit& u) { return {v, u}; }

} // namespace exd::engine::core
