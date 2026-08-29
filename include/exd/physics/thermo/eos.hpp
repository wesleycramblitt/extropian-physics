#pragma once

#include <exd/physics/model_status.hpp>

#include <memory>
#include <string_view>

namespace exd::physics::thermo {

// ─────────────────────────────────────────────────────
// Equations of state. First variant: ideal gas.
// Consumed by engines (in-cylinder), compressors (stage
// thermodynamics), and thermal/chmistry domains later.
// ─────────────────────────────────────────────────────

class IEos
{
public:
    virtual ~IEos() = default;
    virtual std::string_view name() const = 0;

    /// Pressure from density and temperature (Pa).
    virtual double pressure(double rho, double T, exd::physics::ModelStatus& status) const = 0;
    /// Temperature from pressure and density (K).
    virtual double temperature(double p, double rho, exd::physics::ModelStatus& status) const = 0;
    /// Specific internal energy e (J/kg), reference at 0 K.
    virtual double internal_energy(double rho, double T, exd::physics::ModelStatus& status) const = 0;
    /// Specific enthalpy h (J/kg), reference at 0 K.
    virtual double enthalpy(double rho, double T, exd::physics::ModelStatus& status) const = 0;

    /// Ratio of specific heats γ = cp/cv (constant for ideal gas).
    virtual double gamma() const = 0;
    /// Specific gas constant R (J/(kg·K)).
    virtual double gas_constant() const = 0;
    /// Specific heats (J/(kg·K)).
    virtual double specific_heat_cv() const = 0;
    virtual double specific_heat_cp() const = 0;
};

struct IdealGasConfig
{
    double R = 287.05;   // specific gas constant, J/(kg·K) (air default)
    double gamma = 1.4;  // ratio of specific heats (> 1)
};

/// p = ρ·R·T; e = cv·T; h = cp·T.
/// Returns nullptr when R <= 0 or gamma <= 1 (invalid config).
std::unique_ptr<IEos> make_ideal_gas(const IdealGasConfig& config);

} // namespace exd::physics::thermo
