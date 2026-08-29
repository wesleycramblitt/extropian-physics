// Ideal gas equation of state. The concrete model is hidden in an anonymous
// namespace and reachable only through the make_ideal_gas factory function.

#include <exd/physics/thermo/eos.hpp>

#include <memory>

namespace exd::physics::thermo
{

namespace
{

// ── IdealGasEos ───────────────────────────────────────────────
// Thermally perfect ideal gas: p = rho·R·T with constant specific heats
// cv = R/(gamma - 1), cp = gamma·cv. e and h reference 0 K.

class IdealGasEos final : public IEos
{
public:
    IdealGasEos(double R, double gamma)
        : R_(R), gamma_(gamma), cv_(R / (gamma - 1.0)),
          cp_(gamma * R / (gamma - 1.0))
    {
    }

    std::string_view name() const override { return "ideal_gas"; }

    double pressure(double rho, double T, ModelStatus& status) const override
    {
        if (rho <= 0.0 || T <= 0.0)
        {
            status.ok = false;
            status.error = "ideal gas: density and temperature must be positive";
            return 0.0;
        }
        return rho * R_ * T;
    }

    double temperature(double p, double rho, ModelStatus& status) const override
    {
        if (rho <= 0.0 || p < 0.0)
        {
            status.ok = false;
            status.error = "ideal gas: density must be positive and pressure non-negative";
            return 0.0;
        }
        return p / (rho * R_);
    }

    double internal_energy(double, double T, ModelStatus& status) const override
    {
        if (T <= 0.0)
        {
            status.ok = false;
            status.error = "ideal gas: temperature must be positive";
            return 0.0;
        }
        return cv_ * T;
    }

    double enthalpy(double, double T, ModelStatus& status) const override
    {
        if (T <= 0.0)
        {
            status.ok = false;
            status.error = "ideal gas: temperature must be positive";
            return 0.0;
        }
        return cp_ * T;
    }

    double gamma() const override { return gamma_; }
    double gas_constant() const override { return R_; }
    double specific_heat_cv() const override { return cv_; }
    double specific_heat_cp() const override { return cp_; }

private:
    double R_;
    double gamma_;
    double cv_;
    double cp_;
};

} // anonymous namespace

// ── Factory ───────────────────────────────────────────────────

std::unique_ptr<IEos> make_ideal_gas(const IdealGasConfig& config)
{
    if (config.R <= 0.0 || config.gamma <= 1.0)
        return nullptr;
    return std::make_unique<IdealGasEos>(config.R, config.gamma);
}

} // namespace exd::physics::thermo