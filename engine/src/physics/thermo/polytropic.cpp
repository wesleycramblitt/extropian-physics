// Polytropic process primitives — stagnation family (see header).

#include <exd/engine/physics/thermo/polytropic.hpp>

#include <algorithm>
#include <cmath>

namespace exd::engine::physics::thermo::polytropic {

namespace {

// Shared domain validation for the exponent forms.
bool valid_arguments(double gamma, double eta_poly, exd::engine::core::ModelStatus& status)
{
    if (gamma <= 1.0)
    {
        status.ok = false;
        status.error = "polytropic: gamma must be > 1";
        return false;
    }
    if (eta_poly <= 0.0 || eta_poly > 1.0)
    {
        status.ok = false;
        status.error = "polytropic: polytropic efficiency must be in (0, 1]";
        return false;
    }
    return true;
}

bool valid_pi(double pi, exd::engine::core::ModelStatus& status)
{
    if (pi <= 0.0)
    {
        status.ok = false;
        status.error = "polytropic: pressure ratio must be positive";
        return false;
    }
    return true;
}

bool valid_tau(double tau, exd::engine::core::ModelStatus& status)
{
    if (tau <= 0.0)
    {
        status.ok = false;
        status.error = "polytropic: temperature ratio must be positive";
        return false;
    }
    return true;
}

} // anonymous namespace

double temp_ratio_compression(double pi, double gamma, double eta_poly,
                              exd::engine::core::ModelStatus& status)
{
    if (!(valid_arguments(gamma, eta_poly, status) && valid_pi(pi, status)))
        return 0.0;
    const double exponent = (gamma - 1.0) / (gamma * eta_poly);
    return std::pow(pi, exponent);
}

double temp_ratio_expansion(double pi, double gamma, double eta_poly,
                            exd::engine::core::ModelStatus& status)
{
    if (!(valid_arguments(gamma, eta_poly, status) && valid_pi(pi, status)))
        return 0.0;
    const double exponent = eta_poly * (gamma - 1.0) / gamma;
    return std::pow(pi, exponent);
}

double pressure_ratio_compression(double tau, double gamma, double eta_poly,
                                  exd::engine::core::ModelStatus& status)
{
    if (!(valid_arguments(gamma, eta_poly, status) && valid_tau(tau, status)))
        return 0.0;
    const double exponent = gamma * eta_poly / (gamma - 1.0);
    return std::pow(tau, exponent);
}

double pressure_ratio_expansion(double tau, double gamma, double eta_poly,
                                exd::engine::core::ModelStatus& status)
{
    if (!(valid_arguments(gamma, eta_poly, status) && valid_tau(tau, status)))
        return 0.0;
    const double exponent = gamma / ((gamma - 1.0) * eta_poly);
    return std::pow(tau, exponent);
}

double polytropic_efficiency_compression(double pi, double tau, double gamma,
                                         exd::engine::core::ModelStatus& status)
{
    if (!(gamma > 1.0 && valid_pi(pi, status) && valid_tau(tau, status)))
        return 0.0;
    // τ = π^((γ−1)/(γ·η))  →  η = (γ−1)·ln(π)/(γ·ln(τ))
    double eta = (gamma - 1.0) * std::log(pi) / (gamma * std::log(tau));
    return std::clamp(eta, 0.0, 1.0);
}

double polytropic_efficiency_expansion(double pi, double tau, double gamma,
                                       exd::engine::core::ModelStatus& status)
{
    if (!(gamma > 1.0 && valid_pi(pi, status) && valid_tau(tau, status)))
        return 0.0;
    // τ = π^(η·(γ−1)/γ)  →  η = γ·ln(τ)/((γ−1)·ln(π))
    double eta = gamma * std::log(tau) / ((gamma - 1.0) * std::log(pi));
    return std::clamp(eta, 0.0, 1.0);
}

} // namespace exd::engine::physics::thermo::polytropic
