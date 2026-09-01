// polytropic_test.cpp
// Stagnation-polytrope primitives: exponent round-trips, isentropic limit,
// compression/expansion asymmetry (polytropes are NOT reciprocals), and
// entropy-consistency of the stage temperature/pressure chain.

#include <exd/physics/thermo/polytropic.hpp>

#include <doctest/doctest.h>

#include <cmath>

using namespace exd::physics::thermo;
using exd::physics::ModelStatus;

TEST_CASE("Polytropic: eta_p = 1 reproduces the exact isentropic relations")
{
    ModelStatus status;
    const double gamma = 1.4;
    const double pi = 5.0;

    const double tau = polytropic::temp_ratio_compression(pi, gamma, 1.0, status);
    CHECK(status.ok);
    CHECK(tau == doctest::Approx(std::pow(pi, (gamma - 1.0) / gamma)).epsilon(1e-12));

    const double pi_back = polytropic::pressure_ratio_compression(tau, gamma, 1.0, status);
    CHECK(pi_back == doctest::Approx(pi).epsilon(1e-12));

    // Expansion at eta = 1: same form as compression at eta = 1.
    const double tau_e = polytropic::temp_ratio_expansion(pi, gamma, 1.0, status);
    CHECK(tau_e == doctest::Approx(std::pow(pi, (gamma - 1.0) / gamma)).epsilon(1e-12));
}

TEST_CASE("Polytropic: compression/expansion round-trips at eta_p < 1")
{
    ModelStatus status;
    const double gamma = 1.4;
    const double eta = 0.85;
    const double pi = 2.0;

    const double tau_c = polytropic::temp_ratio_compression(pi, gamma, eta, status);
    CHECK(status.ok);
    // Inverse recovers pi for both senses.
    CHECK(polytropic::pressure_ratio_compression(tau_c, gamma, eta, status) ==
          doctest::Approx(pi).epsilon(1e-12));

    const double tau_e = polytropic::temp_ratio_expansion(pi, gamma, eta, status);
    CHECK(polytropic::pressure_ratio_expansion(tau_e, gamma, eta, status) ==
          doctest::Approx(pi).epsilon(1e-12));

    // Compression heats more than expansion cools for the same pi (eta < 1).
    CHECK(tau_c > tau_e);

    // Efficiency recovery round-trips.
    CHECK(polytropic::polytropic_efficiency_compression(pi, tau_c, gamma, status) ==
          doctest::Approx(eta).epsilon(1e-12));
    CHECK(polytropic::polytropic_efficiency_expansion(pi, tau_e, gamma, status) ==
          doctest::Approx(eta).epsilon(1e-12));
}

TEST_CASE("Polytropic: compressor and turbine polytropes are NOT reciprocal")
{
    // The classic wrong implementation: pi_turb = 1/pi_comp with symmetric
    // exponents. The correct expansion at eta_p < 1 gives a DIFFERENT result.
    ModelStatus status;
    const double gamma = 1.4;
    const double eta = 0.85;
    const double pi_comp = 1.309; // compression π for Δh0 = 30 kJ/kg at 300 K

    const double tau_c = polytropic::temp_ratio_compression(pi_comp, gamma, eta, status);
    const double tau_e = polytropic::temp_ratio_expansion(1.0 / pi_comp, gamma, eta, status);

    // pi_comp · pi_turb != 1: tau_e != 1/tau_c on a reciprocal implementation.
    // Check the equivalent: symmetric expo would give tau_e = tau_c^(-1).
    CHECK(std::abs(tau_e - 1.0 / tau_c) > 1e-3);
    // And the physically correct relation holds:
    CHECK(tau_e == doctest::Approx(std::pow(1.0 / pi_comp, eta * (gamma - 1.0) / gamma))
                       .epsilon(1e-12));

    // Same |Δh0|, opposite sense keeps entropy production positive both ways:
    // ds/cv policy check: for compression, tau > pi^((gamma-1)/gamma) (isentropic).
    const double tau_isen = std::pow(pi_comp, (gamma - 1.0) / gamma);
    CHECK(tau_c > tau_isen);   // lossy compression heats MORE than isentropic
    // Lossy expansion cools LESS than isentropic (tau_e > isen limit); both
    // senses produce positive entropy (checked explicitly below).
    CHECK(tau_e > std::pow(1.0 / pi_comp, (gamma - 1.0) / gamma));
}

TEST_CASE("Polytropic: entropy production is positive for real compression/expansion")
{
    ModelStatus status;
    const double gamma = 1.4;
    const double eta = 0.85;
    const double R = 287.05, cp = 1004.5;
    const double T1 = 288.15, p1 = 101325.0;

    // Compression: ds = cp·ln(tau) − R·ln(pi) > 0 at eta < 1.
    const double pi = 2.0;
    const double tau = polytropic::temp_ratio_compression(pi, gamma, eta, status);
    const double ds_c = cp * std::log(tau) - R * std::log(pi);
    CHECK(ds_c > 0.0);

    // Expansion across pi_out = 1/2: reversed flow, same magnitude of loss.
    const double tau_e = polytropic::temp_ratio_expansion(0.5, gamma, eta, status);
    const double ds_e = cp * std::log(tau_e) - R * std::log(0.5);
    CHECK(ds_e > 0.0);
}

TEST_CASE("Polytropic: domain validation reports errors, never throws")
{
    ModelStatus status;

    // gamma <= 1
    const double bad = polytropic::temp_ratio_compression(2.0, 1.0, 0.85, status);
    CHECK(!status.ok);
    CHECK(bad == 0.0);

    // eta_poly out of (0, 1]
    status.ok = true; status.error.clear();
    polytropic::temp_ratio_compression(2.0, 1.4, 1.5, status);
    CHECK(!status.ok);

    // non-positive ratios
    status.ok = true; status.error.clear();
    polytropic::temp_ratio_compression(-1.0, 1.4, 0.85, status);
    CHECK(!status.ok);
}
