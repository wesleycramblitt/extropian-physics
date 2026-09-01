// eos_test.cpp
// Unit tests for the ideal gas equation of state (eos.hpp) and the
// thermophysical property helpers (properties.hpp).

#include <exd/physics/thermo/eos.hpp>
#include <exd/physics/thermo/properties.hpp>

#include <doctest/doctest.h>

#include <cmath>

using namespace exd::physics::thermo;

TEST_CASE("Ideal gas EOS: standard sea-level air state relations")
{
    const IdealGasConfig config; // air defaults: R = 287.05, gamma = 1.4
    auto eos = make_ideal_gas(config);
    REQUIRE(eos);
    CHECK(eos->name() == "ideal_gas");

    exd::physics::ModelStatus status;
    const double rho = 1.225;
    const double T = 288.15;

    // p = rho·R·T = 101323.99 Pa, a hair under one standard atmosphere.
    const double p = eos->pressure(rho, T, status);
    CHECK(status.ok);
    CHECK(p == doctest::Approx(101325.0).epsilon(0.001));

    // T = p / (rho·R) recovers the recorded atmospheric temperature.
    const double T_back = eos->temperature(101325.0, rho, status);
    CHECK(status.ok);
    CHECK(T_back == doctest::Approx(288.15).epsilon(1e-3));

    // Round trip: p -> T -> p is internally consistent.
    const double T_rt = eos->temperature(p, rho, status);
    const double p_rt = eos->pressure(rho, T_rt, status);
    CHECK(p_rt == doctest::Approx(p).epsilon(1e-12));
    CHECK(T_rt == doctest::Approx(T).epsilon(1e-9));
}

TEST_CASE("Ideal gas EOS: specific heats follow cv = R/(gamma-1) and cp = gamma*cv")
{
    const IdealGasConfig config; // R = 287.05, gamma = 1.4
    auto eos = make_ideal_gas(config);
    REQUIRE(eos);

    CHECK(eos->gamma() == doctest::Approx(1.4).epsilon(1e-12));
    CHECK(eos->gas_constant() == doctest::Approx(287.05).epsilon(1e-12));
    CHECK(eos->specific_heat_cv() == doctest::Approx(717.6).epsilon(1e-3));
    CHECK(eos->specific_heat_cp() == doctest::Approx(1004.5).epsilon(1e-3));
    CHECK(eos->specific_heat_cp() ==
          doctest::Approx(1.4 * eos->specific_heat_cv()).epsilon(1e-12));
}

TEST_CASE("Ideal gas EOS: internal energy and enthalpy scale linearly with temperature")
{
    const IdealGasConfig config; // R = 287.05, gamma = 1.4
    auto eos = make_ideal_gas(config);
    REQUIRE(eos);

    exd::physics::ModelStatus status;
    const double T = 300.0;
    const double cv = eos->specific_heat_cv();
    const double cp = eos->specific_heat_cp();

    CHECK(eos->internal_energy(1.225, T, status) == doctest::Approx(cv * T).epsilon(1e-12));
    CHECK(eos->enthalpy(1.225, T, status) == doctest::Approx(cp * T).epsilon(1e-12));
    CHECK(status.ok);

    // Density is not part of the thermally perfect ideal gas e(T)/h(T).
    CHECK(eos->internal_energy(2.5, T, status) == doctest::Approx(cv * T).epsilon(1e-12));
    CHECK(eos->enthalpy(0.1, T, status) == doctest::Approx(cp * T).epsilon(1e-12));
    CHECK(status.ok);
}

TEST_CASE("Ideal gas EOS: invalid state inputs fail the status channel")
{
    auto eos = make_ideal_gas(IdealGasConfig{});
    REQUIRE(eos);
    exd::physics::ModelStatus status;

    SUBCASE("pressure rejects non-positive density")
    {
        status = {};
        CHECK(eos->pressure(0.0, 288.15, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK(status.error == "ideal gas: density and temperature must be positive");

        status = {};
        CHECK(eos->pressure(-1.225, 288.15, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("pressure rejects non-positive temperature")
    {
        status = {};
        CHECK(eos->pressure(1.225, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);

        status = {};
        CHECK(eos->pressure(1.225, -300.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("temperature rejects non-positive density")
    {
        status = {};
        CHECK(eos->temperature(101325.0, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);

        status = {};
        CHECK(eos->temperature(101325.0, -1.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("temperature rejects negative pressure")
    {
        status = {};
        CHECK(eos->temperature(-1.0, 1.225, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("internal_energy rejects non-positive temperature")
    {
        status = {};
        CHECK(eos->internal_energy(1.225, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);

        status = {};
        CHECK(eos->internal_energy(1.225, -100.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("enthalpy rejects non-positive temperature")
    {
        status = {};
        CHECK(eos->enthalpy(1.225, 0.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);

        status = {};
        CHECK(eos->enthalpy(1.225, -100.0, status) == doctest::Approx(0.0));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }
}

TEST_CASE("Ideal gas EOS: factory rejects invalid configurations")
{
    SUBCASE("default air configuration is valid")
    {
        REQUIRE(make_ideal_gas(IdealGasConfig{}));
    }

    SUBCASE("R <= 0 returns nullptr")
    {
        IdealGasConfig config;
        config.R = 0.0;
        CHECK_FALSE(make_ideal_gas(config));

        config.R = -287.05;
        CHECK_FALSE(make_ideal_gas(config));
    }

    SUBCASE("gamma <= 1 returns nullptr")
    {
        IdealGasConfig config;
        config.gamma = 1.0;
        CHECK_FALSE(make_ideal_gas(config));

        config.gamma = 0.9;
        CHECK_FALSE(make_ideal_gas(config));
    }
}

TEST_CASE("Sutherland viscosity: returns the reference viscosity at T_ref")
{
    // Defaults are air: T_ref = 273.15 K, mu_ref = 1.716e-5 Pa·s, S = 110.4 K.
    CHECK(sutherland_viscosity(273.15) == doctest::Approx(1.716e-5).epsilon(1e-9));
    CHECK(sutherland_viscosity(273.15, 273.15, 1.716e-5, 110.4) ==
          doctest::Approx(1.716e-5).epsilon(1e-9));
}

TEST_CASE("Sutherland viscosity: air at 288.15 K matches the published value")
{
    const double mu = sutherland_viscosity(288.15);
    CHECK(mu == doctest::Approx(1.789e-5).epsilon(0.02)); // within 2%
    CHECK(mu > 1.716e-5); // viscosity increases with temperature
}

TEST_CASE("Sutherland viscosity: non-positive temperature is a degenerate input returning zero")
{
    CHECK(sutherland_viscosity(0.0) == doctest::Approx(0.0));
    CHECK(sutherland_viscosity(-50.0) == doctest::Approx(0.0));
}
TEST_CASE("Ideal gas EOS: density(p,T) inverts pressure(rho,T)")
{
    const IdealGasConfig config; // air defaults
    auto eos = make_ideal_gas(config);
    REQUIRE(eos);

    exd::physics::ModelStatus status;
    const double rho = 1.225;
    const double T = 288.15;
    const double p = eos->pressure(rho, T, status);
    REQUIRE(status.ok);

    const double rho_back = eos->density(p, T, status);
    CHECK(status.ok);
    CHECK(rho_back == doctest::Approx(rho).epsilon(1e-12));
}

TEST_CASE("Ideal gas EOS: specific_entropy matches analytic differences")
{
    const IdealGasConfig config; // R = 287.05, gamma = 1.4 -> cp = 1004.5
    auto eos = make_ideal_gas(config);
    REQUIRE(eos);

    exd::physics::ModelStatus status;
    // Isentropic compression p2/p1 = (T2/T1)^(gamma/(gamma-1)) must give ds = 0.
    const double T1 = 288.15, p1 = 101325.0;
    const double tau = 1.2;
    const double pi = std::pow(tau, 1.4 / (1.4 - 1.0));
    const double s1 = eos->specific_entropy(p1, T1, status);
    REQUIRE(status.ok);
    const double s2 = eos->specific_entropy(pi * p1, tau * T1, status);
    REQUIRE(status.ok);
    CHECK(s2 - s1 == doctest::Approx(0.0).epsilon(1e-9));

    // Same temperature, higher pressure: entropy must DECREASE (ds < 0).
    const double s3 = eos->specific_entropy(2.0 * p1, T1, status);
    REQUIRE(status.ok);
    CHECK(s3 - s1 < 0.0);

    // Analytic value at fixed (p,T): s = cp·ln(T/Tref) − R·ln(p/pref).
    const double expect = 1004.5 * std::log(T1 / 298.15) - 287.05 * std::log(p1 / 101325.0);
    CHECK(s1 == doctest::Approx(expect).epsilon(1e-3));
}
