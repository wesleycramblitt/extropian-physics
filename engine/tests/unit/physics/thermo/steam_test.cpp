// steam_test.cpp
// Engineering saturation model: anchor, inverses, monotonicity, enthalpies.

#include <exd/engine/physics/thermo/steam.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <initializer_list>

using namespace exd::engine::physics::thermo;

TEST_CASE("steam: saturation anchor and inverse round-trip")
{
    const SteamConstants c;
    // 100 °C ↔ 101.325 kPa by construction.
    CHECK(satured_pressure(c.T0) == doctest::Approx(c.p0).epsilon(1e-9));
    CHECK(saturation_temperature(c.p0) == doctest::Approx(c.T0).epsilon(1e-9));
    // Inverse of forward over the operating band.
    for (double T : {300.0, 340.0, 373.15, 420.0, 450.0})
    {
        const double p = satured_pressure(T);
        CHECK(saturation_temperature(p) == doctest::Approx(T).epsilon(1e-6));
    }
}

TEST_CASE("steam: monotonic saturation line")
{
    double prev = 0.0;
    for (double T = 280.0; T <= 480.0; T += 10.0)
    {
        const double p = satured_pressure(T);
        CHECK(p > prev);
        prev = p;
    }
    // Sanity ranges: ~1 kPa at 280 K, ~1.5 MPa at 470 K (engineering band).
    CHECK(satured_pressure(280.0) > 500.0);
    CHECK(satured_pressure(280.0) < 3000.0);
    CHECK(satured_pressure(470.0) > 1.0e6);
    CHECK(satured_pressure(470.0) < 2.0e6);
}

TEST_CASE("steam: enthalpies are energy-consistent")
{
    const SteamConstants c;
    for (double T : {323.15, 373.15, 423.15})
    {
        CHECK(h_g(T) - h_f(T) == doctest::Approx(h_fg(T)).epsilon(1e-9));
        CHECK(h_f(T) > 0.0);
    }
    // Dryness fraction blends liquid and vapor enthalpies.
    const double h_wet_val = h_wet(SteamConstants{}.T0, 0.5);
    CHECK(h_wet_val == doctest::Approx(0.5 * (h_f(SteamConstants{}.T0) + h_g(SteamConstants{}.T0))).epsilon(1e-9));
}

TEST_CASE("steam: saturated vapor density is positive and finite")
{
    const SteamConstants c;
    const double rho = rho_g(c.p0, c.T0, c);
    CHECK(rho == doctest::Approx(c.p0 / (c.r_steam * c.T0)).epsilon(1e-12));
    CHECK(rho > 0.55); // ~0.59 kg/m³ at 100 °C
    CHECK(rho < 0.65);
}
