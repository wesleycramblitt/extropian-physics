// units_test.cpp — spec §25: dimensional units and validation.
#include <exd/engine/core/units.hpp>
#include <doctest/doctest.h>

using namespace exd::engine::core;

TEST_CASE("Units: SI base exponents and derived quantities")
{
    CHECK(units::dimensionless.is_dimensionless());
    CHECK(!units::meter.is_dimensionless());
    CHECK(units::newton == units::kilogram * units::meter / units::second.pow(2));
    CHECK(units::pascal == units::newton / units::meter.pow(2));
    CHECK(units::joule == units::newton * units::meter);
    CHECK(units::watt == units::joule / units::second);
    CHECK(units::velocity == units::meter / units::second);
    CHECK(units::force == units::newton);
    CHECK(units::torque == units::newton * units::meter);
    CHECK(units::dynamic_viscosity == units::pascal * units::second);
    CHECK(units::kinematic_viscosity == units::meter.pow(2) / units::second);
}

TEST_CASE("Units: dimensional consistency and rejection")
{
    // force = mass × acceleration is dimensionally valid
    CHECK(units_compatible(units::force, units::kilogram * units::acceleration));
    // temperature → force is NOT valid without a physical transformation
    CHECK(!units_compatible(units::kelvin, units::force));
    CHECK(!units_compatible(units::kelvin, units::newton / units::meter.pow(3)));
}

TEST_CASE("Units: to_string is sane")
{
    // emission order follows the stored exponent order: m, kg, s, A, K, mol, cd
    CHECK(to_string(units::pascal) == "m⁻¹·kg·s⁻²");
    CHECK(to_string(units::dimensionless) == "—");
    CHECK(to_string(units::watt) == "m²·kg·s⁻³");
}
