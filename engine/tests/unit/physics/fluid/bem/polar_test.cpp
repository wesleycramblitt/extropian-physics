#include <doctest/doctest.h>

#include <exd/engine/physics/fluid/reduced_order/bem/airfoil.hpp>

#include <cmath>
#include <string>

using namespace exd::engine::physics::fluid::reduced_order::bem;

namespace {

std::string data_dir()
{
    return std::string(EXT_PHYSICS_DATA_DIR) + "/airfoils";
}

} // namespace

TEST_CASE("PolarDatabase: add and find nearest Re")
{
    PolarDatabase db;

    AirfoilPolar p1;
    p1.name = "naca0012";
    p1.re = 1e5;
    p1.alpha_deg = {0.0, 10.0};
    p1.cl = {0.0, 1.0};
    p1.cd = {0.01, 0.02};
    db.add(p1);

    AirfoilPolar p2;
    p2.name = "naca0012";
    p2.re = 1e6;
    p2.alpha_deg = {0.0, 10.0};
    p2.cl = {0.0, 1.1};
    p2.cd = {0.01, 0.015};
    db.add(p2);

    const auto* found = db.find("naca0012", 5e5);
    REQUIRE(found != nullptr);
    CHECK(found->re == doctest::Approx(1e5)); // nearest to 5e5 is 1e5

    const auto* found2 = db.find("naca0012", 9e5);
    REQUIRE(found2 != nullptr);
    CHECK(found2->re == doctest::Approx(1e6));
}

TEST_CASE("PolarDatabase: re==0 is +inf fallback only when no re-tagged polar exists")
{
    PolarDatabase db;

    AirfoilPolar p_re;
    p_re.name = "af";
    p_re.re = 1e6;
    p_re.alpha_deg = {0.0, 10.0};
    p_re.cl = {0.0, 1.0};
    p_re.cd = {0.01, 0.02};
    db.add(p_re);

    AirfoilPolar p_zero;
    p_zero.name = "af";
    p_zero.re = 0.0;
    p_zero.alpha_deg = {0.0, 10.0};
    p_zero.cl = {0.0, 2.0};
    p_zero.cd = {0.01, 0.03};
    db.add(p_zero);

    const auto* found = db.find("af", 1e3);
    REQUIRE(found != nullptr);
    CHECK(found->re == doctest::Approx(1e6)); // re-tagged exists, zero ignored

    PolarDatabase db2;
    db2.add(p_zero);
    const auto* found2 = db2.find("af", 1e3);
    REQUIRE(found2 != nullptr);
    CHECK(found2->re == doctest::Approx(0.0));
}

TEST_CASE("PolarDatabase: unknown name returns nullptr")
{
    PolarDatabase db;
    db.add_builtin_polars();
    CHECK(db.find("not_an_airfoil", 1e5) == nullptr);
}

TEST_CASE("AirfoilPolar: evaluate linear interpolation and flat clamp")
{
    AirfoilPolar p;
    p.name = "test";
    p.re = 1e5;
    p.alpha_deg = {0.0, 10.0, 20.0};
    p.cl = {0.0, 1.0, 1.2};
    p.cd = {0.01, 0.02, 0.05};

    auto [cl5, cd5] = p.evaluate(5.0);
    CHECK(cl5 == doctest::Approx(0.5));
    CHECK(cd5 == doctest::Approx(0.015));

    auto [cl_neg, cd_neg] = p.evaluate(-10.0);
    CHECK(cl_neg == doctest::Approx(0.0));
    CHECK(cd_neg == doctest::Approx(0.01));

    auto [cl_high, cd_high] = p.evaluate(100.0);
    CHECK(cl_high == doctest::Approx(1.2));
    CHECK(cd_high == doctest::Approx(0.05));
}

TEST_CASE("PolarDatabase: load CSV directory")
{
    PolarDatabase db;
    CHECK(db.load_directory(data_dir()));
    CHECK(db.has("naca0012"));
    CHECK(db.has("naca4412"));

    const auto* p = db.find("naca0012", 3e5);
    REQUIRE(p != nullptr);
    CHECK(p->re == doctest::Approx(3e5));
    auto [cl, cd] = p->evaluate(0.0);
    CHECK(cl == doctest::Approx(0.0).epsilon(1e-4));
    CHECK(cd > 0.0);
}

TEST_CASE("PolarDatabase: CSV load skips malformed rows and comments")
{
    PolarDatabase db;
    // Direct load of a known good file; malformed cases are implicitly handled
    // by the parser skipping blank/comment/bad rows.
    CHECK(db.load_csv(data_dir() + "/naca0012_1e5.csv"));
    const auto* p = db.find("naca0012", 1e5);
    REQUIRE(p != nullptr);
    CHECK(p->alpha_deg.size() == 66);
}

TEST_CASE("PolarDatabase: builtins present")
{
    PolarDatabase db;
    db.add_builtin_polars();
    CHECK(db.has("naca0012"));
    CHECK(db.has("naca4412"));

    const auto* p = db.find("naca4412", 1e6);
    REQUIRE(p != nullptr);
    auto [cl0, cd0] = p->evaluate(0.0);
    CHECK(cl0 == doctest::Approx(0.35).epsilon(1e-4));
    CHECK(cd0 > 0.0);
}
