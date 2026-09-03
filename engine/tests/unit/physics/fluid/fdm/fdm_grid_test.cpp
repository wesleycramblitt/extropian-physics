#include <doctest/doctest.h>
#include <exd/engine/physics/fluid/fdm/fdm_config.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_result.hpp>
#include "../../../../../src/physics/fluid/fdm/fdm_internal.hpp"

using namespace exd::engine::physics::fluid::fdm;

TEST_CASE("FDMGrid allocation and indexing") {
    FDMGrid g;
    g.allocate(10, 8);
    g.dx = 0.1;
    g.dy = 0.125;

    CHECK(g.nx == 10);
    CHECK(g.ny == 8);
    CHECK(g.stride() == 12);
    CHECK(g.idx(0, 0) == 0);
    CHECK(g.idx(1, 0) == 1);
    CHECK(g.idx(0, 1) == 12);

    size_t expected = 12 * 10;
    CHECK(g.u.size() == expected);
    CHECK(g.v.size() == expected);
    CHECK(g.p.size() == expected);
}

TEST_CASE("FDMGrid initialization") {
    FDMGrid g;
    g.allocate(4, 4);
    g.dx = 0.25;
    g.dy = 0.25;
    g.initialize(1.0, 0.5, 101325.0);

    CHECK(g.u[g.idx(0, 0)] == doctest::Approx(1.0));
    CHECK(g.v[g.idx(2, 2)] == doctest::Approx(0.5));
    CHECK(g.p[g.idx(3, 3)] == doctest::Approx(101325.0));
}