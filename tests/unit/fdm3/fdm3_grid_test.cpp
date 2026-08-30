#include <doctest/doctest.h>
#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>
#include "../../../../src/fluid/fdm3/fdm3_internal.hpp"

#include <string>
#include <vector>

using namespace exd::physics::fluid::fdm3;

namespace {

// A valid all-around configuration for validation tests.
FDM3Config valid_config() {
    FDM3Config c;
    c.nx = 8;
    c.ny = 8;
    c.nz = 8;
    c.boundary_conditions = {
        {BoundaryFace::XMin, FDMBoundaryType::Inlet, 1.0, 0.0, 0.0},
        {BoundaryFace::XMax, FDMBoundaryType::Outlet},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::ZMax, FDMBoundaryType::Symmetry},
    };
    return c;
}

} // anonymous namespace

TEST_CASE("FDM3Grid allocation and indexing") {
    FDM3Grid g;
    g.allocate(10, 8, 12);
    g.dx = 0.1;
    g.dy = 0.125;
    g.dz = 1.0 / 12.0;

    CHECK(g.nx == 10);
    CHECK(g.ny == 8);
    CHECK(g.nz == 12);
    CHECK(g.sx() == 12);   // nx + 2
    CHECK(g.sy() == 10);   // ny + 2
    CHECK(g.idx(0, 0, 0) == 0);
    CHECK(g.idx(1, 0, 0) == 1);
    CHECK(g.idx(0, 1, 0) == 12);
    CHECK(g.idx(0, 0, 1) == 12 * 10);
    CHECK(g.idx(11, 9, 13) == g.total() - 1);

    const size_t expected = static_cast<size_t>(12) * 10 * 14;
    CHECK(g.u.size() == expected);
    CHECK(g.v.size() == expected);
    CHECK(g.w.size() == expected);
    CHECK(g.p.size() == expected);
    CHECK(g.u_old.size() == expected);
    CHECK(g.w_old.size() == expected);
    CHECK(g.p_prime.size() == expected);
    CHECK(g.rhs.size() == expected);
}

TEST_CASE("FDM3Grid initialization fills interior and zero ghosts") {
    FDM3Grid g;
    g.allocate(4, 4, 4);
    g.dx = 0.25;
    g.dy = 0.25;
    g.dz = 0.25;
    g.initialize(1.0, 0.5, -2.0, 101325.0);

    // Interior cells take the initial conditions.
    CHECK(g.u[g.idx(2, 2, 2)] == doctest::Approx(1.0));
    CHECK(g.v[g.idx(1, 3, 2)] == doctest::Approx(0.5));
    CHECK(g.w[g.idx(3, 3, 3)] == doctest::Approx(-2.0));
    CHECK(g.p[g.idx(4, 4, 4)] == doctest::Approx(101325.0));

    // Ghost cells stay zero.
    CHECK(g.u[g.idx(0, 2, 2)] == doctest::Approx(0.0));
    CHECK(g.v[g.idx(2, 0, 2)] == doctest::Approx(0.0));
    CHECK(g.w[g.idx(2, 2, 0)] == doctest::Approx(0.0));
    CHECK(g.p[g.idx(5, 2, 2)] == doctest::Approx(0.0));
    CHECK(g.u[g.idx(0, 0, 0)] == doctest::Approx(0.0));
}

TEST_CASE("initialize_grid uses config for dims, spacings and IC") {
    FDM3Config config;
    config.nx = 6;
    config.ny = 7;
    config.nz = 8;
    config.lx = 2.0;
    config.ly = 3.0;
    config.lz = 4.0;
    config.initial_u = 2.5;
    config.initial_v = -1.0;
    config.initial_w = 0.25;
    config.initial_p = 5.0;

    FDM3Grid g;
    initialize_grid(g, config);
    CHECK(g.nx == 6);
    CHECK(g.dx == doctest::Approx(2.0 / 6.0));
    CHECK(g.dy == doctest::Approx(3.0 / 7.0));
    CHECK(g.dz == doctest::Approx(4.0 / 8.0));
    CHECK(g.u[g.idx(3, 4, 5)] == doctest::Approx(2.5));
    CHECK(g.v[g.idx(1, 1, 1)] == doctest::Approx(-1.0));
    CHECK(g.w[g.idx(3, 4, 5)] == doctest::Approx(0.25));
    CHECK(g.p[g.idx(1, 2, 3)] == doctest::Approx(5.0));
    CHECK(g.u[g.idx(0, 4, 5)] == doctest::Approx(0.0));  // ghost stays zero
}

TEST_CASE("validation: valid config passes cleanly") {
    auto c = valid_config();
    std::string error;
    std::vector<std::string> warnings;
    CHECK(c.validate(error, warnings));
    CHECK(error.empty());
    CHECK(warnings.empty());
}

TEST_CASE("validation: grid too small is an error") {
    auto c = valid_config();
    c.nx = 3;
    SUBCASE("nx") {
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
    c = valid_config();
    c.ny = 2;
    SUBCASE("ny") {
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
    c = valid_config();
    c.nz = 1;
    SUBCASE("nz") {
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
}

TEST_CASE("validation: non-positive dimensions are errors") {
    auto c = valid_config();
    c.lx = 0.0;
    std::string error;
    std::vector<std::string> warnings;
    CHECK_FALSE(c.validate(error, warnings));
    CHECK_FALSE(error.empty());
}

TEST_CASE("validation: non-positive density, dt, max_steps are errors") {
    SUBCASE("rho") {
        auto c = valid_config();
        c.rho = 0.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
    SUBCASE("dt") {
        auto c = valid_config();
        c.dt = 0.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
    SUBCASE("max_steps") {
        auto c = valid_config();
        c.max_steps = 0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
}

TEST_CASE("validation: SOR omega must be in (1, 2)") {
    SUBCASE("omega = 1.0") {
        auto c = valid_config();
        c.sor_omega = 1.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
        CHECK_FALSE(error.empty());
    }
    SUBCASE("omega = 2.0") {
        auto c = valid_config();
        c.sor_omega = 2.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
    }
    SUBCASE("omega in range passes") {
        auto c = valid_config();
        c.sor_omega = 1.5;
        std::string error;
        std::vector<std::string> warnings;
        CHECK(c.validate(error, warnings));
    }
}

TEST_CASE("validation: relaxation factors must be in (0, 1]") {
    SUBCASE("velocity 0") {
        auto c = valid_config();
        c.velocity_under_relaxation = 0.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
    }
    SUBCASE("velocity > 1") {
        auto c = valid_config();
        c.velocity_under_relaxation = 1.5;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
    }
    SUBCASE("pressure > 1") {
        auto c = valid_config();
        c.pressure_under_relaxation = 2.0;
        std::string error;
        std::vector<std::string> warnings;
        CHECK_FALSE(c.validate(error, warnings));
    }
}

TEST_CASE("validation: at least one boundary condition required") {
    FDM3Config c;
    c.nx = 8;
    c.ny = 8;
    c.nz = 8;
    std::string error;
    std::vector<std::string> warnings;
    CHECK_FALSE(c.validate(error, warnings));
    CHECK_FALSE(error.empty());
}

TEST_CASE("validation: duplicate face boundary condition is an error") {
    auto c = valid_config();
    c.boundary_conditions.push_back(
        {BoundaryFace::XMin, FDMBoundaryType::Wall}
    );
    std::string error;
    std::vector<std::string> warnings;
    CHECK_FALSE(c.validate(error, warnings));
    CHECK_FALSE(error.empty());
}

TEST_CASE("validation: periodic pairing warnings") {
    SUBCASE("natural opposite pairing passes without warnings") {
        auto c = valid_config();
        c.boundary_conditions = {
            {BoundaryFace::XMin, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::XMax},
            {BoundaryFace::XMax, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::XMin},
            {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
            {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
            {BoundaryFace::ZMin, FDMBoundaryType::Symmetry},
            {BoundaryFace::ZMax, FDMBoundaryType::Symmetry},
        };
        std::string error;
        std::vector<std::string> warnings;
        CHECK(c.validate(error, warnings));
        CHECK(warnings.empty());
    }
    SUBCASE("non-opposite periodic pair emits a warning") {
        auto c = valid_config();
        c.boundary_conditions = {
            {BoundaryFace::XMin, FDMBoundaryType::Periodic, 0.0, 0.0, 0.0, 0.0, BoundaryFace::YMax},
            {BoundaryFace::XMax, FDMBoundaryType::Symmetry},
            {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
            {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
            {BoundaryFace::ZMin, FDMBoundaryType::Symmetry},
            {BoundaryFace::ZMax, FDMBoundaryType::Symmetry},
        };
        std::string error;
        std::vector<std::string> warnings;
        CHECK(c.validate(error, warnings));
        CHECK_FALSE(warnings.empty());
    }
}

TEST_CASE("natural_opposite maps faces correctly") {
    CHECK(natural_opposite(BoundaryFace::XMin) == BoundaryFace::XMax);
    CHECK(natural_opposite(BoundaryFace::XMax) == BoundaryFace::XMin);
    CHECK(natural_opposite(BoundaryFace::YMin) == BoundaryFace::YMax);
    CHECK(natural_opposite(BoundaryFace::YMax) == BoundaryFace::YMin);
    CHECK(natural_opposite(BoundaryFace::ZMin) == BoundaryFace::ZMax);
    CHECK(natural_opposite(BoundaryFace::ZMax) == BoundaryFace::ZMin);
}