#include <doctest/doctest.h>
#include <exd/physics/fluid/fdm3/fdm3_config.hpp>
#include <exd/physics/fluid/fdm3/fdm3_result.hpp>
#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>

#include <cmath>
#include <vector>

using namespace exd::physics;
using namespace exd::physics::fluid::fdm3;

namespace {

// Fixed actuator disk: a decelerating axial body-force slab inside an
// otherwise uniform inflow.  Momentum theory with induction a=0.2 predicts a
// far-wake deficit factor (1-2a) = 0.6.
FDM3Config disk_config() {
    FDM3Config c;
    c.nx = 32; c.ny = 16; c.nz = 32;
    c.lx = 1.0; c.ly = 0.5; c.lz = 1.0;
    c.rho = 1.225;
    c.mu = 1e-5;
    c.dt = 0.005;
    c.max_steps = 450;
    c.time_integration = TimeIntegration::Heun;
    c.advection_scheme = AdvectionScheme::Upwind;
    c.initial_u = 1.0;
    c.pressure_max_iterations = 100;
    c.pressure_tolerance = 1e-5;
    c.sor_omega = 1.5;
    c.velocity_under_relaxation = 0.7;
    c.pressure_under_relaxation = 0.3;
    c.convergence_window = 1000;   // never triggers early
    c.convergence_tolerance = 1e-12;
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

// Fill a decelerating axial body force over a disk slab.  The force is the
// acceleration experienced by the FLUID (negative in x decelerates it).
// f0 is chosen so that the slab carries the momentum-theory thrust
// T = 2*rho*A*u^2*a*(1-a):  f0*t = 2*u^2*a*(1-a) with t the slab thickness.
std::vector<double> disk_force_x(const FDM3Config& config, double u, double a,
                                 int i_mid, int radius_cells, int half_thick) {
    const double t = static_cast<double>(2 * half_thick + 1) * config.dx();
    const double f0 = 2.0 * u * u * a * (1.0 - a) / t;
    const double yc = (config.ny - 1) / 2.0;       // mid j in cell units
    const double zc = (config.nz - 1) / 2.0;       // mid k in cell units
    const double r2 = static_cast<double>(radius_cells * radius_cells);

    std::vector<double> fx(static_cast<size_t>(config.nx) * config.ny * config.nz, 0.0);
    for (int k = 0; k < config.nz; ++k)
        for (int j = 0; j < config.ny; ++j)
            for (int i = 0; i < config.nx; ++i) {
                const double dy = static_cast<double>(j) - yc;
                const double dz = static_cast<double>(k) - zc;
                if (dy * dy + dz * dz <= r2 && std::abs(i - i_mid) <= half_thick)
                    fx[static_cast<size_t>(i) + config.nx * (static_cast<size_t>(j) +
                        static_cast<size_t>(config.ny) * k)] = -f0;
            }
    return fx;
}

double centerline_u(const FDM3Result& result, int i) {
    return result.field.u[result.field.index(i, result.field.ny / 2, result.field.nz / 2)];
}

} // anonymous namespace

TEST_CASE("fdm3: decelerating actuator disk reduces downstream velocity") {
    auto config = disk_config();

    FDM3Solver solver;
    ModelStatus status;
    REQUIRE(solver.initialize(config, status));
    REQUIRE(status.ok);

    // Disk centered at i=16 (x ~ 0.5), radius 4 cells, 3 cells thick.
    const int i_mid = 16;
    const int radius = 4;
    const int half_thick = 1;
    const double a = 0.2;
    auto fx = disk_force_x(config, 1.0, a, i_mid, radius, half_thick);
    std::vector<double> fy(fx.size(), 0.0);
    std::vector<double> fz(fx.size(), 0.0);

    REQUIRE(solver.set_body_force(fx, fy, fz, status));
    REQUIRE(status.ok);

    for (int s = 0; s < config.max_steps; ++s) {
        REQUIRE(solver.step(config.dt, status));
        REQUIRE(status.ok);
    }

    const auto& field = solver.field();

    // Sample the centerline before and after the disk.
    const double u_us = field.u[field.index(10, field.ny / 2, field.nz / 2)];
    const double u_disk = field.u[field.index(i_mid, field.ny / 2, field.nz / 2)];
    const double u_ds = field.u[field.index(20, field.ny / 2, field.nz / 2)];
    const double u_far = field.u[field.index(26, field.ny / 2, field.nz / 2)];

    MESSAGE("u_us=", u_us, " u_disk=", u_disk, " u_ds=", u_ds, " u_far=", u_far);

    // The decelerating disk must reduce the streamwise velocity.
    CHECK(u_ds < u_us);
    CHECK(u_disk < u_us);
    CHECK(u_far < u_us);
    CHECK(u_us > 0.5);          // upstream is essentially the freestream
    CHECK(u_ds > 0.0);          // no flow reversal on the centerline

    // Momentum-theory far-wake deficit factor (1-2a) = 0.6, with a generous
    // engineering tolerance.  Blockage and the collocated projector shift the
    // measured value, so this is deliberately loose.
    const double deficit = u_ds / u_us;
    CHECK(deficit < 0.95);
    CHECK(deficit > 0.35);
    CHECK(std::abs(deficit - 0.6) < 0.30);
}

TEST_CASE("fdm3: body force span validation") {
    auto config = disk_config();
    FDM3Solver solver;
    ModelStatus status;
    REQUIRE(solver.initialize(config, status));
    REQUIRE(status.ok);

    // Clearing the force with empty spans is allowed.
    std::vector<double> empty;
    REQUIRE(solver.set_body_force(empty, empty, empty, status));
    REQUIRE(status.ok);

    // Wrong-sized arrays are rejected.
    std::vector<double> short_vec(3, 0.0);
    CHECK_FALSE(solver.set_body_force(short_vec, short_vec, short_vec, status));
    CHECK_FALSE(status.ok);
}