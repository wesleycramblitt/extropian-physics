// fdm3_output_driver_test.cpp
// run_fdm3_simulation: field stamping at a cadence + scheduler path.

#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/physics/io/field_writer.hpp>
#include <exd/physics/io/output_policy.hpp>

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace exd::physics::fluid::fdm3;
using namespace exd::physics::io;
using exd::physics::ModelStatus;

namespace
{
FDM3Config small_config()
{
    FDM3Config cfg;
    cfg.nx = 8;
    cfg.ny = 6;
    cfg.nz = 8;
    cfg.lx = 1.6;
    cfg.ly = 1.2;
    cfg.lz = 1.6;
    cfg.dt = 0.005;
    cfg.max_steps = 60;
    cfg.boundary_conditions = {
        {BoundaryFace::ZMax, FDMBoundaryType::Inlet, 0.0, 0.0, -1.0},
        {BoundaryFace::ZMin, FDMBoundaryType::Outlet},
        {BoundaryFace::XMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::XMax, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMin, FDMBoundaryType::Symmetry},
        {BoundaryFace::YMax, FDMBoundaryType::Symmetry},
    };
    return cfg;
}

std::string temp_dir()
{
    auto base = std::filesystem::temp_directory_path();
    auto d = base / ("exd_fdm3_out_" + std::to_string(::getpid()));
    std::error_code ec;
    std::filesystem::remove_all(d, ec);
    std::filesystem::create_directories(d, ec);
    return d.string();
}
} // anonymous namespace

TEST_CASE("run_fdm3_simulation: stamps fields at the configured interval")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(writer);

    auto cfg = small_config();
    cfg.field_stamp_interval = 20; // stamps at steps 0, 20, 40 → 3 files

    auto r = run_fdm3_simulation(cfg, writer.get(), nullptr, &status);
    REQUIRE(r.valid);

    const std::string timeline = dir + "/timeline.txt";
    std::ifstream tl(timeline);
    REQUIRE(tl.good());
    std::string line;
    int stamps = 0;
    while (std::getline(tl, line)) ++stamps;
    CHECK(stamps == 3); // 0, 20, 40 (steps 0..59)
    CHECK(std::filesystem::exists(dir + "/step_00000000.fld"));
    CHECK(std::filesystem::exists(dir + "/step_00000040.fld"));
}

TEST_CASE("run_fdm3_simulation: scheduler controls the cadence")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(writer);

    auto cfg = small_config();
    // Deterministic injected clock: sim time IS the wall clock here
    // (0.005 s/step), so a 0.01 s cadence emits every 2nd step.
    OutputScheduler sched(OutputPolicy{0, 0.01}, 0.0);
    auto r = run_fdm3_simulation(cfg, writer.get(), &sched, &status);
    REQUIRE(r.valid);

    std::ifstream tl(dir + "/timeline.txt");
    REQUIRE(tl.good());
    std::string line;
    int stamps = 0;
    while (std::getline(tl, line)) ++stamps;
    CHECK(stamps >= 25); // ≈ 60 steps / 2 per cadence
    CHECK(stamps <= 31);
}

TEST_CASE("run_fdm3_simulation: no writer equals solve_fdm3 behavior")
{
    auto cfg = small_config();
    ModelStatus status;
    auto r = run_fdm3_simulation(cfg, nullptr, nullptr, &status);
    REQUIRE(r.valid);
    CHECK(r.steps_taken == cfg.max_steps);
    CHECK(r.field.u.size() == static_cast<size_t>(cfg.nx) * cfg.ny * cfg.nz);
}
