// steady_conduction_preset_test.cpp — the CLASS-level heat-conduction preset:
// the same preset entry runs every problem in the class (here: a chiplet-like
// board as one fixture).
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/presets/thermal/steady_conduction.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::mesh;
using namespace exd::engine::presets::thermal;

TEST_CASE("steady_conduction preset: class-level entry, chiplet-like fixture")
{
    SteadyConductionConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.02, 0.02, 0.001}, {41, 41, 3});
    cfg.base_conductivity = 20.0;
    for (int a = 0; a < 6; ++a)
        cfg.face_kind[static_cast<size_t>(a)] = exd::engine::physics::thermal::ConductionFaceKind::FixedValue;
    cfg.face_kind[4] = exd::engine::physics::thermal::ConductionFaceKind::Insulated;
    cfg.face_kind[5] = exd::engine::physics::thermal::ConductionFaceKind::Insulated;
    exd::engine::physics::thermal::SourceRegion chip;
    chip.center = {0.01, 0.01, 0.0};
    chip.half_extents = {0.002, 0.002, 0.001};
    chip.volumetric_source = 5.0 / (0.004 * 0.004 * 0.001);
    cfg.sources = {chip};
    const double cell_vol = cfg.grid.cell_volume();
    double expected = 0.0;
    for (int k = 0; k < cfg.grid.dims[2]; ++k)
        for (int j = 0; j < cfg.grid.dims[1]; ++j)
            for (int i = 0; i < cfg.grid.dims[0]; ++i)
            {
                const double x = cfg.grid.origin[0] + i * cfg.grid.spacing[0];
                const double y = cfg.grid.origin[1] + j * cfg.grid.spacing[1];
                const double z = cfg.grid.origin[2] + k * cfg.grid.spacing[2];
                if (std::fabs(x - 0.01) <= 0.002 && std::fabs(y - 0.01) <= 0.002 &&
                    std::fabs(z - 0.0) <= 0.001)
                    expected += chip.volumetric_source * cell_vol;
            }

    const auto r = run_steady_conduction(cfg);
    REQUIRE(r.ok);
    CHECK(r.total_power == doctest::Approx(expected).epsilon(1e-9));
    CHECK(r.sink_flux == doctest::Approx(expected).epsilon(0.02));
    CHECK(r.peak_temperature > 300.0);
}
