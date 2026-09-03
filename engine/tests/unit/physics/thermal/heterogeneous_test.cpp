// heterogeneous_test.cpp — the GENERAL heterogeneous steady-conduction class
// (W15).  The "chiplet board" is just one fixture of this class — the same
// config shape serves heat sinks, reactor walls, building envelopes, etc.
// Verification:
//   * a single heat source → energy balance (source == sink flux) + hot spot
//   * doubling the source → doubling the peak rise (linearity)
//   * a high-k spreader region lowers the peak
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/physics/thermal/heterogeneous.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::mesh;
using namespace exd::engine::physics::thermal;

namespace {
/// A general box-region conduction problem on a thin plate (the "chiplet
/// board" shape: two heat sources plus an optional spreader).
/// The region contract is W/m³; the DISCRETE integrated power is the sum
/// over the nodes the box covers (a surface-counting effect at the box
/// boundary is expected and accounted for in the fixture).
struct PlateFixture
{
    HeterogeneousConductionConfig cfg;
    double expected_power = 0.0;
};

PlateFixture plate_with_sources(double power1, double power2, bool with_spreader)
{
    PlateFixture fx;
    auto& cfg = fx.cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.02, 0.02, 0.001}, {41, 41, 3});
    cfg.base_conductivity = 20.0;
    // board BCs: the in-plane edges are the sinks, the top/bottom faces are
    // INSULATED (a thin board conducts in-plane — the honest chiplet case)
    for (int a = 0; a < 6; ++a)
        cfg.face_kind[static_cast<size_t>(a)] = ConductionFaceKind::FixedValue;
    cfg.face_kind[4] = ConductionFaceKind::Insulated;
    cfg.face_kind[5] = ConductionFaceKind::Insulated;
    const double vol = 0.004 * 0.004 * 0.001;
    SourceRegion s1;
    s1.center = {0.005, 0.005, 0.0};
    s1.half_extents = {0.002, 0.002, 0.001};
    s1.volumetric_source = power1 / vol;
    SourceRegion s2 = s1;
    s2.center = {0.015, 0.015, 0.0};
    s2.volumetric_source = power2 / vol;
    cfg.sources = {s1, s2};
    // the integrated power: count the nodes the box covers × q × cell_vol
    const double cell_vol = cfg.grid.cell_volume();
    auto count = [&](const SourceRegion& s) {
        int n = 0;
        const auto& g = cfg.grid;
        for (int k = 0; k < g.dims[2]; ++k)
            for (int j = 0; j < g.dims[1]; ++j)
                for (int i = 0; i < g.dims[0]; ++i)
                {
                    const double x = g.origin[0] + i * g.spacing[0];
                    const double y = g.origin[1] + j * g.spacing[1];
                    const double z = g.origin[2] + k * g.spacing[2];
                    if (std::fabs(x - s.center[0]) <= s.half_extents[0] &&
                        std::fabs(y - s.center[1]) <= s.half_extents[1] &&
                        std::fabs(z - s.center[2]) <= s.half_extents[2])
                        ++n;
                }
        return n;
    };
    fx.expected_power = (count(s1) * s1.volumetric_source +
                         count(s2) * s2.volumetric_source) * cell_vol;
    if (with_spreader)
    {
        ConductionRegion sp;
        sp.center = {0.01, 0.01, 0.0};
        sp.half_extents = {0.0075, 0.0075, 0.001};
        sp.conductivity = 400.0;               // copper
        cfg.materials = {sp};
    }
    return fx;
}
} // namespace

TEST_CASE("Heterogeneous conduction: energy balance and hot spot (chiplet fixture)")
{
    const auto fx = plate_with_sources(3.0, 3.0, false);
    const auto r = solve_heterogeneous_conduction(fx.cfg);
    REQUIRE(r.ok);
    CHECK(r.total_power == doctest::Approx(fx.expected_power).epsilon(1e-9));
    CHECK(r.sink_flux == doctest::Approx(r.total_power).epsilon(0.02));
    CHECK(r.peak_temperature > 300.0);
    // the hot spot sits at one of the sources
    const double d1 = std::fabs(r.peak_x - 0.015);
    const double d2 = std::fabs(r.peak_x - 0.005);
    const double dmin = std::min(d1, d2);
    CHECK(dmin < 0.004);
}

TEST_CASE("Heterogeneous conduction: peak rise scales linearly with power")
{
    const auto r1 = solve_heterogeneous_conduction(plate_with_sources(5.0, 5.0, false).cfg);
    REQUIRE(r1.ok);
    const auto r2 = solve_heterogeneous_conduction(plate_with_sources(10.0, 10.0, false).cfg);
    REQUIRE(r2.ok);
    const double rise1 = r1.peak_temperature - 300.0;
    const double rise2 = r2.peak_temperature - 300.0;
    CHECK(rise1 > 0.0);
    CHECK(rise2 == doctest::Approx(2.0 * rise1).epsilon(0.01));
}

TEST_CASE("Heterogeneous conduction: the spreader region lowers the peak")
{
    const auto fx = plate_with_sources(3.0, 3.0, false);
    const auto r_plain = solve_heterogeneous_conduction(fx.cfg);
    REQUIRE(r_plain.ok);
    const auto r_spread = solve_heterogeneous_conduction(plate_with_sources(3.0, 3.0, true).cfg);
    REQUIRE(r_spread.ok);
    CHECK(r_spread.peak_temperature < r_plain.peak_temperature);
    CHECK(r_spread.sink_flux == doctest::Approx(fx.expected_power).epsilon(0.02));
}

TEST_CASE("Heterogeneous conduction: data-driven fields (the CAD contract)")
{
    // the user CAN skip regions entirely and hand over per-node k and q
    HeterogeneousConductionConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.01, 0.01, 0.001}, {21, 21, 3});
    cfg.base_conductivity = 10.0;
    cfg.conductivity_field.assign(cfg.grid.node_count(), 10.0);
    cfg.source_field.assign(cfg.grid.node_count(), 0.0);
    // heated slab in the middle, copper everywhere else
    const size_t nx = static_cast<size_t>(cfg.grid.dims[0]);
    const size_t ny = static_cast<size_t>(cfg.grid.dims[1]);
    for (size_t k = 0; k < static_cast<size_t>(cfg.grid.dims[2]); ++k)
        for (size_t j = 0; j < ny; ++j)
            for (size_t i = 0; i < nx; ++i)
            {
                const size_t id = i + nx * (j + ny * k);
                if (i >= 9 && i <= 11 && j >= 9 && j <= 11)
                {
                    cfg.conductivity_field[id] = 400.0;
                    cfg.source_field[id] = 1e9;   // W/m³
                }
            }
    const auto r = solve_heterogeneous_conduction(cfg);
    REQUIRE(r.ok);
    CHECK(r.total_power > 0.0);
    CHECK(r.sink_flux == doctest::Approx(r.total_power).epsilon(0.05));
    CHECK(r.peak_temperature > 300.0);
}

TEST_CASE("Heterogeneous conduction: insulated faces admit a clean 1D profile")
{
    // two fixed-temperature x faces, insulated y/z: the exact linear bridge
    HeterogeneousConductionConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {1, 0.1, 0.1}, {41, 3, 3});
    cfg.base_conductivity = 50.0;
    cfg.face_value = {300, 400, 0, 0, 0, 0};
    for (int a = 0; a < 6; ++a)
        cfg.face_kind[static_cast<size_t>(a)] = (a < 2)
            ? ConductionFaceKind::FixedValue
            : ConductionFaceKind::Insulated;
    const auto r = solve_heterogeneous_conduction(cfg);
    REQUIRE(r.ok);
    double err = 0.0;
    for (int k = 0; k < 3; ++k)
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 41; ++i)
            {
                const double exact = 300.0 + 100.0 * static_cast<double>(i) / 40.0;
                const size_t id = static_cast<size_t>(i) + 41ull *
                    (static_cast<size_t>(j) + 3ull * k);
                err = std::max(err, std::fabs(r.temperature.values[id] - exact));
            }
    CHECK(err < 1e-6);                          // exact discrete linear bridge
}
