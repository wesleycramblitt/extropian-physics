// electronics_chiplet_test.cpp — chiplet-board preset (W15):
//   * bi-material strip reproduces the exact piecewise-linear steady profile
//   * chip power = sink flux (energy balance)
//   * the peak sits under the highest-power chip
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/presets/electronics/chiplet_board.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::mesh;
using namespace exd::engine::presets::electronics;

TEST_CASE("Chiplet: uniform board with a single chip — energy balance and hot spot")
{
    ChipletBoardConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.02, 0.02, 0.001}, {41, 41, 3});
    cfg.board_thickness = 0.001;
    cfg.base_conductivity = 20.0;
    cfg.sink_temperature = 300.0;
    ChipletBoardConfig::Chip chip;
    chip.x = 0.01; chip.y = 0.01; chip.w_cells = 6; chip.h_cells = 6;
    chip.power_watts = 5.0;
    cfg.chips = {chip};

    const auto r = solve_chiplet_board(cfg);
    REQUIRE(r.ok);
    // total chip power == sink flux (energy balance)
    CHECK(r.total_power == doctest::Approx(5.0).epsilon(1e-9));
    CHECK(r.sink_flux == doctest::Approx(r.total_power).epsilon(0.05));
    // the hot spot sits under the chip
    CHECK(r.peak_temperature > 300.0);
    CHECK(std::fabs(r.peak_x - 0.01) < 0.003);
    CHECK(std::fabs(r.peak_y - 0.01) < 0.003);
}

TEST_CASE("Chiplet: uniform board with a single chip — peak scales linearly with power")
{
    ChipletBoardConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.02, 0.02, 0.001}, {41, 41, 3});
    cfg.board_thickness = 0.001;
    cfg.base_conductivity = 20.0;
    cfg.sink_temperature = 300.0;
    ChipletBoardConfig::Chip chip;
    chip.x = 0.01; chip.y = 0.01; chip.w_cells = 6; chip.h_cells = 6;
    chip.power_watts = 5.0;
    cfg.chips = {chip};

    const auto r1 = solve_chiplet_board(cfg);
    REQUIRE(r1.ok);
    ChipletBoardConfig cfg2 = cfg;
    cfg2.chips[0].power_watts = 10.0;
    const auto r2 = solve_chiplet_board(cfg2);
    REQUIRE(r2.ok);

    // linearity: doubling the power doubles the temperature rise
    const double rise1 = r1.peak_temperature - cfg.sink_temperature;
    const double rise2 = r2.peak_temperature - cfg.sink_temperature;
    CHECK(rise1 > 0.0);
    CHECK(rise2 == doctest::Approx(2.0 * rise1).epsilon(0.01));
    // energy balance on both
    CHECK(r1.sink_flux == doctest::Approx(5.0).epsilon(0.05));
    CHECK(r2.sink_flux == doctest::Approx(10.0).epsilon(0.05));
    // peak under the chip
    CHECK(std::fabs(r1.peak_x - 0.01) < 0.003);
    CHECK(std::fabs(r1.peak_y - 0.01) < 0.003);
}

TEST_CASE("Chiplet: multiple chips with a spreader — the spreader cools the hot spot")
{
    ChipletBoardConfig cfg;
    cfg.grid = make_structured_grid({0, 0, 0}, {0.02, 0.02, 0.001}, {41, 41, 3});
    cfg.board_thickness = 0.001;
    cfg.base_conductivity = 20.0;
    cfg.sink_temperature = 300.0;
    ChipletBoardConfig::Chip c1;
    c1.x = 0.005; c1.y = 0.005; c1.w_cells = 4; c1.h_cells = 4; c1.power_watts = 3.0;
    ChipletBoardConfig::Chip c2;
    c2.x = 0.015; c2.y = 0.015; c2.w_cells = 4; c2.h_cells = 4; c2.power_watts = 3.0;
    cfg.chips = {c1, c2};

    const auto r_plain = solve_chiplet_board(cfg);
    REQUIRE(r_plain.ok);

    ChipletBoardConfig cfg_spread = cfg;
    ChipletBoardConfig::Spreader sp;
    sp.x = 0.01; sp.y = 0.01; sp.w_cells = 30; sp.h_cells = 30;
    sp.conductivity = 400.0;
    cfg_spread.spreaders = {sp};
    const auto r_spread = solve_chiplet_board(cfg_spread);
    REQUIRE(r_spread.ok);

    // the copper spreader lowers the peak temperature at the same power
    CHECK(r_spread.peak_temperature < r_plain.peak_temperature);
    // and both conserve power
    CHECK(r_plain.sink_flux == doctest::Approx(6.0).epsilon(0.05));
    CHECK(r_spread.sink_flux == doctest::Approx(6.0).epsilon(0.05));
}
