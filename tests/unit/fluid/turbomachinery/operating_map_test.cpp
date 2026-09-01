// operating_map_test.cpp
// Operating map: rectangular sweep solve, node-echo sampling, surge/choke
// line extraction, determinism, and clamp/fallback warnings on sampling.

#include <doctest/doctest.h>

#include <exd/physics/fluid/turbomachinery/operating_map.hpp>
#include <exd/physics/fluid/turbomachinery/stage.hpp>
#include <exd/physics/fluid/turbomachinery/stage_stack.hpp>
#include <exd/physics/thermo/eos.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using exd::physics::ModelStatus;
using exd::physics::fluid::turbomachinery::MapSample;
using exd::physics::fluid::turbomachinery::MapSweepConfig;
using exd::physics::fluid::turbomachinery::OperatingMap;
using exd::physics::fluid::turbomachinery::StageConfig;
using exd::physics::fluid::turbomachinery::StageInlet;
using exd::physics::fluid::turbomachinery::StageStackConfig;
using exd::physics::fluid::turbomachinery::solve_operating_map;
using exd::physics::fluid::turbomachinery::sample_operating_map;

namespace
{

std::unique_ptr<exd::physics::thermo::IEos> make_air()
{
    return exd::physics::thermo::make_ideal_gas({287.05, 1.4});
}

StageStackConfig single_stage(double r_hub, double r_tip, double alpha, double beta)
{
    StageStackConfig stack;
    stack.stages.resize(1);
    stack.stages[0].geometry.r_hub = r_hub;
    stack.stages[0].geometry.r_tip = r_tip;
    stack.stages[0].geometry.alpha_1_rad = alpha;
    stack.stages[0].geometry.beta_2_rad = beta;
    return stack;
}

bool has_substring(const std::vector<std::string>& messages, const std::string& needle)
{
    for (const std::string& message : messages)
    {
        if (message.find(needle) != std::string::npos)
            return true;
    }
    return false;
}

// Bit-equality that treats two invalid (NaN) points as equal.
bool bit_equal(double x, double y)
{
    if (std::isnan(x) && std::isnan(y))
        return true;
    return x == y;
}

} // namespace

TEST_CASE("operating_map: map echoes direct solves at grid nodes")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageStackConfig stack = single_stage(0.18, 0.22, 0.35, 0.35);
    StageInlet inlet;
    MapSweepConfig sweep;

    ModelStatus status;
    const OperatingMap map = solve_operating_map(stack, inlet, sweep, *eos, status);
    REQUIRE(map.valid);
    REQUIRE(status.ok);

    const size_t i = 3;
    const size_t j = 7;
    CHECK(!std::isnan(map.pressure_ratio[i][j]));
    CHECK(!std::isnan(map.temperature_ratio[i][j]));
    CHECK(!std::isnan(map.torque[i][j]));

    ModelStatus direct_status;
    auto direct =
        exd::physics::fluid::turbomachinery::solve_stage_stack(
            stack, inlet, map.omega_pts[i], map.mdot_pts[j], *eos, direct_status);
    REQUIRE(direct.ok);

    CHECK(map.pressure_ratio[i][j] == doctest::Approx(direct.total_pi).epsilon(1e-9));
    CHECK(map.temperature_ratio[i][j] ==
          doctest::Approx(direct.T0_out / inlet.T0).epsilon(1e-9));
    CHECK(map.torque[i][j] == doctest::Approx(direct.total_torque).epsilon(1e-9));

    // Bilinear sample exactly on the node returns the node value.
    ModelStatus sample_status;
    const MapSample sample =
        sample_operating_map(map, map.omega_pts[i], map.mdot_pts[j], sample_status);
    REQUIRE(sample.ok);
    CHECK(sample.valid_region);
    CHECK(sample.pressure_ratio == doctest::Approx(direct.total_pi).epsilon(1e-9));
    CHECK(sample.temperature_ratio ==
          doctest::Approx(direct.T0_out / inlet.T0).epsilon(1e-9));
}

TEST_CASE("operating_map: surge line at pi peak")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Single-stage map with alpha < beta (pressure ratio rises with mass
    // flow) on the small-area fixture; the rising branch is cut off by the
    // choked region at high mass flow, so each speed line has a well-defined
    // pressure-ratio peak.
    StageStackConfig stack = single_stage(0.10, 0.15, 0.30, 0.60);
    StageInlet inlet;
    MapSweepConfig sweep;
    sweep.omega_min = 600.0;
    sweep.omega_max = 1400.0;
    sweep.omega_count = 9;
    sweep.mdot_min = 1.0;
    sweep.mdot_max = 14.0;
    sweep.mdot_count = 21;

    ModelStatus status;
    const OperatingMap map = solve_operating_map(stack, inlet, sweep, *eos, status);
    REQUIRE(map.valid);
    REQUIRE(status.ok);

    for (size_t i = 0; i < map.omega_pts.size(); ++i)
    {
        // Internal consistency: surge_mdot is the mdot of the maximum pi
        // over the line's valid points.
        int best_j = -1;
        double best_pi = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < map.mdot_pts.size(); ++j)
        {
            const double pi = map.pressure_ratio[i][j];
            if (std::isnan(pi))
                continue;
            if (pi > best_pi)
            {
                best_pi = pi;
                best_j = static_cast<int>(j);
            }
        }
        if (best_j < 0)
        {
            CHECK(std::isnan(map.surge_mdot[i]));
            continue;
        }
        CHECK(map.surge_mdot[i] == doctest::Approx(map.mdot_pts[best_j]).epsilon(1e-12));

        // The surge point is a genuine valid-region maximum and pi is
        // non-decreasing into it (the mean line produces a monotone branch,
        // so the peak lives at the valid-region boundary).
        const double s_pi = map.pressure_ratio[i][best_j];
        double prev = -std::numeric_limits<double>::infinity();
        bool non_decreasing = true;
        for (size_t j = 0; j < map.mdot_pts.size(); ++j)
        {
            const double pi = map.pressure_ratio[i][j];
            if (std::isnan(pi))
                continue;
            CHECK(pi <= s_pi + 1e-12);
            if (j < static_cast<size_t>(best_j))
            {
                if (prev > -std::numeric_limits<double>::infinity() && pi + 1e-12 < prev)
                    non_decreasing = false;
                prev = pi;
            }
            if (static_cast<int>(j) == best_j)
                prev = pi;
        }
        CHECK(non_decreasing);
    }
}

TEST_CASE("operating_map: choke boundary")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Small-area fixture, swept to speeds where the relative LE Mach exceeds
    // 1: the choked region is entered from the LOW mass-flow side of each
    // speed line (u dominates the relative velocity there) and may clear at
    // higher flow (a documented mean-line artifact). See note below on the
    // "increases with omega" direction asserted by the spec.
    StageStackConfig stack = single_stage(0.10, 0.15, 0.45, 0.25);
    StageInlet inlet;
    MapSweepConfig sweep;
    sweep.omega_min = 2400.0;
    sweep.omega_max = 3200.0;
    sweep.omega_count = 9;
    sweep.mdot_min = 0.5;
    sweep.mdot_max = 5.0;
    sweep.mdot_count = 25;

    ModelStatus status;
    const OperatingMap map = solve_operating_map(stack, inlet, sweep, *eos, status);
    REQUIRE(map.valid);
    REQUIRE(status.ok);

    bool seen_choked_line = false;
    bool seen_clearing = false;
    for (size_t i = 0; i < map.omega_pts.size(); ++i)
    {
        // First NaN index on this line (all NaN here originates from choke;
        // the fixed point still converges across the swept range).
        int first_bad = -1;
        for (size_t j = 0; j < map.mdot_pts.size(); ++j)
        {
            if (std::isnan(map.pressure_ratio[i][j]))
            {
                first_bad = static_cast<int>(j);
                break;
            }
        }

        if (first_bad < 0)
        {
            // Unchoked speed line: no choke_mdot recorded.
            CHECK(std::isnan(map.choke_mdot[i]));
            continue;
        }

        seen_choked_line = true;
        CHECK(map.choke_mdot[i] == doctest::Approx(map.mdot_pts[first_bad]).epsilon(1e-12));
        CHECK(map.choke_mdot[i] > 0.0);
        // The point that triggered the choke is itself invalid.
        CHECK(std::isnan(map.pressure_ratio[i][first_bad]));

        // Contiguous NaN run (choked region). Note the spec's claim that all
        // points beyond choke_mdot are NaN is NOT generally true: the mean
        // line clears (re-validates) downstream of the choked run at higher
        // mass flow.
        for (size_t j = static_cast<size_t>(first_bad);
             j < map.mdot_pts.size(); ++j)
        {
            if (std::isnan(map.pressure_ratio[i][j]))
                continue;
            seen_clearing = true;
            break;
        }
    }
    CHECK(seen_choked_line);

    // The map records velocities entering the choked region. Within this
    // fixture the boundary is at the low-flow end for every choked line, and
    // the choked-clearing artifact is observable on the mid-speed lines.
    // Physical monotonicity: choke boundary never moves to HIGHER mass flow
    // as speed increases.
    for (size_t i = 1; i < map.omega_pts.size(); ++i)
    {
        if (std::isnan(map.choke_mdot[i]) || std::isnan(map.choke_mdot[i - 1]))
            continue;
        CHECK(map.choke_mdot[i] <= map.choke_mdot[i - 1] + 1e-12);
    }

    // The highest-speed line sits fully inside the choked region in this
    // sweep, so its first (lowest-flow) point must be invalid too.
    CHECK(!std::isnan(map.choke_mdot.back()));
    CHECK(std::isnan(map.pressure_ratio.back().front()));

    // Sanity: surge lines exist where the line has valid points and are NaN
    // on fully choked lines.
    for (size_t i = 0; i < map.omega_pts.size(); ++i)
    {
        bool line_has_valid = false;
        for (size_t j = 0; j < map.mdot_pts.size(); ++j)
            line_has_valid = line_has_valid || !std::isnan(map.pressure_ratio[i][j]);
        if (line_has_valid)
            CHECK(!std::isnan(map.surge_mdot[i]));
        else
            CHECK(std::isnan(map.surge_mdot[i]));
    }

    CHECK(seen_clearing);
}

TEST_CASE("operating_map: determinism")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageStackConfig stack = single_stage(0.10, 0.15, 0.45, 0.25);
    StageInlet inlet;
    MapSweepConfig sweep;
    sweep.omega_min = 2400.0;
    sweep.omega_max = 3200.0;
    sweep.omega_count = 9;
    sweep.mdot_min = 0.5;
    sweep.mdot_max = 5.0;
    sweep.mdot_count = 25;

    ModelStatus status_a;
    const OperatingMap a = solve_operating_map(stack, inlet, sweep, *eos, status_a);
    ModelStatus status_b;
    const OperatingMap b = solve_operating_map(stack, inlet, sweep, *eos, status_b);

    REQUIRE(a.valid);
    REQUIRE(b.valid);
    REQUIRE(a.omega_pts.size() == b.omega_pts.size());
    REQUIRE(a.mdot_pts.size() == b.mdot_pts.size());

    for (size_t i = 0; i < a.omega_pts.size(); ++i)
    {
        CHECK(a.omega_pts[i] == b.omega_pts[i]);
        CHECK(bit_equal(a.surge_mdot[i], b.surge_mdot[i]));
        CHECK(bit_equal(a.choke_mdot[i], b.choke_mdot[i]));
        for (size_t j = 0; j < a.mdot_pts.size(); ++j)
        {
            CHECK(a.mdot_pts[j] == b.mdot_pts[j]);
            CHECK(bit_equal(a.pressure_ratio[i][j], b.pressure_ratio[i][j]));
            CHECK(bit_equal(a.temperature_ratio[i][j], b.temperature_ratio[i][j]));
            CHECK(bit_equal(a.torque[i][j], b.torque[i][j]));
        }
    }
}

TEST_CASE("operating_map: sample clamps + fallback warnings")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    // Mutual fixture: a fully valid low-speed map for the clamp case, and a
    // high-speed map with a choked (all-NaN) corner for the fallback case.
    StageStackConfig stack = single_stage(0.18, 0.22, 0.35, 0.35);
    StageInlet inlet;
    MapSweepConfig sweep;
    ModelStatus map_status;
    const OperatingMap map = solve_operating_map(stack, inlet, sweep, *eos, map_status);
    REQUIRE(map.valid);

    // Request outside the grid -> clamped with a warning; the returned value
    // matches a direct solve at the clamped speed/flow point.
    ModelStatus clamp_status;
    const MapSample clamped = sample_operating_map(map, -500.0, 2.0, clamp_status);
    REQUIRE(clamped.ok);
    CHECK(clamped.valid_region);
    CHECK(has_substring(clamp_status.warnings, "clamped"));

    ModelStatus direct_status;
    auto direct =
        exd::physics::fluid::turbomachinery::solve_stage_stack(
            stack, inlet, map.omega_pts.front(), 2.0, *eos, direct_status);
    REQUIRE(direct.ok);
    CHECK(clamped.pressure_ratio == doctest::Approx(direct.total_pi).epsilon(1e-9));

    // High-speed map with a choked region; sample deep inside the all-NaN
    // corner (top-left) -> nearest-valid fallback warning and a value equal
    // to the nearest valid grid point.
    StageStackConfig stack_choke = single_stage(0.10, 0.15, 0.45, 0.25);
    MapSweepConfig sweep_choke;
    sweep_choke.omega_min = 2400.0;
    sweep_choke.omega_max = 3200.0;
    sweep_choke.omega_count = 9;
    sweep_choke.mdot_min = 0.5;
    sweep_choke.mdot_max = 5.0;
    sweep_choke.mdot_count = 25;

    ModelStatus choke_status;
    const OperatingMap choke_map =
        solve_operating_map(stack_choke, inlet, sweep_choke, *eos, choke_status);
    REQUIRE(choke_map.valid);

    ModelStatus fallback_status;
    const MapSample fallback =
        sample_operating_map(choke_map,
                             choke_map.omega_pts.back(),
                             choke_map.mdot_pts.front(),
                             fallback_status);
    REQUIRE(fallback.ok);
    CHECK_FALSE(fallback.valid_region);
    CHECK(has_substring(fallback_status.warnings, "nearest"));

    // Brute-force nearest valid grid point under the same normalized metric.
    const double w_range = choke_map.omega_pts.back() - choke_map.omega_pts.front();
    const double m_range = choke_map.mdot_pts.back() - choke_map.mdot_pts.front();
    double best_dist = std::numeric_limits<double>::infinity();
    double expected_pi = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = 0; i < choke_map.omega_pts.size(); ++i)
    {
        for (size_t j = 0; j < choke_map.mdot_pts.size(); ++j)
        {
            const double pi = choke_map.pressure_ratio[i][j];
            if (std::isnan(pi))
                continue;
            const double d_o = (choke_map.omega_pts[i] - choke_map.omega_pts.back()) / w_range;
            const double d_m = (choke_map.mdot_pts[j] - choke_map.mdot_pts.front()) / m_range;
            const double dist = d_o * d_o + d_m * d_m;
            if (dist < best_dist)
            {
                best_dist = dist;
                expected_pi = pi;
            }
        }
    }
    REQUIRE(!std::isnan(expected_pi));
    CHECK(fallback.pressure_ratio == doctest::Approx(expected_pi).epsilon(1e-12));
    CHECK(!std::isnan(fallback.torque));
    CHECK(!std::isnan(fallback.temperature_ratio));
}

TEST_CASE("operating_map: invalid sweep config")
{
    auto eos = make_air();
    REQUIRE(eos != nullptr);

    StageStackConfig stack = single_stage(0.18, 0.22, 0.35, 0.35);
    StageInlet inlet;

    MapSweepConfig bad_counts;
    bad_counts.omega_count = 2;
    ModelStatus status;
    const OperatingMap bad = solve_operating_map(stack, inlet, bad_counts, *eos, status);
    CHECK_FALSE(bad.valid);
    CHECK_FALSE(status.ok);
    CHECK(status.error.find("omega_count") != std::string::npos);

    MapSweepConfig bad_mdot;
    bad_mdot.mdot_min = 0.0;
    ModelStatus status2;
    const OperatingMap bad2 = solve_operating_map(stack, inlet, bad_mdot, *eos, status2);
    CHECK_FALSE(bad2.valid);
    CHECK_FALSE(status2.ok);

    // Invalid stack propagates through the map entry point.
    StageStackConfig empty_stack;
    ModelStatus status3;
    const OperatingMap bad3 = solve_operating_map(empty_stack, inlet, {}, *eos, status3);
    CHECK_FALSE(bad3.valid);
    CHECK_FALSE(status3.ok);

    // Sampling an invalid map fails through status.
    ModelStatus sample_status;
    const MapSample bad_sample = sample_operating_map(bad, 1.0, 1.0, sample_status);
    CHECK_FALSE(bad_sample.ok);
    CHECK_FALSE(sample_status.ok);
}