#pragma once

// Operating map for a stage stack: rectangular sweep over shaft speed and
// mass flow, per-point solve, plus derived surge/choke lines. The vectors
// that hold the grid are deliberately plain row-major config-size data
// (this is not a hot path).

#include <exd/engine/physics/fluid/turbomachinery/stage_stack.hpp>

#include <string>
#include <vector>

namespace exd::engine::physics::fluid::turbomachinery
{

struct MapSweepConfig
{
    double omega_min = 50.0;   // rad/s, >= 0
    double omega_max = 1000.0; // rad/s
    int omega_count = 21;      // >= 3
    double mdot_min = 0.1;     // kg/s, > 0
    double mdot_max = 5.0;     // kg/s
    int mdot_count = 25;       // >= 3
};

bool validate_map_sweep_config(const MapSweepConfig& config,
                               std::string& error,
                               std::vector<std::string>& warnings);

struct OperatingMap
{
    bool valid = false;

    std::vector<double> omega_pts; // speed line centers/points, rad/s
    std::vector<double> mdot_pts;  // mass-flow points, kg/s

    // Grids indexed [omega][mdot]; invalid (failed or choked) points are NaN.
    std::vector<std::vector<double>> pressure_ratio;
    std::vector<std::vector<double>> temperature_ratio;
    std::vector<std::vector<double>> torque;

    // Per speed line: mdot where total-pressure ratio is maximal over the
    // valid points (d pi/d mdot changes sign there). NaN when the line has
    // no valid points.
    std::vector<double> surge_mdot;

    // Per speed line: smallest mdot whose point is choked (the choked region
    // is entered from low mass flow). Choked "clearing" at higher flow is a
    // known mean-line artifact; only the first choke point is recorded.
    // NaN when the line never chokes.
    std::vector<double> choke_mdot;
};

struct MapSample
{
    bool ok = false;
    bool valid_region = true; // false when a nearest-valid fallback was used
    double pressure_ratio = 0.0;
    double temperature_ratio = 0.0;
    double torque = 0.0;
};

/// Rectangular sweep over (omega, mdot). Per point solve_stage_stack; store
/// NaN when the solve fails or any stage chokes. Surge and choke lines are
/// derived per speed line. All stage warnings are accumulated.
OperatingMap solve_operating_map(const StageStackConfig& stack_config,
                                 const StageInlet& inlet,
                                 const MapSweepConfig& sweep,
                                 const exd::engine::physics::thermo::IEos& eos,
                                 exd::engine::core::ModelStatus& status);

/// Bilinear sample of the map in (omega, mdot). Requested points outside the
/// grid extents are clamped (warning). If any of the four cell corners is
/// invalid, fall back to the nearest valid grid point (warning) and set
/// valid_region = false.
MapSample sample_operating_map(const OperatingMap& map,
                               double omega,
                               double mdot,
                               exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::turbomachinery