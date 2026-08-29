#pragma once

#include "fdm_config.hpp"
#include "fdm_result.hpp"

namespace exd::physics::fluid::fdm {

/// Solve 2D incompressible Navier-Stokes using finite differences.
///
/// @param config  Solver configuration (grid, physics, BCs, time integration)
/// @return        Result containing final fields, convergence history, diagnostics
///
/// This is a Tier 1 (CPU single-threaded) solver. The implementation is designed
/// for future multithreading (element-level parallelism) and GPU porting (flat
/// array layouts, trivially copyable result types).
FDMResult solve_fdm(const FDMConfig& config);

/// Advance the solver by a single time step.
/// Useful for external time-stepping loops or coupling scenarios.
///
/// @param config   Solver configuration
/// @param t        Current time (updated in-place)
/// @param step     Current step number (updated in-place)
/// @param field    Current field state (updated in-place)
/// @return         Step diagnostics
FDMStepResult step_fdm(const FDMConfig& config, double& t, int& step, FDMFieldData& field);

} // namespace exd::physics::fluid::fdm
