#pragma once

#include <exd/geometry/turbine.hpp>

#include "airfoil.hpp"
#include "bem_config.hpp"
#include "bem_result.hpp"

namespace exd::physics::fluid::reduced_order::bem
{

/// Solve a single-rotor duct/hull-coupled BEM turbine model.
///
/// The duct acceleration, hull drag split, Gaussian wake, and Bernoulli
/// pressure field are analytical / engineering estimates, not a CFD solution.
///
/// @param turbine    Parametric turbine definition from exd-geometry.
/// @param conditions Freestream operating conditions.
/// @param polars     Airfoil polar database (must include requested airfoils).
/// @param config     Solver configuration.
/// @return           TurbineResult with rotor/system loads, radial stations,
///                   and an engineering-approximate flow field.
TurbineResult solve_turbine(const exd::geometry::TurbineDefinition& turbine,
                            const OperatingConditions& conditions,
                            const PolarDatabase& polars,
                            const BEMSolverConfig& config = {});

} // namespace exd::physics::fluid::reduced_order::bem
