#pragma once

// Umbrella header for extropian-physics.
#include <exd/physics/solver/plugin_interface.hpp>

// Level-3 duct-coupled BEM turbine solver (dependency order).
#include <exd/physics/fluid/reduced_order/bem/bem_config.hpp>
#include <exd/physics/fluid/reduced_order/bem/airfoil.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_result.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>
// Additional headers will be added as public types are promoted from src/ to include/
