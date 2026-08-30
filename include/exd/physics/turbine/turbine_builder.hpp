#pragma once

// ─────────────────────────────────────────────────────
// Parametric turbine builder.
//
// Constructs an exd-geometry TurbineDefinition from a
// handful of engineering parameters so real runs don't
// require hand-writing control-point fixtures. Uniform
// chord, linear twist (hub → tip), constant-radius hub
// and shroud; 90° / edge-on sections allowed by the
// geometry module, resolved here to plain values.
// ─────────────────────────────────────────────────────

#include <exd/geometry/turbine.hpp>
#include <exd/physics/model_status.hpp>

#include <string>

namespace exd::physics::turbine {

struct TurbineBuilderConfig
{
    // Rotor geometry
    double hub_radius = 0.4;       // m, > 0
    double tip_radius = 2.0;       // m, > hub_radius
    double chord = 0.5;            // m, > 0 (uniform)
    double twist_hub_deg = 0.0;    // stagger at the hub
    double twist_tip_deg = -2.0;   // stagger at the tip (linear in span)
    double rpm = 30.0;             // informational design speed (> 0)
    int blade_count = 3;           // ≥ 1
    int section_count = 3;         // ≥ 1 (sections along the span)

    // Axial placement of the rotor plane
    double leading_edge_z = 0.3;   // m (rotor plane ≈ LE + chord/2)
    double duct_length = 2.0;      // m, ≥ LE_z + chord + margin (flow path)

    // Duct / shroud
    double shroud_radius = 0.0;    // 0 → tip_radius (open rotor)

    std::string default_airfoil = "naca0012"; // naming hint for the apps
};

/// Build a TurbineDefinition from the parameters. Sets
/// `status.ok = false` with a message on invalid input.
/// The built definition is a single rotor row; all other
/// fields keep exd-geometry defaults.
exd::geometry::TurbineDefinition make_turbine_definition(
    const TurbineBuilderConfig& config,
    ModelStatus& status);

} // namespace exd::physics::turbine
