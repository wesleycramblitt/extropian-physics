// turbine_builder.cpp
// Parametric TurbineDefinition construction (Wave 7: real-run ergonomics).

#include <exd/engine/presets/turbine/turbine_builder.hpp>

#include <cmath>

namespace exd::engine::presets::turbine {

using exd::geometry::BladeRow;
using exd::geometry::BladeRowType;
using exd::geometry::BladeSection;
using exd::geometry::TurbineDefinition;

exd::geometry::TurbineDefinition make_turbine_definition(
    const TurbineBuilderConfig& config, ModelStatus& status)
{
    if (config.hub_radius <= 0.0) { status.ok = false; status.error = "hub_radius must be > 0"; return {}; }
    if (config.tip_radius <= config.hub_radius) { status.ok = false; status.error = "tip_radius must be > hub_radius"; return {}; }
    if (config.chord <= 0.0) { status.ok = false; status.error = "chord must be > 0"; return {}; }
    if (config.rpm <= 0.0) { status.ok = false; status.error = "rpm must be > 0"; return {}; }
    if (config.blade_count < 1) { status.ok = false; status.error = "blade_count must be >= 1"; return {}; }
    if (config.section_count < 1) { status.ok = false; status.error = "section_count must be >= 1"; return {}; }
    if (config.leading_edge_z < 0.0) { status.ok = false; status.error = "leading_edge_z must be >= 0"; return {}; }
    if (config.duct_length < config.leading_edge_z + config.chord + 0.1)
    {
        status.ok = false;
        status.error = "duct_length must cover the blade row (>= LE_z + chord + 0.1)";
        return {};
    }
    if (config.shroud_radius > 0.0 && config.shroud_radius < config.tip_radius)
    {
        status.ok = false;
        status.error = "shroud_radius must be >= tip_radius (or 0 for an open rotor)";
        return {};
    }

    const double shroud = config.shroud_radius > 0.0 ? config.shroud_radius
                                                     : config.tip_radius;
    const double te_z = config.leading_edge_z + config.chord;

    TurbineDefinition t;
    t.flow_path.hub_points = {{0.0f, static_cast<float>(config.hub_radius)},
                              {static_cast<float>(config.duct_length),
                               static_cast<float>(config.hub_radius)}};
    t.flow_path.shroud_points = {{0.0f, static_cast<float>(shroud)},
                                 {static_cast<float>(config.duct_length),
                                  static_cast<float>(shroud)}};

    BladeRow row;
    row.type = BladeRowType::Rotor;
    row.blade_count.value = static_cast<float>(config.blade_count);
    row.rotational_speed.value = static_cast<float>(config.rpm);
    row.leading_edge_hub = {static_cast<float>(config.leading_edge_z),
                            static_cast<float>(config.hub_radius)};
    row.leading_edge_shroud = {static_cast<float>(config.leading_edge_z),
                               static_cast<float>(config.tip_radius)};
    row.trailing_edge_hub = {static_cast<float>(te_z),
                             static_cast<float>(config.hub_radius)};
    row.trailing_edge_shroud = {static_cast<float>(te_z),
                                static_cast<float>(config.tip_radius)};

    const int n = config.section_count;
    for (int i = 0; i < n; ++i)
    {
        BladeSection s;
        const double f = n == 1 ? 0.5 : static_cast<double>(i) / (n - 1);
        s.span = static_cast<float>(f);
        s.chord.value = static_cast<float>(config.chord);
        s.stagger.value =
            static_cast<float>(config.twist_hub_deg
                               + f * (config.twist_tip_deg - config.twist_hub_deg));
        row.sections.push_back(s);
    }
    t.blade_rows = {row};
    return t;
}

} // namespace exd::engine::presets::turbine
