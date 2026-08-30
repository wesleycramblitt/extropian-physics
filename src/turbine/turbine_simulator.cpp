// turbine_simulator.cpp
// Turbine application assembly: maps an exd-geometry TurbineDefinition onto the
// generic fluid::forces + mechanics modules and wires them into a coupled
// field -> blade surface -> forces -> moments -> generator load -> rotational
// dynamics loop. This is the only turbine-specific layer of the physics stack;
// the same generic modules drive propellers, pumps and fans without it.
//
// Geometry convention follows exd-geometry: rotation axis = +Z, meridional
// points are (z, r), and the rotor plane sits at the mid LE/TE axial station.

#include <exd/physics/turbine/turbine_simulator.hpp>

#include <exd/geometry/turbine.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/fluid/forces/force_evaluator.hpp>
#include <exd/physics/fluid/forces/flow_types.hpp>
#include <exd/physics/mechanics/dynamics.hpp>
#include <exd/physics/mechanics/moment_model.hpp>
#include <exd/physics/mechanics/rotational_state.hpp>
#include <exd/physics/mechanics/status.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace exd::physics::turbine
{

namespace
{

constexpr double kPi = 3.14159265358979323846;

// Purely meridional LE/TE helpers; Vec2f is the exd-geometry (z, r) point.
struct Vec2 { double x = 0.0; double y = 0.0; };

Vec2 to_vec(const math::Vec2f& v)
{
    return {static_cast<double>(v.x), static_cast<double>(v.y)};
}

double length(const Vec2& a, const Vec2& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

Vec2 lerp(const Vec2& a, const Vec2& b, double f)
{
    return {a.x + (b.x - a.x) * f, a.y + (b.y - a.y) * f};
}

/// Interpolate one per-section scalar by span, clamped flat beyond the section
/// table ends (mirrors stagger_at_span in the reduced-order BEM blade
/// geometry). Extrapolation beyond the ends emits a warning.
template <typename Getter>
double section_value_at_span(const exd::geometry::BladeRow& row, double span,
                             Getter getter, std::vector<std::string>& warnings,
                             const char* what)
{
    std::vector<std::pair<double, double>> pts;
    pts.reserve(row.sections.size());
    for (const auto& s : row.sections)
        pts.emplace_back(static_cast<double>(s.span), static_cast<double>(getter(s)));
    std::sort(pts.begin(), pts.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b)
              { return a.first < b.first; });
    if (pts.empty())
        return 0.0; // defensive; callers guard on empty sections themselves

    if (span <= pts.front().first)
    {
        if (span < pts.front().first)
            warnings.push_back(std::string(what) + " extrapolated at low span");
        return pts.front().second;
    }
    if (span >= pts.back().first)
    {
        if (span > pts.back().first)
            warnings.push_back(std::string(what) + " extrapolated at high span");
        return pts.back().second;
    }

    auto it = std::lower_bound(pts.begin(), pts.end(), span,
        [](const std::pair<double, double>& a, double v) { return a.first < v; });
    const std::size_t i = static_cast<std::size_t>(it - pts.begin());
    const std::size_t im1 = i - 1;
    const double f = (span - pts[im1].first) / (pts[i].first - pts[im1].first);
    return pts[im1].second + f * (pts[i].second - pts[im1].second);
}

} // anonymous namespace

fluid::forces::BladeGeometry make_blade_geometry(const exd::geometry::TurbineDefinition& turbine,
                                                 int element_count,
                                                 const std::string& default_airfoil,
                                                 std::vector<std::string>& warnings)
{
    fluid::forces::BladeGeometry geo;
    if (element_count < 4)
    {
        warnings.push_back("turbine: element_count < 4; blade geometry disabled");
        return geo; // default: empty stations, r_hub = r_tip = 0
    }

    // Locate the single Rotor row; every other row type is ignored.
    const exd::geometry::BladeRow* rotor = nullptr;
    int rotor_count = 0;
    for (const auto& row : turbine.blade_rows)
    {
        if (row.type == exd::geometry::BladeRowType::Rotor)
        {
            ++rotor_count;
            rotor = &row;
        }
        else
        {
            warnings.push_back("non-rotor blade row ignored");
        }
    }
    if (rotor_count == 0)
    {
        warnings.push_back("turbine: no Rotor row found");
        return geo;
    }
    if (rotor_count > 1)
    {
        warnings.push_back("turbine: " + std::to_string(rotor_count) +
                           " Rotor rows; single-rotor solver");
        return geo;
    }

    // LE/TE hubs and shrouds are (z, r): radius is .y, axial station is .x.
    const Vec2 le_hub = to_vec(rotor->leading_edge_hub);
    const Vec2 te_hub = to_vec(rotor->trailing_edge_hub);
    const Vec2 le_shroud = to_vec(rotor->leading_edge_shroud);
    const Vec2 te_shroud = to_vec(rotor->trailing_edge_shroud);

    const double r_hub = 0.5 * (le_hub.y + te_hub.y);
    const double r_tip = 0.5 * (le_shroud.y + te_shroud.y);
    const double z_rotor = 0.5 * (0.5 * (le_hub.x + le_shroud.x) +
                                  0.5 * (te_hub.x + te_shroud.x));

    if (!(r_tip > r_hub) || r_hub <= 0.0)
    {
        warnings.push_back("turbine: R_tip <= R_hub or non-positive hub radius");
        return geo;
    }

    const int blade_count = static_cast<int>(rotor->blade_count.value);
    if (blade_count < 1)
    {
        warnings.push_back("turbine: blade_count < 1");
        return geo;
    }

    geo.r_hub = r_hub;
    geo.r_tip = r_tip;
    geo.z_rotor = z_rotor;
    geo.blade_count = blade_count;
    geo.stations.reserve(static_cast<std::size_t>(element_count));

    const bool have_sections = !rotor->sections.empty();
    if (!have_sections)
        warnings.push_back("BladeRow.sections empty; zero-stagger defaults");

    const double dr = (r_tip - r_hub) / static_cast<double>(element_count);
    for (int i = 0; i < element_count; ++i)
    {
        const double r = r_hub + (static_cast<double>(i) + 0.5) * dr;
        const double span = (r - r_hub) / (r_tip - r_hub);

        // Chord follows the reference computation: |te(span) - le(span)|.
        const Vec2 le = lerp(le_hub, le_shroud, span);
        const Vec2 te = lerp(te_hub, te_shroud, span);
        const double chord = length(le, te);

        fluid::forces::BladeStation st;
        st.r = r;
        st.dr = dr;
        st.chord = chord;
        st.twist_deg = have_sections
            ? section_value_at_span(*rotor, span,
                  [](const exd::geometry::BladeSection& s)
                  { return static_cast<double>(s.stagger.value); },
                  warnings, "stagger")
            : 0.0;
        st.thickness_ratio = have_sections
            ? section_value_at_span(*rotor, span,
                  [](const exd::geometry::BladeSection& s)
                  { return static_cast<double>(s.max_thickness.value); },
                  warnings, "thickness")
            : 0.12;
        st.airfoil = default_airfoil;
        geo.stations.push_back(std::move(st));
    }
    return geo;
}

mechanics::CurveMomentConfig make_generator_curve(double rated_power,
                                                  double efficiency,
                                                  double min_omega,
                                                  int points)
{
    mechanics::CurveMomentConfig curve;
    if (rated_power <= 0.0 || efficiency <= 0.0 || min_omega <= 0.0 || points < 2)
        return curve; // empty curve = no load

    const double w_ref = std::max(min_omega, rated_power / (efficiency * min_omega));
    const double w_hi = 4.0 * w_ref;

    curve.omega_pts.reserve(static_cast<std::size_t>(points));
    curve.torque_pts.reserve(static_cast<std::size_t>(points));
    for (int i = 0; i < points; ++i)
    {
        const double f = static_cast<double>(i) / static_cast<double>(points - 1);
        const double w = min_omega + (w_hi - min_omega) * f;
        curve.omega_pts.push_back(w);
        // Opposing moment, positive magnitude per IMomentModel convention.
        curve.torque_pts.push_back(rated_power / (efficiency * w));
    }
    return curve;
}

TurbineStepResult step_turbine(mechanics::RotationalState& state,
                               const exd::geometry::TurbineDefinition& turbine,
                               const coupling::IFlowField3D& flow,
                               const TurbineConfig& config,
                               control::IController* governor)
{
    TurbineStepResult result;
    auto& status = result.status;

    if (config.element_count < 4)
    {
        status.ok = false;
        status.error = "invalid turbine config: element_count < 4";
        return result;
    }
    if (config.inertia <= 0.0)
    {
        status.ok = false;
        status.error = "invalid turbine config: inertia must be positive";
        return result;
    }
    if (config.dt <= 0.0)
    {
        status.ok = false;
        status.error = "invalid turbine config: dt must be positive";
        return result;
    }

    std::vector<std::string> warnings;
    const fluid::forces::BladeGeometry blade =
        make_blade_geometry(turbine, config.element_count, config.default_airfoil, warnings);
    status.warnings.insert(status.warnings.end(), warnings.begin(), warnings.end());

    const bool blade_invalid = blade.stations.empty() ||
                               !(blade.r_tip > blade.r_hub) || blade.r_hub <= 0.0;
    if (blade_invalid)
    {
        status.ok = false;
        status.error = warnings.empty()
            ? "invalid blade geometry"
            : "invalid blade geometry: " + warnings.back();
        return result;
    }

    // exd-geometry machine convention: rotation axis = +Z at the origin.
    const mechanics::RotationAxis axis;

    const fluid::forces::BladeSurface surface =
        fluid::forces::build_blade_surface(blade, state.angle_rad, axis);

    // Pressure reference is irrelevant here: pressure integration over the
    // paired upper/lower surfaces cancels any constant offset, and the
    // coefficient-based (momentum/table) models never read pressure.
    const fluid::forces::SurfaceFlow sampled =
        coupling::sample_flow(flow, surface, /*p_ref*/ 0.0);

    std::unique_ptr<fluid::forces::IForceEvaluator> evaluator =
        fluid::forces::make_force_evaluator(config.force, status);
    if (!evaluator)
    {
        status.ok = false;
        status.error = "failed to create force evaluator";
        return result;
    }

    std::vector<mechanics::ElementForce3D> per_element;
    evaluator->compute(blade, sampled, state.omega, axis, per_element, status);
    if (!status.ok)
        return result; // state unchanged

    const mechanics::MomentResult moments = mechanics::integrate_moment(per_element, axis);
    if (!moments.valid)
    {
        status.ok = false;
        status.error = "axis degenerate";
        return result;
    }

    // Generator load opposes rotation; an empty curve is a documented no-load.
    std::unique_ptr<mechanics::IMomentModel> load =
        mechanics::make_curve_moment(config.generator);
    double external = load->moment(state, status);
    if (!status.ok)
        return result; // state unchanged

    // Speed governor: exactly one update per step. Measurement is
    // (setpoint − ω): positive when the rotor runs FAST → effort up →
    // generator load fraction up → decelerates; below setpoint the
    // fraction drops so the rotor accelerates. Clamped to [min, max].
    double throttle = 1.0;
    if (governor)
    {
        ModelStatus gst;
        const double measurement = config.governor.setpoint_omega - state.omega;
        throttle = governor->update(0.0, measurement, config.dt, gst);
        if (!gst.ok)
        {
            status.ok = false;
            status.error = "governor: " + gst.error;
            return result; // state unchanged
        }
        throttle = std::clamp(throttle, config.governor.throttle_min,
                              config.governor.throttle_max);
        external *= throttle;
    }

    const double omega_before = state.omega;
    const double net = moments.torque - external;

    std::unique_ptr<mechanics::IRotationalDynamics> dynamics =
        mechanics::make_rigid_rotor_dynamics({config.inertia, config.integration});
    const mechanics::RotationalState new_state =
        dynamics->advance(config.dt, net, state, status);
    if (!status.ok)
        return result; // state unchanged

    state = new_state; // caller owns the state

    result.ok = status.ok;
    result.state = state;
    result.aero = moments;
    result.external_moment = external;
    result.power = moments.torque * omega_before; // aero power (W)
    result.throttle = throttle;
    result.per_element = std::move(per_element);
    return result;
}

TurbineSimResult simulate_turbine(const exd::geometry::TurbineDefinition& turbine,
                                  const fluid::forces::Freestream& freestream,
                                  const TurbineConfig& config)
{
    TurbineSimResult result;

    const double v = std::sqrt(freestream.velocity[0] * freestream.velocity[0] +
                               freestream.velocity[1] * freestream.velocity[1] +
                               freestream.velocity[2] * freestream.velocity[2]);
    if (v <= 0.0)
    {
        result.error = "freestream velocity must be non-zero";
        return result;
    }
    if (freestream.rho <= 0.0)
    {
        result.error = "freestream density must be positive";
        return result;
    }
    if (freestream.mu <= 0.0)
    {
        result.error = "freestream viscosity must be positive";
        return result;
    }
    if (config.dt <= 0.0)
    {
        result.error = "config.dt must be positive";
        return result;
    }
    if (config.max_steps <= 0)
    {
        result.error = "config.max_steps must be positive";
        return result;
    }
    if (config.history_interval < 1)
    {
        result.error = "config.history_interval must be >= 1";
        return result;
    }
    if (config.inertia <= 0.0)
    {
        result.error = "config.inertia must be positive";
        return result;
    }
    if (config.element_count < 4)
    {
        result.error = "config.element_count must be >= 4";
        return result;
    }
    if (config.governor.enabled && config.generator.omega_pts.empty())
    {
        result.warnings.push_back(
            "governor enabled with an empty generator curve: nothing to scale");
    }

    coupling::UniformFieldConfig field_config;
    field_config.velocity = freestream.velocity;
    field_config.rho = freestream.rho;
    field_config.mu = freestream.mu;
    field_config.p_ref = freestream.p_ref;
    std::unique_ptr<coupling::IFlowField3D> field =
        coupling::make_uniform_field(field_config);

    // Governor (optional): PI state lives here, one update per step.
    std::unique_ptr<control::IController> governor;
    if (config.governor.enabled)
    {
        governor = control::make_pi_controller(config.governor.pi);
        if (!governor)
        {
            result.error = "invalid governor config (PI)";
            return result;
        }
    }

    mechanics::RotationalState state{config.initial_omega, 0.0};

    TurbineStepResult final_step{};
    double total_energy = 0.0;
    bool failed = false;

    for (int step = 0; step < config.max_steps; ++step)
    {
        TurbineStepResult step_result =
            step_turbine(state, turbine, *field, config, governor.get());
        result.warnings.insert(result.warnings.end(),
                               step_result.status.warnings.begin(),
                               step_result.status.warnings.end());
        total_energy += step_result.power * config.dt;
        if (config.record_history && (step % config.history_interval == 0))
            result.history.push_back(step_result);
        final_step = std::move(step_result);
        if (!final_step.ok)
        {
            result.error = final_step.status.error;
            failed = true;
            break;
        }
    }

    result.final_step = final_step;
    result.total_energy = total_energy;

    // Rotor radius for the dimensionless coefficients; recomputing the blade is
    // cheap and keeps make_blade_geometry the single source of r_tip.
    std::vector<std::string> geo_warnings;
    const fluid::forces::BladeGeometry blade =
        make_blade_geometry(turbine, config.element_count, config.default_airfoil,
                            geo_warnings);
    result.warnings.insert(result.warnings.end(), geo_warnings.begin(), geo_warnings.end());

    const double R = blade.r_tip;
    const double A = kPi * R * R;
    const double power_denom = 0.5 * freestream.rho * A * v * v * v;
    const double force_denom = 0.5 * freestream.rho * A * v * v;
    // Thrust acts downwind (-e_z); axial_force is negative, so -axial_force is
    // the thrust magnitude. Clamped at zero against numerical noise.
    const double thrust_mag = std::max(0.0, -(final_step.aero.axial_force));

    result.cp = (power_denom > 0.0) ? final_step.power / power_denom : 0.0;
    result.ct = (force_denom > 0.0) ? thrust_mag / force_denom : 0.0;
    result.tsr = (v > 0.0) ? final_step.state.omega * R / v : 0.0;

    result.valid = final_step.ok && !failed;
    return result;
}

} // namespace exd::physics::turbine