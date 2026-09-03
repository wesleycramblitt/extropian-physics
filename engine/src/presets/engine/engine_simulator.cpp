// engine_simulator.cpp
// Slider-crank engine over the shared integrator module:
//   θ̇ = ω
//   J(θ)·ω̇ = M_gas(θ) − M_load(ω) − ½·(dJ/dθ)·ω²
// The governor (optional PI) is updated ONCE per step — the
// resulting throttle is held constant inside the integrator
// so RK4 stage evaluations never mutate controller state.

#include <exd/engine/presets/engine/engine_simulator.hpp>

#include <exd/engine/output/series_writer.hpp>
#include <exd/engine/physics/thermo/steam.hpp>

#include <cmath>
#include <memory>

#include "engine_internal.hpp"

namespace exd::engine::presets::engine {

namespace
{
constexpr double PI = 3.14159265358979323846;
constexpr double DEG = PI / 180.0;

double clamp(double v, double lo, double hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

double cycle_period(const EngineConfig& config)
{
    return config.thermo.cycle == EngineCycleType::Steam ? 2.0 * PI : 4.0 * PI;
}

double engine_volume_at_cutoff(const EngineConfig& config, double cutoff_deg)
{
    return cylinder_volume(cutoff_deg * (PI / 180.0), config.geometry);
}

double p_back_pressure(const EngineThermoConfig& t)
{
    // Crankcase back pressure for the gas force.
    if (t.cycle == EngineCycleType::Steam) return t.p_condenser;
    return t.p_exhaust;
}
} // anonymous namespace

void engine_rhs(const EngineConfig& config, double throttle,
                std::span<const double> state, std::span<double> dstate)
{
    const double theta = state[0];
    const double w = state[1];
    const CrankKinematics k = crank_kinematics(theta, config.geometry);
    EngineThermoConfig thermo = config.thermo;
    if (config.thermo.cycle != EngineCycleType::Steam) thermo.q_in_cycle *= throttle;
    double p_cyl = 0.0, T_cyl = 0.0;
    cylinder_state(theta, config.geometry, thermo, p_cyl, T_cyl);
    const double f_gas = -(p_cyl - p_back_pressure(config.thermo))
                         * piston_area(config.geometry);
    const double m_gas = f_gas * k.dx_dtheta;
    const double m_load = load_moment(config.load, w);
    dstate[0] = w;
    dstate[1] = (m_gas - m_load - 0.5 * k.dj_dtheta * w * w) / k.j_eq;
}

bool validate_engine_config(const EngineConfig& config,
                            std::string& error,
                            std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();
    const auto& g = config.geometry;
    const auto& t = config.thermo;

    if (g.crank_radius <= 0.0) { error = "crank_radius must be > 0"; return false; }
    if (g.rod_length <= g.crank_radius) { error = "rod_length must be > crank_radius"; return false; }
    if (g.bore <= 0.0) { error = "bore must be > 0"; return false; }
    if (g.clearance_volume <= 0.0) { error = "clearance_volume must be > 0 (keeps pV^γ finite)"; return false; }
    if (g.piston_mass < 0.0) { error = "piston_mass must be >= 0"; return false; }
    if (g.flywheel_inertia <= 0.0) { error = "flywheel_inertia must be > 0"; return false; }

    if (t.r_gas <= 0.0) { error = "r_gas must be > 0"; return false; }
    if (t.p_intake <= 0.0) { error = "p_intake must be > 0"; return false; }
    if (t.T_intake <= 0.0) { error = "T_intake must be > 0"; return false; }
    if (t.q_in_cycle < 0.0) { error = "q_in_cycle must be >= 0"; return false; }
    if (t.wiebe_burn_duration_deg <= 0.0) { error = "wiebe_burn_duration_deg must be > 0"; return false; }

    if (t.cycle == EngineCycleType::Otto)
    {
        if (t.gamma_compression <= 1.0 || t.gamma_expansion <= 1.0)
        { error = "gamma_compression/gamma_expansion must be > 1"; return false; }
        if (std::fabs(t.gamma_expansion - t.gamma_compression) > 1e-9)
        {
            warnings.push_back(
                "gamma_expansion != gamma_compression: the polytropic cycle does "
                "net work per revolution with zero heat release (heat-transfer "
                "stand-in). Keep them equal unless you model that transfer.");
        }
        if (t.p_exhaust <= 0.0) { error = "p_exhaust must be > 0"; return false; }
        if (t.wiebe_ignition_deg < 0.0 || t.wiebe_ignition_deg >= 180.0)
        { error = "wiebe_ignition_deg must be in [0, 180)"; return false; }
        if (config.governor.enabled && t.q_in_cycle <= 0.0)
        { warnings.push_back("governor enabled with q_in_cycle = 0 (no heat to throttle)"); }
    }
    else
    {
        if (t.p_boiler <= t.p_condenser) { error = "p_boiler must be > p_condenser"; return false; }
        if (t.p_condenser <= 0.0) { error = "p_condenser must be > 0"; return false; }
        if (t.steam_cutoff_deg <= 0.0 || t.steam_cutoff_deg >= 180.0)
        { error = "steam_cutoff_deg must be in (0, 180)"; return false; }
        if (!(t.steam_quality_cutoff > 0.0 && t.steam_quality_cutoff <= 1.0))
        { error = "steam_quality_cutoff must be in (0, 1]"; return false; }
        if (t.steam_gamma <= 1.0) { error = "steam_gamma must be > 1"; return false; }
        if (config.governor.enabled)
        { warnings.push_back("governor ignored for steam cycles (fixed admission, placeholder model)"); }
    }

    if (config.load.generator_enabled)
    {
        const auto& wp = config.load.generator_omega_pts;
        const auto& tp = config.load.generator_torque_pts;
        if (wp.size() != tp.size() || wp.size() < 2)
        { error = "generator curve needs >= 2 matching omega/torque points"; return false; }
        for (size_t i = 1; i < wp.size(); ++i)
            if (!(wp[i] > wp[i - 1]))
            { error = "generator_omega_pts must be strictly increasing"; return false; }
    }

    if (config.dt <= 0.0) { error = "dt must be > 0"; return false; }
    if (config.max_steps == 0) { error = "max_steps must be > 0"; return false; }
    if (config.history_interval == 0) { error = "history_interval must be >= 1"; return false; }

    // Sanity: crank should not advance more than ~5° per step at the
    // governor setpoint (or a nominal 100 rad/s when ungoverned).
    const double w_nom = config.governor.enabled ? config.governor.setpoint_omega : 100.0;
    if (config.dt * w_nom > 5.0 * DEG)
    {
        warnings.push_back("dt is large: crank advances > 5° per step at nominal omega");
    }
    return true;
}

EngineStepResult step_engine(EngineState& state,
                             double t,
                             const EngineConfig& config,
                             exd::engine::physics::control::IController* governor,
                             ModelStatus& status)
{
    EngineStepResult result;
    result.status = status;
    result.t = t;

    std::string error;
    if (!validate_engine_config(config, error, result.status.warnings))
    {
        status.ok = false;
        status.error = error;
        result.status = status;
        result.ok = false;
        return result;
    }

    // ── Governor: exactly one update per step ─────────
    double throttle = 1.0;
    if (governor)
    {
        ModelStatus cst;
        throttle = governor->update(config.governor.setpoint_omega, state.omega,
                                    config.dt, cst);
        if (!cst.ok)
        {
            status.ok = false;
            status.error = "governor: " + cst.error;
            result.status = status;
            return result;
        }
        throttle = clamp(throttle, config.governor.throttle_min, config.governor.throttle_max);
    }

    // Heat release seen by the RHS is scaled by the throttle (Otto;
    // steam ignores throttle — admission is config-fixed).

    exd::engine::numerics::DerivativeFn rhs = [&](std::span<const double> st, std::span<double> dst, double)
    {
        engine_rhs(config, throttle, st, dst);
    };

    double state_arr[2] = {state.theta_rad, state.omega};
    double dt_used = config.dt;
    ModelStatus ist;
    const bool ok = exd::engine::numerics::integrate_step(config.integration, t, config.dt,
                                           std::span<double>(state_arr, 2), rhs,
                                           ist, &dt_used);
    if (!ok || !ist.ok)
    {
        status.ok = false;
        status.error = "engine: integration failed: " + ist.error;
        result.status = status;
        result.ok = false;
        return result;
    }

    state.theta_rad = state_arr[0];
    state.omega = state_arr[1];
    state.cycles = static_cast<uint64_t>(
        std::max(0.0, std::floor(state.theta_rad / cycle_period(config))));

    // ── Result at the new state ────────────────────────
    // dt_used is only reported by AdaptiveRK45; fixed-step
    // methods leave it 0 → fall back to config.dt.
    if (!(dt_used > 0.0)) dt_used = config.dt;
    result.dt_used = dt_used;
    result.t = t + dt_used;
    result.ok = true;
    result.state = state;
    result.throttle = throttle;

    const CrankKinematics k = crank_kinematics(state.theta_rad, config.geometry);
    double p_cyl = 0.0, T_cyl = 0.0;
    {
        EngineThermoConfig rthermo = config.thermo;
        if (config.thermo.cycle != EngineCycleType::Steam) rthermo.q_in_cycle *= throttle;
        cylinder_state(state.theta_rad, config.geometry, rthermo, p_cyl, T_cyl);
    }
    result.piston_x = k.x;
    result.piston_v = state.omega * k.dx_dtheta;
    result.p_cyl = p_cyl;
    result.T_cyl = T_cyl;
    result.V_cyl = cylinder_volume(state.theta_rad, config.geometry);
    result.gas_force = -(p_cyl - p_back_pressure(config.thermo)) * piston_area(config.geometry);
    result.indicated_moment = result.gas_force * k.dx_dtheta;
    result.load_moment = load_moment(config.load, state.omega);
    result.power = (result.indicated_moment - result.load_moment) * state.omega;
    return result;
}

EngineSimResult simulate_engine(const EngineConfig& config, ModelStatus& status)
{
    EngineSimResult result;
    std::string error;
    if (!validate_engine_config(config, error, result.warnings))
    {
        result.valid = false;
        result.error = error;
        status = {false, error, result.warnings};
        return result;
    }

    std::unique_ptr<exd::engine::physics::control::IController> governor;
    if (config.governor.enabled)
    {
        governor = exd::engine::physics::control::make_pi_controller(config.governor.pi);
        if (!governor)
        {
            result.valid = false;
            result.error = "engine: invalid governor config";
            status = {false, result.error, result.warnings};
            return result;
        }
    }

    bool csv_ok = config.csv_path.empty();
    std::unique_ptr<exd::engine::output::CsvSeriesWriter> csv;
    if (!config.csv_path.empty())
    {
        csv = std::make_unique<exd::engine::output::CsvSeriesWriter>(
            config.csv_path,
            std::vector<std::string>{"theta_rad", "omega_rad_s", "piston_x_m",
                                     "piston_v_m_s", "p_cyl_pa", "T_cyl_K",
                                     "indicated_moment_Nm", "load_moment_Nm",
                                     "power_W", "throttle", "cycles"},
            true /* flush per row: real-time visible */, &status);
        csv_ok = status.ok;
        if (!csv_ok)
        {
            result.valid = false;
            result.error = status.error;
            return result;
        }
    }

    // Steam boiler heat per cycle: raise the trapped mass from liquid at
    // the condenser saturation temperature to saturated vapor at the boiler.
    // m = ρ_g(p_b)·V_cut·x_cut. (Rankine-lite; documented engineering model.)
    const double steam_heat_per_cycle = [&]
    {
        if (config.thermo.cycle != EngineCycleType::Steam) return 0.0;
        const double t_sat_b = exd::engine::physics::thermo::saturation_temperature(config.thermo.p_boiler);
        const double t_sat_c = exd::engine::physics::thermo::saturation_temperature(config.thermo.p_condenser);
        const double m = exd::engine::physics::thermo::rho_g(config.thermo.p_boiler, t_sat_b)
                         * engine_volume_at_cutoff(config, config.thermo.steam_cutoff_deg)
                         * config.thermo.steam_quality_cutoff;
        const double h_g = exd::engine::physics::thermo::h_g(t_sat_b);
        const double h_f_c = exd::engine::physics::thermo::h_f(t_sat_c);
        return m * (h_g - h_f_c);
    }();

    EngineState state{config.initial_theta_rad, config.initial_omega, 0};
    double t = 0.0;
    double total_work = 0.0;
    double omega_sum = 0.0;
    uint64_t omega_count = 0;
    double throttle_sum = 0.0;
    double heat_sum = 0.0;             // exact released heat: Σ q_in·u·Δθ/4π
    uint64_t completed = 0;
    bool failed = false;
    EngineStepResult last;

    for (uint64_t step = 0; step < config.max_steps && !failed; ++step)
    {
        const double w_before = state.omega;
        const double theta_before = state.theta_rad;
        EngineStepResult r = step_engine(state, t, config, governor.get(), status);
        if (!r.ok || !status.ok)
        {
            failed = true;
            last = r;
            break;
        }
        total_work += r.indicated_moment * 0.5 * (w_before + state.omega) * r.dt_used;
        completed = state.cycles;
        t = r.t;
        last = r;
        omega_sum += state.omega;
        ++omega_count;
        throttle_sum += r.throttle;
        {
            const double dtheta = std::max(0.0, state.theta_rad - theta_before);
            if (config.thermo.cycle == EngineCycleType::Otto)
            {
                heat_sum += config.thermo.q_in_cycle * r.throttle * dtheta / (4.0 * PI);
            }
            else
            {
                heat_sum += steam_heat_per_cycle * dtheta / (2.0 * PI);
            }
        }

        if (config.record_history && (step % config.history_interval == 0))
            result.history.push_back(r);

        if (csv)
            csv->write_row(r.t, std::vector<double>{r.state.theta_rad, r.state.omega,
                                                    r.piston_x, r.piston_v, r.p_cyl,
                                                    r.T_cyl, r.indicated_moment,
                                                    r.load_moment, r.power, r.throttle,
                                                    static_cast<double>(r.state.cycles)});
    }

    if (csv) csv->close();

    if (failed)
    {
        result.valid = false;
        result.error = status.error.empty() ? "engine: step failed" : status.error;
        return result;
    }

    result.valid = true;
    result.final_step = last;
    result.total_time = t;
    result.total_indicated_work = total_work;
    result.mean_omega = omega_count > 0 ? omega_sum / static_cast<double>(omega_count) : 0.0;
    result.cycles_completed = static_cast<double>(completed);
    result.mean_indicated_power = t > 0.0 ? total_work / t : 0.0;
    result.mean_throttle = omega_count > 0 ? throttle_sum / omega_count : 1.0;
    if (heat_sum > 0.0)
        result.efficiency_estimate = total_work / heat_sum;

    return result;
}

} // namespace exd::engine::presets::engine
