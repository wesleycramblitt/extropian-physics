// compression_system.cpp
// Motor-driven compressor + lumped plenum + throttle (with optional speed
// governing). The compressor row is solved through the same stage-stack path
// as the standalone plenum module so the physics is identical; the motor and
// governor are updated once per step (engine convention).

#include <exd/engine/physics/fluid/turbomachinery/compression_system.hpp>

#include <exd/engine/output/series_writer.hpp>
#include <exd/engine/physics/thermo/eos.hpp>

#include <cmath>
#include <memory>
#include <string>

namespace exd::engine::physics::fluid::turbomachinery {

namespace
{

double clamp_value(double v, double lo, double hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

bool stage_geometry_ok(const StageGeometryConfig& g)
{
    if (!(g.r_hub > 0.0)) return false;
    if (!(g.r_tip > g.r_hub)) return false;
    if (!(g.alpha_1_rad > -3.0 && g.alpha_1_rad < 3.0)) return false;
    if (!(g.beta_2_rad > -3.0 && g.beta_2_rad < 3.0)) return false;
    return true;
}

} // anonymous namespace

bool validate_compression_system_config(const CompressionSystemConfig& config,
                                        std::string& error,
                                        std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.compressor_stages.stages.empty())
    {
        error = "compression system: compressor_stages.stages must be non-empty";
        return false;
    }
    for (const auto& s : config.compressor_stages.stages)
    {
        if (!stage_geometry_ok(s.geometry))
        {
            error = "compression system: stage geometry invalid "
                    "(r_hub > 0, r_tip > r_hub, angles in (-3, 3) rad)";
            return false;
        }
        const double eta = s.loss.polytropic_efficiency;
        if (!(eta > 0.0 && eta <= 1.0))
        {
            error = "compression system: stage polytropic_efficiency must be in (0, 1]";
            return false;
        }
    }

    {
        std::string perr;
        if (!lumped::validate_plenum_config(config.plenum, perr, warnings))
        {
            error = "compression system: " + perr;
            return false;
        }
    }

    if (!(config.motor.R > 0.0)) { error = "motor.R must be > 0"; return false; }
    if (!(config.motor.L >= 0.0)) { error = "motor.L must be >= 0"; return false; }
    if (!(config.motor.kt > 0.0)) { error = "motor.kt must be > 0 (torque-producing motor)"; return false; }
    if (!(config.motor.ke >= 0.0)) { error = "motor.ke must be >= 0"; return false; }
    if (!(config.motor.v_supply >= 0.0)) { error = "motor.v_supply must be >= 0"; return false; }

    if (!(config.shaft_inertia > 0.0)) { error = "shaft_inertia must be > 0"; return false; }

    if (!(config.throttle_gain >= 0.0)) { error = "throttle_gain must be >= 0"; return false; }
    if (!(config.throttle_gain_min >= 0.0)) { error = "throttle_gain_min must be >= 0"; return false; }
    if (!(config.throttle_gain_max >= config.throttle_gain_min))
    {
        error = "throttle_gain_max must be >= throttle_gain_min";
        return false;
    }
    if (config.throttle_gain < config.throttle_gain_min
        || config.throttle_gain > config.throttle_gain_max)
    {
        if (config.governor_enabled)
        {
            warnings.push_back("throttle_gain outside [min, max]; clamped at simulation start");
        }
        else
        {
            error = "throttle_gain must lie in [throttle_gain_min, throttle_gain_max]";
            return false;
        }
    }

    if (config.governor_enabled)
    {
        const double sp = config.governor_setpoint_omega;
        if (!(sp > 0.0))
        {
            error = "governor_setpoint_omega must be > 0";
            return false;
        }
        if (config.governor_pi.kp == 0.0 && config.governor_pi.ki == 0.0)
            warnings.push_back("governor enabled with kp = ki = 0: throttle gain is frozen");
        if (config.throttle_gain_max - config.throttle_gain_min < 1e-15)
            warnings.push_back("governor gain range [min, max] is collapsed; no regulation authority");
        // NOTE: the PI gains are NOT restricted to be positive here. On this
        // plant the throttle gain opens flow, which INCREASES the compressor
        // load (torque ~ mdot) and therefore LOWERS shaft speed, so a
        // stabilizing loop needs the reflected error sense (negative gains).
        // Both senses are allowed; see the governing test for the tuning.
    }

    if (!(config.dt > 0.0)) { error = "dt must be > 0"; return false; }
    if (config.max_steps == 0) { error = "max_steps must be > 0"; return false; }
    if (config.history_interval == 0) { error = "history_interval must be >= 1"; return false; }

    if (!(config.initial_omega >= 0.0)) { error = "initial_omega must be >= 0"; return false; }
    if (!(config.initial_p_plenum >= 0.0))
    { error = "initial_p_plenum must be >= 0 (0 means p_ambient)"; return false; }
    if (!std::isfinite(config.initial_mdot_duct)) { error = "initial_mdot_duct must be finite"; return false; }

    // Sanity: the shaft should not rotate more than ~1 rad per step at the
    // nominal speed (mirrors the engine's > 5 deg/step warning).
    const double w_nom = config.governor_enabled ? config.governor_setpoint_omega : 100.0;
    if (config.dt * w_nom > 1.0)
    {
        warnings.push_back("dt is large: > 1 rad of shaft rotation per step at nominal omega");
    }
    return true;
}

CompressionSystemResult simulate_compression_system(const CompressionSystemConfig& config,
                                                    exd::engine::core::ModelStatus& status)
{
    CompressionSystemResult result;
    std::string error;
    if (!validate_compression_system_config(config, error, result.warnings))
    {
        result.valid = false;
        result.error = error;
        status = {false, error, result.warnings};
        return result;
    }

    std::unique_ptr<exd::engine::physics::thermo::IEos> eos = exd::engine::physics::thermo::make_ideal_gas(exd::engine::physics::thermo::IdealGasConfig{});
    if (!eos)
    {
        result.valid = false;
        result.error = "compression system: failed to create ideal-gas EOS";
        status = {false, result.error, result.warnings};
        return result;
    }

    exd::engine::physics::electromagnetics::DcMotorModel motor(config.motor);
    std::unique_ptr<exd::engine::physics::control::IController> governor;
    if (config.governor_enabled)
    {
        governor = exd::engine::physics::control::make_pi_controller(config.governor_pi);
        if (!governor)
        {
            result.valid = false;
            result.error = "compression system: invalid governor config";
            status = {false, result.error, result.warnings};
            return result;
        }
    }

    std::unique_ptr<exd::engine::output::CsvSeriesWriter> csv;
    if (!config.csv_path.empty())
    {
        csv = std::make_unique<exd::engine::output::CsvSeriesWriter>(
            config.csv_path,
            std::vector<std::string>{"omega_rad_s",
                                     "p_plenum_pa",
                                     "mdot_duct_kg_s",
                                     "pressure_ratio",
                                     "torque_compressor_Nm",
                                     "torque_motor_Nm",
                                     "throttle_gain_kg_s_pa05"},
            true /* flush per row: real-time */, &status);
        if (!status.ok)
        {
            result.valid = false;
            result.error = status.error;
            return result;
        }
    }

    double omega = config.initial_omega;
    double p = config.initial_p_plenum > 0.0 ? config.initial_p_plenum
                                              : config.plenum.p_ambient;
    double mdot = config.initial_mdot_duct;
    double gain = clamp_value(config.throttle_gain, config.throttle_gain_min,
                              config.throttle_gain_max);

    double t = 0.0;
    double total_motor_work = 0.0;
    double total_compressor_work = 0.0;
    double omega_sum = 0.0;
    uint64_t step_count = 0;

    bool failed = false;
    std::string fail_error;
    CompressionStepResult last;

    exd::engine::numerics::IntegratorConfig integration;
    integration.method = exd::engine::numerics::IntegrationMethod::RK4;

    for (uint64_t step = 0; step < config.max_steps; ++step)
    {
        // ---- Governor: exactly one update per step (engine convention) ----
        if (governor)
        {
            ModelStatus gst;
            const double g = governor->update(config.governor_setpoint_omega,
                                              omega, config.dt, gst);
            if (!gst.ok)
            {
                failed = true;
                fail_error = "governor: " + gst.error;
                break;
            }
            gain = clamp_value(g, config.throttle_gain_min, config.throttle_gain_max);
        }

        // ---- Motor: explicit coupling once per step ----
        ModelStatus mst;
        const double torque_motor = motor.step(config.dt, omega, mst);
        if (!mst.ok)
        {
            failed = true;
            fail_error = "motor: " + mst.error;
            break;
        }

        const double omega_before = omega;

        // ---- Full-state derivative ----
        // Rows 1,2 reuse plenum_derivative with the stage stack as the
        // compressor characteristic (returning p0_out - p_ambient).
        ModelStatus substatus;
        bool sub_failed = false;
        exd::engine::numerics::DerivativeFn deriv = [&](std::span<const double> st,
                                         std::span<double> dst, double)
        {
            if (sub_failed)
            {
                dst[0] = dst[1] = dst[2] = 0.0;
                return;
            }
            const double w = st[0];
            const double pw = st[1];
            const double mc = st[2];

            // The stage stack requires mdot > 0. At zero/reverse flow evaluate
            // at a tiny flow floor: this recovers the blocked-flow pressure
            // rise (the impeller still does Euler work at zero through-flow)
            // with negligible torque, so a drive can raise plenum pressure
            // from rest (initial mdot = 0) and develop flow.
            const double mc_solve = mc > 0.0 ? mc : 1.0e-9;
            StageStackResult sres = solve_stage_stack(
                config.compressor_stages,
                StageInlet{config.plenum.p_ambient, config.plenum.T_ambient, 0.0},
                w, mc_solve, *eos, substatus);
            if (!sres.ok || !substatus.ok)
            {
                sub_failed = true;
                if (substatus.error.empty())
                    substatus.error = "compressor stack solve failed in system derivative";
                dst[0] = dst[1] = dst[2] = 0.0;
                return;
            }

            const double torque_comp = sres.total_torque;   // positive resisting
            const double dp_c = sres.p0_out - config.plenum.p_ambient;
            dst[0] = (torque_motor - torque_comp) / config.shaft_inertia;

            lumped::PlenumState ps;
            ps.p_plenum = pw;
            ps.mdot_duct = mc;
            lumped::PlenumCharacteristic comp_char = [&](double) -> double
            {
                return dp_c;
            };
            lumped::PlenumCharacteristic th_char = [&](double pp) -> double
            {
                const double dpp = pp - config.plenum.p_ambient;
                return gain * std::sqrt(dpp > 0.0 ? dpp : 0.0);
            };
            lumped::PlenumDerivative pd = lumped::plenum_derivative(
                config.plenum, ps, comp_char, th_char, *eos, substatus);
            if (!pd.ok || !substatus.ok)
            {
                sub_failed = true;
                if (substatus.error.empty())
                    substatus.error = "plenum derivative failed in system derivative";
                dst[0] = dst[1] = dst[2] = 0.0;
                return;
            }
            dst[1] = pd.dp_dt;
            dst[2] = pd.dmdot_dt;
        };

        double st_arr[3] = {omega, p, mdot};
        double dt_used = config.dt;
        ModelStatus ist;
        const bool iok = exd::engine::numerics::integrate_step(integration, t, config.dt,
                                                std::span<double>(st_arr, 3),
                                                deriv, ist, &dt_used);
        if (!iok || !ist.ok || sub_failed)
        {
            failed = true;
            fail_error = "compression system: integration failed: "
                         + (sub_failed ? (substatus.error.empty() ? "state evaluation failed"
                                                                  : substatus.error)
                                       : ist.error);
            break;
        }
        if (!(dt_used > 0.0)) dt_used = config.dt;

        omega = st_arr[0];
        p = st_arr[1];
        mdot = st_arr[2];

        if (!std::isfinite(omega) || !std::isfinite(p) || !std::isfinite(mdot))
        {
            failed = true;
            fail_error = "compression system: non-finite state after step";
            break;
        }
        if (p <= 0.0)
        {
            failed = true;
            fail_error = "compression system: p_plenum non-positive after step";
            break;
        }

        // ---- Final-state stack evaluation for the recorded step ----
        const double mc_solve = mdot > 0.0 ? mdot : 1.0e-9;
        ModelStatus fst;
        StageStackResult fres = solve_stage_stack(
            config.compressor_stages,
            StageInlet{config.plenum.p_ambient, config.plenum.T_ambient, 0.0},
            omega, mc_solve, *eos, fst);
        if (!fres.ok || !fst.ok)
        {
            failed = true;
            fail_error = "compression system: compressor stack solve failed at "
                         "the step end state: "
                         + (fst.error.empty() ? "mdot out of solve domain (choked?)"
                                              : fst.error);
            break;
        }
        const double torque_compressor = fres.total_torque;
        const double pressure_ratio = fres.p0_out / config.plenum.p_ambient;

        t += dt_used;

        // ---- Energy bookkeeping ----
        // Rectangle rule in time with the step-average shaft speed; using the
        // same omega_avg for motor and compressor makes the stored kinetic
        // energy balance exact for constant torques (RK4 error is next-order).
        const double omega_avg = 0.5 * (omega_before + omega);
        total_motor_work += torque_motor * omega_avg * dt_used;
        total_compressor_work += torque_compressor * omega_avg * dt_used;
        omega_sum += omega;
        ++step_count;

        CompressionStepResult step_result;
        step_result.ok = true;
        step_result.status.ok = true;
        step_result.t = t;
        step_result.dt_used = dt_used;
        step_result.omega = omega;
        step_result.p_plenum = p;
        step_result.mdot_duct = mdot;
        step_result.pressure_ratio = pressure_ratio;
        step_result.torque_compressor = torque_compressor;
        step_result.torque_motor = torque_motor;
        step_result.power_compressor = torque_compressor * omega;
        step_result.power_motor = torque_motor * omega;
        step_result.throttle_gain = gain;
        last = step_result;

        if (config.record_history && (step % config.history_interval == 0))
            result.history.push_back(step_result);

        if (csv)
        {
            csv->write_row(t, std::vector<double>{omega, p, mdot, pressure_ratio,
                                                  torque_compressor, torque_motor,
                                                  gain});
        }
    }

    if (csv) csv->close();

    if (failed)
    {
        result.valid = false;
        result.error = fail_error;
        status = {false, fail_error, result.warnings};
        return result;
    }

    result.valid = true;
    result.error.clear();
    result.final_step = last;
    result.total_time = t;
    result.mean_omega = step_count > 0 ? omega_sum / static_cast<double>(step_count) : 0.0;

    // Settled-window averages over the last history entries.
    const auto& hist = result.history;
    if (!hist.empty())
    {
        const size_t n = hist.size();
        const size_t window = n > 200 ? 200 : n;
        double pr_sum = 0.0;
        double md_sum = 0.0;
        for (size_t i = n - window; i < n; ++i)
        {
            pr_sum += hist[i].pressure_ratio;
            md_sum += hist[i].mdot_duct;
        }
        result.settle_pressure_ratio = pr_sum / static_cast<double>(window);
        result.settle_mdot = md_sum / static_cast<double>(window);
    }
    else
    {
        result.settle_pressure_ratio = last.pressure_ratio;
        result.settle_mdot = last.mdot_duct;
    }

    result.total_motor_work = total_motor_work;
    result.total_compressor_work = total_compressor_work;
    if (total_motor_work > 1e-9)
        result.efficiency_estimate = total_compressor_work / total_motor_work;

    const double dke = 0.5 * config.shaft_inertia
                       * (omega * omega - config.initial_omega * config.initial_omega);
    const double imbalance = total_motor_work - total_compressor_work - dke;
    result.energy_balance_closed = total_motor_work > 1e-9
                                   && std::fabs(imbalance) / total_motor_work < 0.05;

    status.ok = true;
    status.error.clear();
    status.warnings = result.warnings;
    return result;
}

} // namespace exd::engine::physics::fluid::turbomachinery