// plenum.cpp
// Lumped Greitzer-style surge cell with isentropic plenum mapping (see the
// header for the model and the assumption write-up).

#include <exd/engine/physics/fluid/lumped/plenum.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace exd::engine::physics::fluid::lumped {

bool validate_plenum_config(const PlenumModelConfig& config,
                            std::string& error,
                            std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.volume <= 0.0) { error = "plenum: volume must be > 0"; return false; }
    if (config.duct_area <= 0.0) { error = "plenum: duct_area must be > 0"; return false; }
    if (config.duct_length <= 0.0) { error = "plenum: duct_length must be > 0"; return false; }
    if (config.p_ambient <= 0.0) { error = "plenum: p_ambient must be > 0"; return false; }
    if (config.T_ambient <= 0.0) { error = "plenum: T_ambient must be > 0"; return false; }
    return true;
}

PlenumDerivative plenum_derivative(const PlenumModelConfig& config,
                                   const PlenumState& state,
                                   const PlenumCharacteristic& compressor,
                                   const PlenumCharacteristic& throttle,
                                   const exd::engine::physics::thermo::IEos& eos,
                                   exd::engine::core::ModelStatus& status)
{
    PlenumDerivative out;

    std::string verror;
    if (!validate_plenum_config(config, verror, out.status.warnings))
    {
        out.status = exd::engine::core::ModelStatus{false, verror, out.status.warnings};
        status = out.status;
        return out;
    }

    if (!std::isfinite(state.p_plenum) || state.p_plenum <= 0.0)
    {
        const std::string err = "plenum_derivative: p_plenum must be finite and > 0";
        out.status = exd::engine::core::ModelStatus{false, err, out.status.warnings};
        status = out.status;
        return out;
    }
    if (!std::isfinite(state.mdot_duct))
    {
        const std::string err = "plenum_derivative: mdot_duct must be finite";
        out.status = exd::engine::core::ModelStatus{false, err, out.status.warnings};
        status = out.status;
        return out;
    }

    const double gamma = eos.gamma();
    const double R = eos.gas_constant();
    if (!(gamma > 1.0) || !(R > 0.0))
    {
        const std::string err = "plenum_derivative: EOS must provide gamma > 1 and gas_constant > 0";
        out.status = exd::engine::core::ModelStatus{false, err, out.status.warnings};
        status = out.status;
        return out;
    }

    // Isentropic plenum mapping: T_p = T_amb * (p/p_amb)^((gamma-1)/gamma).
    const double inertia = config.duct_length / config.duct_area;
    const double T_p = config.T_ambient
                       * std::pow(state.p_plenum / config.p_ambient,
                                  (gamma - 1.0) / gamma);
    const double a_sq = gamma * R * T_p;

    const double dp_c = compressor(state.mdot_duct);
    const double m_th = throttle(state.p_plenum);
    if (!std::isfinite(dp_c) || !std::isfinite(m_th))
    {
        const std::string err = "plenum_derivative: characteristic returned a non-finite value";
        out.status = exd::engine::core::ModelStatus{false, err, out.status.warnings};
        status = out.status;
        return out;
    }

    out.inertia = inertia;
    out.speed_of_sound_sq = a_sq;
    out.dmdot_dt = (dp_c - (state.p_plenum - config.p_ambient)) / inertia;
    out.dp_dt = (a_sq / config.volume) * (state.mdot_duct - m_th);

    out.ok = true;
    out.status = exd::engine::core::ModelStatus{true, "", out.status.warnings};
    status = out.status;
    return out;
}

PlenumState step_plenum(double dt,
                        const PlenumModelConfig& config,
                        const PlenumState& state,
                        const PlenumCharacteristic& compressor,
                        const PlenumCharacteristic& throttle,
                        const exd::engine::physics::thermo::IEos& eos,
                        const exd::engine::numerics::IntegratorConfig& integration,
                        exd::engine::core::ModelStatus& status)
{
    if (dt <= 0.0)
    {
        status.ok = false;
        status.error = "step_plenum: dt must be > 0";
        return state;
    }

    std::string verror;
    if (!validate_plenum_config(config, verror, status.warnings))
    {
        status.ok = false;
        status.error = verror;
        return state;
    }

    if (!std::isfinite(state.p_plenum) || state.p_plenum <= 0.0)
    {
        status.ok = false;
        status.error = "step_plenum: p_plenum must be finite and > 0";
        return state;
    }
    if (!std::isfinite(state.mdot_duct))
    {
        status.ok = false;
        status.error = "step_plenum: mdot_duct must be finite";
        return state;
    }

    const double gamma = eos.gamma();
    const double R = eos.gas_constant();
    if (!(gamma > 1.0) || !(R > 0.0))
    {
        status.ok = false;
        status.error = "step_plenum: EOS must provide gamma > 1 and gas_constant > 0";
        return state;
    }

    // Integration sanity: warn when dt covers a meaningful fraction of the
    // lumped acoustic period (omega_n = sqrt(a^2/(V*I))).
    const double inertia = config.duct_length / config.duct_area;
    const double T_p = config.T_ambient
                       * std::pow(state.p_plenum / config.p_ambient,
                                  (gamma - 1.0) / gamma);
    const double omega_n = std::sqrt((gamma * R * T_p) / (config.volume * inertia));
    const double phase = dt * omega_n;
    if (phase > 1.0)
    {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "step_plenum: dt = %g covers a significant "
                                        "fraction of the plenum acoustic period "
                                        "(dt*omega_n = %.3g); use a smaller dt for "
                                        "an explicit method", dt, phase);
        status.warnings.push_back(buf);
    }

    // The derivative closure keeps characteristic failures contained in a
    // per-call status so a bad evaluation fails the step cleanly instead of
    // silently poisoning the integrator.
    exd::engine::core::ModelStatus substatus;
    bool sub_failed = false;
    exd::engine::numerics::DerivativeFn deriv = [&](std::span<const double> st,
                                     std::span<double> dst, double)
    {
        if (sub_failed)
        {
            dst[0] = 0.0;
            dst[1] = 0.0;
            return;
        }
        PlenumState ps;
        ps.p_plenum = st[0];
        ps.mdot_duct = st[1];
        PlenumDerivative pd = plenum_derivative(config, ps, compressor, throttle,
                                                eos, substatus);
        if (!pd.ok || !substatus.ok)
        {
            sub_failed = true;
            dst[0] = 0.0;
            dst[1] = 0.0;
            return;
        }
        dst[0] = pd.dp_dt;
        dst[1] = pd.dmdot_dt;
    };

    double st_arr[2] = {state.p_plenum, state.mdot_duct};
    double dt_used = 0.0;
    exd::engine::core::ModelStatus ist;
    const bool ok = exd::engine::numerics::integrate_step(integration, 0.0, dt,
                                           std::span<double>(st_arr, 2), deriv,
                                           ist, &dt_used);
    if (!ok || !ist.ok || sub_failed)
    {
        const std::string err = "step_plenum: integration failed: "
                                + (ist.error.empty() ? "derivative evaluation failed"
                                                     : ist.error);
        status.ok = false;
        status.error = err;
        return state;
    }

    PlenumState out;
    out.p_plenum = st_arr[0];
    out.mdot_duct = st_arr[1];
    status.ok = true;
    status.error.clear();
    return out;
}

} // namespace exd::engine::physics::fluid::lumped