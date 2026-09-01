// stage.cpp
// Single-stage mean-line closure. See stage.hpp for the derivation summary.

#include <exd/physics/fluid/turbomachinery/stage.hpp>
#include <exd/physics/thermo/polytropic.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace exd::physics::fluid::turbomachinery
{

bool validate_stage_config(const StageConfig& config, std::string& error,
                           std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.geometry.r_hub <= 0.0)
    {
        error = "stage: r_hub must be > 0";
        return false;
    }
    if (config.geometry.r_tip <= 0.0)
    {
        error = "stage: r_tip must be > 0";
        return false;
    }
    if (config.geometry.r_hub >= config.geometry.r_tip)
    {
        error = "stage: r_hub must be < r_tip";
        return false;
    }
    if (!(config.loss.polytropic_efficiency > 0.0 &&
          config.loss.polytropic_efficiency <= 1.0))
    {
        error = "stage: polytropic_efficiency must be in (0, 1]";
        return false;
    }

    // Mean-line validity envelope: strong free-vortex twist is not captured.
    if (config.geometry.hub_tip_ratio() < 0.5)
    {
        warnings.push_back(
            "stage: hub/tip ratio < 0.5 exceeds the mean-line validity envelope");
    }
    return true;
}

namespace
{
constexpr int kMaxFixedPointIterations = 50;
constexpr double kFixedPointTolerance = 1e-10;
} // anonymous namespace

StageResult solve_stage(const StageConfig& config, const StageInlet& inlet,
                        double omega, double mdot, const thermo::IEos& eos,
                        exd::physics::ModelStatus& status)
{
    StageResult result;
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    std::string verror;
    std::vector<std::string> vwarnings;
    if (!validate_stage_config(config, verror, vwarnings))
    {
        status.ok = false;
        status.error = verror;
        result.status = status;
        return result;
    }
    status.warnings.insert(status.warnings.end(), vwarnings.begin(), vwarnings.end());

    if (!(mdot > 0.0))
    {
        status.ok = false;
        status.error = "stage: mdot must be > 0";
        result.status = status;
        return result;
    }
    if (!(omega >= 0.0))
    {
        status.ok = false;
        status.error = "stage: omega must be >= 0";
        result.status = status;
        return result;
    }
    if (!(inlet.p0 > 0.0))
    {
        status.ok = false;
        status.error = "stage: p0 must be > 0";
        result.status = status;
        return result;
    }
    if (!(inlet.T0 > 0.0))
    {
        status.ok = false;
        status.error = "stage: T0 must be > 0";
        result.status = status;
        return result;
    }

    const double cp = eos.specific_heat_cp();
    const double gamma = eos.gamma();
    const double R = eos.gas_constant();
    const double r_mean = config.geometry.r_mean();
    const double area = config.geometry.flow_area();
    const double alpha_1 = config.geometry.alpha_1_rad;
    const double beta_2 = config.geometry.beta_2_rad;
    const double eta_poly = config.loss.polytropic_efficiency;

    // -- Inlet density-velocity fixed point ----------------------------
    // Iterate c_a, static (T1, p1, rho1) with direct substitution; the
    // compressible axial-machine closure converges for axial flow.
    const double rho0 = eos.density(inlet.p0, inlet.T0, status);
    double c_a = mdot / (rho0 * area);
    double c_w1 = 0.0;
    double T1 = inlet.T0;
    double p1 = inlet.p0;
    double rho1 = rho0;
    bool converged = false;
    for (int iter = 0; iter < kMaxFixedPointIterations; ++iter)
    {
        c_w1 = c_a * std::tan(alpha_1);
        const double c1_sq = c_a * c_a + c_w1 * c_w1;
        T1 = inlet.T0 - c1_sq / (2.0 * cp);
        p1 = inlet.p0 * std::pow(T1 / inlet.T0, gamma / (gamma - 1.0));
        rho1 = eos.density(p1, T1, status);
        const double c_a_new = mdot / (rho1 * area);
        if (std::fabs(c_a_new - c_a) / c_a < kFixedPointTolerance)
        {
            c_a = c_a_new;
            converged = true;
            break;
        }
        c_a = c_a_new;
    }
    if (!converged)
    {
        status.ok = false;
        status.error = "stage: density-velocity fixed point did not converge";
        result.status = status;
        return result;
    }

    // -- Euler work (sign is geometry-emergent, no mode flags) ---------
    const double u = omega * r_mean;
    const double w_u1 = c_w1 - u;
    const double w_u2 = c_a * std::tan(beta_2);  // relative exit swirl
    const double c_w2 = u + w_u2;                // absolute exit swirl
    const double delta_c_w = c_w2 - c_w1;
    const double delta_h0 = u * delta_c_w;

    // -- Shaft loads ---------------------------------------------------
    const double torque = mdot * r_mean * delta_c_w;  // momentum at omega = 0 too
    const double power = torque * omega;

    // -- Total-state bookkeeping ---------------------------------------
    const double T02 = inlet.T0 + delta_h0 / cp;
    const double tau = T02 / inlet.T0;
    if (!(tau > 0.0))  // extreme expansion guard
    {
        status.ok = false;
        status.error = "stage: temperature ratio non-positive (extreme expansion)";
        result.status = status;
        return result;
    }

    double pi = 0.0;
    if (delta_h0 > 0.0)
        pi = thermo::polytropic::pressure_ratio_compression(tau, gamma, eta_poly, status);
    else
        pi = thermo::polytropic::pressure_ratio_expansion(tau, gamma, eta_poly, status);
    if (!status.ok)
    {
        result.status = status;
        return result;
    }
    const double p0_out = pi * inlet.p0;

    // -- Static exit ---------------------------------------------------
    const double c2_sq = c_a * c_a + c_w2 * c_w2;
    const double T2 = T02 - c2_sq / (2.0 * cp);
    const double p2 = p0_out * std::pow(T2 / T02, gamma / (gamma - 1.0));
    const double static_rho = eos.density(p2, T2, status);

    // -- Relative Mach at rotor LE (choke/envelope) --------------------
    const double w1_sq = c_a * c_a + w_u1 * w_u1;
    const double a1 = std::sqrt(gamma * R * T1);
    const double mach_rel_le = std::sqrt(w1_sq) / a1;
    result.choked = mach_rel_le >= 1.0;
    if (result.choked)
        status.warnings.push_back("stage: relative Mach >= 1 at rotor LE (choked)");
    else if (mach_rel_le > 0.7)
        status.warnings.push_back("stage: relative Mach exceeds 0.7 validity envelope");

    // -- Non-dimensional coefficients ----------------------------------
    result.work_coefficient = (u != 0.0) ? delta_h0 / (u * u) : 0.0;
    result.flow_coefficient = (u != 0.0) ? c_a / u : 0.0;

    result.p0_out = p0_out;
    result.T0_out = T02;
    result.c_theta_out = c_w2;
    result.static_p = p2;
    result.static_T = T2;
    result.static_rho = static_rho;
    result.delta_h0 = delta_h0;
    result.torque = torque;
    result.power = power;
    result.pi = pi;
    result.tau = tau;
    result.mach_rel_le = mach_rel_le;

    result.status = status;
    result.ok = status.ok;
    return result;
}

} // namespace exd::physics::fluid::turbomachinery