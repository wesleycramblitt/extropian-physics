// reactor.cpp
// 0D well-stirred constant-volume reactor with generic mass-action chemistry
// and optional Arrhenius temperature scaling, integrated with the shared
// solver::integrate_step() integrator.

#include <exd/physics/chemistry/reactor.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace exd::physics::chemistry {

namespace {

constexpr double kGasConstant = 8.314; // J/(mol*K)

/// Upper bound on the number of recorded history snapshots.
constexpr std::size_t kMaxHistory = 20000;

} // namespace

bool validate_chemistry_config(const ChemistryConfig& config,
                               std::string& error,
                               std::vector<std::string>& warnings)
{
    error.clear();
    warnings.clear();

    if (config.species.empty())
    {
        error = "chemistry: at least one species is required";
        return false;
    }
    for (std::size_t i = 0; i < config.species.size(); ++i)
    {
        if (config.species[i].initial_concentration < 0.0)
        {
            error = "chemistry: initial_concentration must be >= 0";
            return false;
        }
    }

    const std::size_t ns = config.species.size();
    for (std::size_t r = 0; r < config.reactions.size(); ++r)
    {
        const ReactionSpec& re = config.reactions[r];
        for (const int32_t idx : re.reactant_indices)
        {
            if (idx < 0 || static_cast<std::size_t>(idx) >= ns)
            {
                error = "chemistry: reaction reactant index out of range";
                return false;
            }
        }
        for (const int32_t idx : re.product_indices)
        {
            if (idx < 0 || static_cast<std::size_t>(idx) >= ns)
            {
                error = "chemistry: reaction product index out of range";
                return false;
            }
        }
        for (const double s : re.reactant_stoich)
        {
            if (!(s > 0.0))
            {
                error = "chemistry: reactant stoichimetric coefficients must be > 0";
                return false;
            }
        }
        for (const double s : re.product_stoich)
        {
            if (!(s > 0.0))
            {
                error = "chemistry: product stoichimetric coefficients must be > 0";
                return false;
            }
        }
        if (re.rate_constant < 0.0)
        {
            error = "chemistry: rate_constant must be >= 0";
            return false;
        }
        if (re.activation_energy < 0.0)
        {
            error = "chemistry: activation_energy must be >= 0";
            return false;
        }
        if (re.pre_exponential < 0.0)
        {
            error = "chemistry: pre_exponential must be >= 0";
            return false;
        }
        if (re.activation_energy > 0.0 && re.pre_exponential == 0.0)
            warnings.push_back("chemistry: reaction " + std::to_string(r)
                               + " uses Arrhenius scaling with zero pre-exponential (rate is zero)");
        if (re.activation_energy == 0.0 && re.rate_constant == 0.0)
            warnings.push_back("chemistry: reaction " + std::to_string(r)
                               + " has zero rate constant (inert)");
    }

    if (!(config.temperature > 0.0))
    {
        error = "chemistry: temperature must be > 0";
        return false;
    }
    if (!(config.dt > 0.0))
    {
        error = "chemistry: dt must be > 0";
        return false;
    }
    if (!(config.end_time > 0.0))
    {
        error = "chemistry: end_time must be > 0";
        return false;
    }
    return true;
}

ChemistryResult solve_chemistry(const ChemistryConfig& config, ModelStatus& status)
{
    ChemistryResult res;
    std::string error;
    std::vector<std::string> warnings;
    if (!validate_chemistry_config(config, error, warnings))
    {
        res.status = ModelStatus{false, error, warnings};
        status = res.status;
        return res;
    }

    const std::size_t ns = config.species.size();

    // ── effective rate constants (Arrhenius when enabled) ──
    std::vector<double> keff(config.reactions.size(), 0.0);
    for (std::size_t r = 0; r < config.reactions.size(); ++r)
    {
        const ReactionSpec& re = config.reactions[r];
        if (re.activation_energy > 0.0)
        {
            keff[r] = re.pre_exponential
                      * std::exp(-re.activation_energy
                                 / (kGasConstant * config.temperature));
        }
        else
        {
            keff[r] = re.rate_constant;
        }
    }

    std::vector<double> c(ns, 0.0);
    for (std::size_t i = 0; i < ns; ++i) c[i] = config.species[i].initial_concentration;
    for (std::size_t i = 0; i < ns; ++i)
    {
        if (!std::isfinite(c[i]))
        {
            res.status = ModelStatus{false, "chemistry: initial concentration must be finite", warnings};
            status = res.status;
            return res;
        }
    }

    solver::DerivativeFn deriv = [&](std::span<const double> st,
                                     std::span<double> dst, double)
    {
        for (std::size_t i = 0; i < ns; ++i) dst[i] = 0.0;
        for (std::size_t r = 0; r < config.reactions.size(); ++r)
        {
            const ReactionSpec& re = config.reactions[r];
            double rate = keff[r];
            for (std::size_t q = 0; q < re.reactant_indices.size(); ++q)
            {
                const std::size_t idx = static_cast<std::size_t>(re.reactant_indices[q]);
                rate *= std::pow(std::max(st[idx], 0.0), re.reactant_stoich[q]);
            }
            for (std::size_t q = 0; q < re.product_indices.size(); ++q)
            {
                const std::size_t idx = static_cast<std::size_t>(re.product_indices[q]);
                dst[idx] += re.product_stoich[q] * rate;
            }
            for (std::size_t q = 0; q < re.reactant_indices.size(); ++q)
            {
                const std::size_t idx = static_cast<std::size_t>(re.reactant_indices[q]);
                dst[idx] -= re.reactant_stoich[q] * rate;
            }
        }
    };

    // ── history recording, capped at kMaxHistory snapshots ──
    const std::uint64_t total_steps = std::max<std::uint64_t>(
        1, static_cast<std::uint64_t>(std::ceil(config.end_time / config.dt)));
    const std::uint64_t record_interval = std::max<std::uint64_t>(
        1, total_steps / static_cast<std::uint64_t>(kMaxHistory));
    const double clamp_tol = 1e-6;

    auto record = [&](double t)
    {
        res.history.push_back(c);
        res.time_history.push_back(t);
    };

    double t = 0.0;
    record(t);
    double max_mass_err = 0.0;
    for (std::uint64_t step = 0; step < total_steps; ++step)
    {
        if (t >= config.end_time - 1e-12) break;
        const double step_dt = std::min(config.dt, config.end_time - t);
        const std::vector<double> c_old = c;
        exd::physics::ModelStatus ist;
        double dt_used = 0.0;
        const bool ok = solver::integrate_step(config.integration, t, step_dt,
                                               std::span<double>(c), deriv,
                                               ist, &dt_used);
        if (!ok || !ist.ok)
        {
            const std::string err = "chemistry: integration failed: "
                                    + (ist.error.empty() ? "integrator step rejected"
                                                         : ist.error);
            res.status = ModelStatus{false, err, warnings};
            status = res.status;
            return res;
        }
        t += (dt_used > 0.0) ? dt_used : step_dt;

        // ── non-negativity enforcement ──
        double max_neg = 0.0;
        for (std::size_t i = 0; i < ns; ++i)
        {
            if (c[i] < -max_neg) max_neg = -c[i];
        }
        if (max_neg > 0.0)
        {
            if (max_neg > clamp_tol)
            {
                res.status = ModelStatus{
                    false,
                    "chemistry: species went strongly negative (max -value "
                        + std::to_string(max_neg) + " > clamp tolerance)",
                    warnings};
                status = res.status;
                return res;
            }
            warnings.push_back("chemistry: species went slightly negative; clamped to 0");
            for (std::size_t i = 0; i < ns; ++i)
            {
                if (c[i] < 0.0) c[i] = 0.0;
            }
        }

        // ── total-mass bookkeeping ──
        double sum_old = 0.0;
        double sum_new = 0.0;
        for (std::size_t i = 0; i < ns; ++i)
        {
            sum_old += c_old[i];
            sum_new += c[i];
        }
        const double me = std::fabs(sum_new - sum_old) / (1.0 + sum_old);
        if (me > max_mass_err) max_mass_err = me;

        if ((step + 1) % record_interval == 0) record(t);
    }
    // Ensure the final state is recorded exactly once.
    if (res.time_history.empty() || res.time_history.back() != t) record(t);

    res.final_concentrations = c;
    res.max_total_mass_error = max_mass_err;
    res.ok = true;
    res.status = ModelStatus{true, "", warnings};
    status = res.status;
    return res;
}

} // namespace exd::physics::chemistry