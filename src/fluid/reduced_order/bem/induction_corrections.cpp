#include "bem_internal.hpp"

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace exd::physics::fluid::reduced_order::bem
{

namespace
{

inline double deg2rad(double d) { return d * M_PI / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / M_PI; }

// ── GlauertIterativeInductionModel ────────────────────────────────
// Iterative application of Glauert's empirical high-induction correction.
// For a <= threshold: momentum theory (same as standard).
// For a > threshold:  Glauert empirical CT curve solved via fixed-point
// iteration with under-relaxation.
//
// Reference: Glauert (1935), implemented per NREL AeroDyn convention.

class GlauertIterativeInductionModel : public InductionModel
{
public:
    StationState solve(const BladeElementInput& elem,
                       double v_rotor, double omega,
                       const OperatingConditions& conditions,
                       const PolarDatabase& polars,
                       const BEMSolverConfig& config,
                       const LossModel& loss_model,
                       std::vector<std::string>& warnings) const override
    {
        StationState s;
        const double tol = config.induction_tolerance;
        const double lambda = config.under_relaxation;
        const double threshold = config.glauert_threshold;
        const double sigma = elem.blade_count * elem.chord / (2.0 * M_PI * elem.r);

        double a = 0.0;
        double a_prime = 0.0;

        for (uint32_t iter = 1; iter <= config.max_iterations; ++iter)
        {
            const double Va = v_rotor * (1.0 - a);
            const double Vt = omega * elem.r * (1.0 + a_prime);
            const double W = std::sqrt(Va * Va + Vt * Vt);
            const double phi = std::atan2(Va, Vt);
            const double sin_phi = std::sin(phi);

            if (std::fabs(sin_phi) < 1e-3)
            {
                s.a = 0.0;
                s.a_prime = 0.0;
                s.phi_rad = phi;
                s.alpha_deg = rad2deg(phi) - elem.beta_deg;
                s.W = W;
                s.Va = Va;
                s.cl = 0.0;
                s.cd = 0.0;
                s.Re = conditions.rho * W * elem.chord / conditions.mu;
                s.F = 1.0;
                s.converged = true;
                s.iterations = 0;
                return s;
            }

            const double alpha_deg = rad2deg(phi) - elem.beta_deg;
            const double Re = conditions.rho * W * elem.chord / conditions.mu;
            const auto* polar = polars.find(elem.airfoil, Re);
            double cl = 0.0, cd = 0.0;
            if (polar)
            {
                auto [clv, cdv] = polar->evaluate(alpha_deg);
                cl = clv;
                cd = cdv;
            }
            else
            {
                warnings.push_back("airfoil '" + elem.airfoil + "' not found in polar database");
            }

            const double cos_phi = std::cos(phi);
            const double cn = cl * cos_phi + cd * sin_phi;
            const double ct = cl * sin_phi - cd * cos_phi;

            const double F = loss_model.loss_factor(sin_phi, elem.blade_count, elem.r,
                                                    elem.r_tip, elem.r_hub, a, elem.chord, sigma);

            // Axial induction factor.
            double a_new;
            if (a <= threshold)
            {
                // Low-induction: standard momentum theory.
                const double denom = 4.0 * F * sin_phi * sin_phi + sigma * cn;
                a_new = sigma * cn / denom;
            }
            else
            {
                // High-induction: Glauert empirical correction, iterated.
                const double sin_phi2 = sin_phi * sin_phi;
                const double CT_blade = sigma * cn * (1.0 - a) * (1.0 - a)
                                        / (F * sin_phi2);

                if (CT_blade < 0.889)
                {
                    // Below Glauert's empirical limit: direct mapping.
                    a_new = CT_blade / (4.0 * F);
                }
                else
                {
                    // Glauert (1935) empirical curve inverted:
                    // CT = 4a(1 - 0.25(5-3a)a) -> solve for a.
                    const double disc = 1.0 - CT_blade / F;
                    a_new = 0.5 * (1.0 + std::sqrt(std::max(0.0, disc)));
                }
            }

            // Tangential induction factor.
            const double denom_ap = 4.0 * F * sin_phi * cos_phi - sigma * ct;
            double a_prime_new;
            if (denom_ap <= 1e-12)
            {
                warnings.push_back("tangential denom <= 0");
                a_prime_new = (ct > 0.0) ? 1.0 : 0.0;
            }
            else
            {
                a_prime_new = sigma * ct / denom_ap;
            }

            // Clamps.
            if (a_new < 0.0)
            {
                warnings.push_back("axial induction clamped to 0");
                a_new = 0.0;
            }
            if (a_new >= 1.0)
            {
                warnings.push_back("axial induction clamped to <1");
                a_new = 0.999999;
            }
            if (a_prime_new < 0.0)
            {
                warnings.push_back("tangential induction clamped to 0");
                a_prime_new = 0.0;
            }
            if (a_prime_new > 1.0)
            {
                warnings.push_back("tangential induction clamped to 1");
                a_prime_new = 1.0;
            }

            // Under-relaxation.
            const double da = a_new - a;
            const double dap = a_prime_new - a_prime;

            a = a + lambda * da;
            a_prime = a_prime + lambda * dap;

            // Convergence check.
            const double abs_da = std::fabs(da);
            const double abs_dap = std::fabs(dap);
            const double ref_a = std::max(a, 1e-2);
            const double ref_ap = std::max(a_prime, 1e-2);

            if (abs_da < tol * ref_a && abs_dap < tol * ref_ap)
            {
                s.a = a;
                s.a_prime = a_prime;
                s.phi_rad = phi;
                s.alpha_deg = alpha_deg;
                s.W = W;
                s.Va = Va;
                s.cl = cl;
                s.cd = cd;
                s.Re = Re;
                s.F = F;
                s.converged = true;
                s.iterations = iter;
                return s;
            }
        }

        warnings.push_back("induction did not converge at radius " + std::to_string(elem.r));
        s.a = a;
        s.a_prime = a_prime;
        s.converged = false;
        s.iterations = config.max_iterations;
        return s;
    }
};

// ── SnelInductionModel ───────────────────────────────────────────
// Snel (1993) exponential blending between momentum theory and the
// empirical high-induction regime.  Smooth transition, no hard threshold.
//
// a_corrected = a_mom * (1 - exp(-4 F sin^2(phi) / (sigma Cn) * a_mom))
//
// For small a_mom the exponential ~ 1 - (...), so a_corrected ~ a_mom.
// For large a_mom the exponential decays, limiting growth.

class SnelInductionModel : public InductionModel
{
public:
    StationState solve(const BladeElementInput& elem,
                       double v_rotor, double omega,
                       const OperatingConditions& conditions,
                       const PolarDatabase& polars,
                       const BEMSolverConfig& config,
                       const LossModel& loss_model,
                       std::vector<std::string>& warnings) const override
    {
        StationState s;
        const double tol = config.induction_tolerance;
        const double lambda = config.under_relaxation;
        const double sigma = elem.blade_count * elem.chord / (2.0 * M_PI * elem.r);

        double a = 0.0;
        double a_prime = 0.0;

        for (uint32_t iter = 1; iter <= config.max_iterations; ++iter)
        {
            const double Va = v_rotor * (1.0 - a);
            const double Vt = omega * elem.r * (1.0 + a_prime);
            const double W = std::sqrt(Va * Va + Vt * Vt);
            const double phi = std::atan2(Va, Vt);
            const double sin_phi = std::sin(phi);

            if (std::fabs(sin_phi) < 1e-3)
            {
                s.a = 0.0;
                s.a_prime = 0.0;
                s.phi_rad = phi;
                s.alpha_deg = rad2deg(phi) - elem.beta_deg;
                s.W = W;
                s.Va = Va;
                s.cl = 0.0;
                s.cd = 0.0;
                s.Re = conditions.rho * W * elem.chord / conditions.mu;
                s.F = 1.0;
                s.converged = true;
                s.iterations = 0;
                return s;
            }

            const double alpha_deg = rad2deg(phi) - elem.beta_deg;
            const double Re = conditions.rho * W * elem.chord / conditions.mu;
            const auto* polar = polars.find(elem.airfoil, Re);
            double cl = 0.0, cd = 0.0;
            if (polar)
            {
                auto [clv, cdv] = polar->evaluate(alpha_deg);
                cl = clv;
                cd = cdv;
            }
            else
            {
                warnings.push_back("airfoil '" + elem.airfoil + "' not found in polar database");
            }

            const double cos_phi = std::cos(phi);
            const double cn = cl * cos_phi + cd * sin_phi;
            const double ct = cl * sin_phi - cd * cos_phi;

            const double F = loss_model.loss_factor(sin_phi, elem.blade_count, elem.r,
                                                    elem.r_tip, elem.r_hub, a, elem.chord, sigma);

            // Standard momentum theory axial induction.
            const double sin_phi2 = sin_phi * sin_phi;
            const double denom_mom = 4.0 * F * sin_phi2 + sigma * cn;
            const double a_mom = sigma * cn / denom_mom;

            // Snel exponential blending correction.
            double a_new;
            if (sigma * cn > 1e-12)
            {
                const double exponent = -4.0 * F * sin_phi2 / (sigma * cn);
                a_new = a_mom * (1.0 - std::exp(exponent * a_mom));
            }
            else
            {
                a_new = 0.0;
            }

            // Tangential induction factor.
            const double denom_ap = 4.0 * F * sin_phi * cos_phi - sigma * ct;
            double a_prime_new;
            if (denom_ap <= 1e-12)
            {
                warnings.push_back("tangential denom <= 0");
                a_prime_new = (ct > 0.0) ? 1.0 : 0.0;
            }
            else
            {
                a_prime_new = sigma * ct / denom_ap;
            }

            // Clamps.
            if (a_new < 0.0)
            {
                warnings.push_back("axial induction clamped to 0");
                a_new = 0.0;
            }
            if (a_new >= 1.0)
            {
                warnings.push_back("axial induction clamped to <1");
                a_new = 0.999999;
            }
            if (a_prime_new < 0.0)
            {
                warnings.push_back("tangential induction clamped to 0");
                a_prime_new = 0.0;
            }
            if (a_prime_new > 1.0)
            {
                warnings.push_back("tangential induction clamped to 1");
                a_prime_new = 1.0;
            }

            // Under-relaxation.
            const double da = a_new - a;
            const double dap = a_prime_new - a_prime;

            a = a + lambda * da;
            a_prime = a_prime + lambda * dap;

            // Convergence check.
            const double abs_da = std::fabs(da);
            const double abs_dap = std::fabs(dap);
            const double ref_a = std::max(a, 1e-2);
            const double ref_ap = std::max(a_prime, 1e-2);

            if (abs_da < tol * ref_a && abs_dap < tol * ref_ap)
            {
                s.a = a;
                s.a_prime = a_prime;
                s.phi_rad = phi;
                s.alpha_deg = alpha_deg;
                s.W = W;
                s.Va = Va;
                s.cl = cl;
                s.cd = cd;
                s.Re = Re;
                s.F = F;
                s.converged = true;
                s.iterations = iter;
                return s;
            }
        }

        warnings.push_back("induction did not converge at radius " + std::to_string(elem.r));
        s.a = a;
        s.a_prime = a_prime;
        s.converged = false;
        s.iterations = config.max_iterations;
        return s;
    }
};

} // namespace

// ── Factory functions ─────────────────────────────────────────────

std::unique_ptr<InductionModel> make_glauert_iterative_induction_model()
{
    return std::make_unique<GlauertIterativeInductionModel>();
}

std::unique_ptr<InductionModel> make_snel_induction_model()
{
    return std::make_unique<SnelInductionModel>();
}

std::unique_ptr<InductionModel> make_induction_model(InductionCorrection type)
{
    switch (type)
    {
        case InductionCorrection::Standard:         return make_standard_induction_model();
        case InductionCorrection::GlauertIterative: return make_glauert_iterative_induction_model();
        case InductionCorrection::Snel:             return make_snel_induction_model();
    }
    return make_standard_induction_model();
}

} // namespace exd::physics::fluid::reduced_order::bem
