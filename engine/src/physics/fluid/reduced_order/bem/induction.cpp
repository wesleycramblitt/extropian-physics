#include "bem_internal.hpp"

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::reduced_order::bem
{

namespace
{

inline double deg2rad(double d) { return d * M_PI / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / M_PI; }

class StandardInductionModel : public InductionModel
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

            const double C_T = sigma * cn * (1.0 - a) * (1.0 - a) / (F * sin_phi * sin_phi);
            const double rad = C_T * (50.0 - 36.0 * F) + 12.0 * F * (3.0 * F - 4.0);
            const bool use_buhl = (a > threshold) && (rad >= 0.0);

            double a_new;
            if (use_buhl)
            {
                a_new = (18.0 * F - 20.0 - 3.0 * std::sqrt(rad)) / (36.0 * F - 50.0);
            }
            else
            {
                const double denom = 4.0 * F * sin_phi * sin_phi + sigma * cn;
                a_new = sigma * cn / denom;
            }

            if (a > threshold && rad < 0.0)
            {
                warnings.push_back("Buhl discriminant < 0; low-induction branch used");
            }

            // Tangential induction.
            double denom_ap = 4.0 * F * sin_phi * cos_phi - sigma * ct;
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

            const double da = a_new - a;
            const double dap = a_prime_new - a_prime;

            a = a + lambda * da;
            a_prime = a_prime + lambda * dap;

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

std::unique_ptr<InductionModel> make_standard_induction_model()
{
    return std::make_unique<StandardInductionModel>();
}

} // namespace exd::engine::physics::fluid::reduced_order::bem
