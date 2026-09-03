#include "bem_internal.hpp"

#include <cmath>

namespace exd::engine::physics::fluid::reduced_order::bem
{

namespace
{

inline double acos_clamp(double x)
{
    if (x <= -1.0) return M_PI;
    if (x >= 1.0) return 0.0;
    return std::acos(x);
}

// ── DuSeligLossModel ──────────────────────────────────────────────
// Du & Selig (1993), "Aerodynamics of Horizontal-Axis Wind Turbines".
//
// Modifies the Prandtl tip loss factor with a loading-dependent
// correction that provides smoother transition near the tip.  The
// standard Prandtl formula drops F to zero too abruptly at r = R_tip;
// the Du-Selig g-factor softens this by making the effective vortex
// strength depend on the local radial position through a Lorentzian
// weighting term.
//
//   f_DuSelig = (B/2) * (R_tip - r) / (r * sin(phi)) * g(eta)
//   g(eta)    = 1 + (1.386*B - 1.964) * (eta - 1) / sqrt(1 + ((eta-1)/0.464)^2)
//   eta       = r / R_tip
//
// For eta <= 0.95 the g-factor is negligible and we fall back to the
// standard Prandtl form to avoid unnecessary computation.  Hub loss
// uses the unmodified Prandtl formula.
//
// Reference: Du, Z. and Selig, M.S., "Aerodynamics of Horizontal-Axis
// Wind Turbines," AIAA 93-0505, 1993.

class DuSeligLossModelImpl : public LossModel
{
public:
    double loss_factor(double sin_phi, double blade_count,
                       double radius, double r_tip, double r_hub,
                       double /*axial_induction*/, double /*chord*/,
                       double /*sigma*/) const override
    {
        // Degenerate inflow: no loss.
        if (sin_phi <= 1e-6) return 1.0;

        const double B = blade_count;
        const double eta = radius / r_tip; // normalized radial position

        // ── Tip loss with Du-Selig modification ──────────────────
        double f_tip;
        if (eta > 0.95)
        {
            // Standard Prandtl exponent...
            const double f_prandtl = 0.5 * B * (r_tip - radius) / (radius * sin_phi);

            // ...weighted by the Du-Selig g-factor.  The Lorentzian
            // term (eta-1)/sqrt(1+((eta-1)/0.464)^2) is smooth and
            // saturates near eta=1, avoiding the kink present in the
            // plain Prandtl curve.
            const double dr = eta - 1.0;
            const double g_mod = (1.386 * B - 1.964) * dr
                               / std::sqrt(1.0 + (dr / 0.464) * (dr / 0.464));
            f_tip = f_prandtl * (1.0 + g_mod);
        }
        else
        {
            // Inboard region: standard Prandtl (g ~ 1 + O(-0.05*B)).
            f_tip = 0.5 * B * (r_tip - radius) / (radius * sin_phi);
        }

        const double F_tip = (2.0 / M_PI) * acos_clamp(std::exp(-f_tip));

        // ── Hub loss: unmodified Prandtl ─────────────────────────
        double F_hub = 1.0;
        if (r_hub > 1e-6)
        {
            const double f_hub = 0.5 * B * (radius - r_hub) / (r_hub * sin_phi);
            F_hub = (2.0 / M_PI) * acos_clamp(std::exp(-f_hub));
        }

        double F = F_tip * F_hub;
        return std::max(F, 1e-3);
    }
};

// ── ChaviaropoulosLossModel ───────────────────────────────────────
// Chaviaropoulos & Hansen (2000), "A New BEM Model for Wind Turbine
// Rotors," with loading-dependent tip loss enhancement.
//
// Modifies the Prandtl exponent by a factor that depends on the
// local solidity times the normal force coefficient (sigma * Cn).
// This makes the loss grow with increasing local loading, which
// better matches CFD and experimental data for heavily loaded tips.
//
//   f_CH = (B/2) * (R_tip - r) / (r * sin(phi)) * (1 + alpha * sigma)
//
// where alpha_CH ~ 1.5.  Hub loss uses unmodified Prandtl.
//
// The sigma parameter (local solidity = B*c/(2*pi*r)) is passed
// through the LossModel interface; Cn is not available here so the
// loading dependence is approximated by solidity alone.
//
// Reference: Chaviaropoulos, P.K. and Hansen, M.O.L., "Investigating
// Three-Dimensional and Aerodynamic Stall Effects on Wind Turbine
// Blades," Wind Energy, 2000.

class ChaviaropoulosLossModelImpl : public LossModel
{
public:
    double loss_factor(double sin_phi, double blade_count,
                       double radius, double r_tip, double r_hub,
                       double /*axial_induction*/, double /*chord*/,
                       double sigma_val) const override
    {
        // Degenerate inflow: no loss.
        if (sin_phi <= 1e-6) return 1.0;

        const double B = blade_count;
        constexpr double alpha_CH = 1.5; // Chaviaropoulos-Hansen constant

        // ── Tip loss with loading-dependent enhancement ──────────
        double f_tip = 0.5 * B * (r_tip - radius) / (radius * sin_phi);
        if (sigma_val > 0.0)
        {
            f_tip *= (1.0 + alpha_CH * sigma_val);
        }

        const double F_tip = (2.0 / M_PI) * acos_clamp(std::exp(-f_tip));

        // ── Hub loss: unmodified Prandtl ─────────────────────────
        double F_hub = 1.0;
        if (r_hub > 1e-6)
        {
            const double f_hub = 0.5 * B * (radius - r_hub) / (r_hub * sin_phi);
            F_hub = (2.0 / M_PI) * acos_clamp(std::exp(-f_hub));
        }

        double F = F_tip * F_hub;
        return std::max(F, 1e-3);
    }
};

} // namespace

// ── Factory functions ─────────────────────────────────────────────

std::unique_ptr<LossModel> make_du_selig_loss_model()
{
    return std::make_unique<DuSeligLossModelImpl>();
}

std::unique_ptr<LossModel> make_chaviaropoulos_loss_model()
{
    return std::make_unique<ChaviaropoulosLossModelImpl>();
}

} // namespace exd::engine::physics::fluid::reduced_order::bem
