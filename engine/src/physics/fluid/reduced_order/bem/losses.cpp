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

/// Classic Prandtl tip/hub loss factor (§6.3 of architecture doc).
/// F = F_tip * F_hub, where each factor = (2/pi)*acos(exp(-f)).
/// The additional parameters (axial_induction, chord, sigma) are unused by
/// this model but present for interface compatibility.
class PrandtlLossModelImpl : public LossModel
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
        const double f_tip = 0.5 * B * (r_tip - radius) / (radius * sin_phi);
        const double f_hub = 0.5 * B * (radius - r_hub) / (r_hub * sin_phi);

        const double F_tip = (2.0 / M_PI) * acos_clamp(std::exp(-f_tip));
        double F_hub = 1.0;
        if (r_hub > 1e-6)
            F_hub = (2.0 / M_PI) * acos_clamp(std::exp(-f_hub));

        double F = F_tip * F_hub;
        if (F < 1e-3) F = 1e-3;
        return F;
    }
};

} // namespace

std::unique_ptr<LossModel> make_prandtl_loss_model()
{
    return std::make_unique<PrandtlLossModelImpl>();
}

// ── Dispatcher ─────────────────────────────────────────────────────
std::unique_ptr<LossModel> make_loss_model(LossCorrection type)
{
    switch (type)
    {
        case LossCorrection::Prandtl:         return make_prandtl_loss_model();
        case LossCorrection::DuSelig:         return make_du_selig_loss_model();
        case LossCorrection::Chaviaropoulos:  return make_chaviaropoulos_loss_model();
    }
    return make_prandtl_loss_model(); // fallback
}

} // namespace exd::engine::physics::fluid::reduced_order::bem
