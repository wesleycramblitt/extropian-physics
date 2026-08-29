// momentum_balance.cpp
// Reduced-order blade-element momentum balance. Mirrors the standard
// induction solve (induction.cpp) and the element force loop
// (bem_solver.cpp) of the BEM turbine solver: same iterative induction
// with Buhl/Glauert high-induction branch, same clamps, same dL/dL/dQ/dT
// formulas and the same sign convention (thrust into the inflow,
// torque in the direction of positive omega).

#include <exd/physics/fluid/forces/force_evaluator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace exd::physics::fluid::forces
{

namespace
{

constexpr double kNearZero = 1e-12;
constexpr double kGlauertThreshold = 0.4; // Buhl/Glauert high-induction guard
constexpr double kSinPhiDegenerate = 1e-3;
constexpr double kTangentialDenomGuard = 1e-12;
constexpr char kFallbackAirfoil[] = "naca0012";

bool normalize3(std::array<double, 3>& v)
{
    const double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len < kNearZero) return false;
    const double inv = 1.0 / len;
    v[0] *= inv;
    v[1] *= inv;
    v[2] *= inv;
    return true;
}

/// Resolve (cl, cd) for one station; falls back to "naca0012" with a
/// single warning when `airfoil` is absent from the database.
std::pair<double, double> resolve_coefficients(const PolarDatabase& db,
                                               const std::string& airfoil,
                                               double alpha_deg, double re,
                                               bool& warned_fallback,
                                               mechanics::ModelStatus& status)
{
    const auto* polar = db.find(airfoil, re);
    if (!polar)
    {
        if (!warned_fallback)
        {
            status.warnings.push_back("airfoil '" + airfoil + "' not found in "
                                     "polar database; using 'naca0012'");
            warned_fallback = true;
        }
        polar = db.find(kFallbackAirfoil, re);
    }
    if (!polar) return {0.0, 0.0};
    return polar->evaluate(alpha_deg);
}

class MomentumBalanceEvaluator final : public IForceEvaluator
{
public:
    MomentumBalanceEvaluator(MomentumBalanceConfig config, const PolarDatabase* polars)
        : config_(config)
    {
        if (polars)
        {
            polars_ = polars;
        }
        else
        {
            fallback_.add_builtin_polars();
            polars_ = &fallback_;
        }
    }

    std::string_view name() const override { return "momentum_balance"; }

    void compute(const BladeGeometry& blade,
                 const SurfaceFlow& flow,
                 double omega,
                 const mechanics::RotationAxis& axis,
                 std::vector<mechanics::ElementForce3D>& per_element,
                 mechanics::ModelStatus& status) const override
    {
        per_element.clear();
        if (!flow.valid())
        {
            status.ok = false;
            status.error = "invalid SurfaceFlow";
            return;
        }

        // Uniform inflow: the first sampled velocity, projected on the axis.
        std::array<double, 3> e_z = axis.direction;
        if (!normalize3(e_z))
        {
            status.ok = false;
            status.error = "axis direction is (near) zero";
            return;
        }
        const double v_ax = flow.velocity[0][0] * e_z[0] +
                            flow.velocity[0][1] * e_z[1] +
                            flow.velocity[0][2] * e_z[2];
        const double v_inf = std::fabs(v_ax); // incoming axial speed
        if (v_inf < kNearZero)
        {
            status.ok = false;
            status.error = "zero axial inflow";
            return;
        }

        const int blade_count = blade.blade_count;

        bool warned_fallback = false;
        for (std::size_t st_i = 0; st_i < blade.stations.size(); ++st_i)
        {
            const BladeStation& st = blade.stations[st_i];
            const double r = std::max(st.r, 1e-9);
            const double sigma = static_cast<double>(blade_count) * st.chord /
                                 (2.0 * M_PI * r);

            // ── Induction solve (mirror StandardInductionModel, F = 1) ──
            double a = 0.0;
            double a_prime = 0.0;
            double phi = 0.0;
            double W = 0.0;
            double cl = 0.0;
            double cd = 0.0;
            bool converged = false;

            for (int iter = 1; iter <= config_.max_iterations; ++iter)
            {
                const double va = v_inf * (1.0 - a);
                const double vt = omega * r * (1.0 + a_prime);
                W = std::sqrt(va * va + vt * vt);
                phi = std::atan2(va, vt);
                const double sin_phi = std::sin(phi);

                // Degenerate near-axial inflow: no blade force.
                if (std::fabs(sin_phi) < kSinPhiDegenerate)
                {
                    a = 0.0;
                    a_prime = 0.0;
                    cl = 0.0;
                    cd = 0.0;
                    converged = true;
                    break;
                }

                const double alpha_deg = phi * 180.0 / M_PI - st.twist_deg;
                const double re = flow.density * W * st.chord / flow.viscosity;
                const auto [clv, cdv] = resolve_coefficients(*polars_, st.airfoil,
                                                             alpha_deg, re,
                                                             warned_fallback, status);
                cl = clv;
                cd = cdv;

                const double cos_phi = std::cos(phi);
                const double cn = cl * cos_phi + cd * sin_phi;
                const double ct = cl * sin_phi - cd * cos_phi;

                // Buhl high-induction branch (loss factor F = 1 here).
                constexpr double F = 1.0;
                const double c_t = sigma * cn * (1.0 - a) * (1.0 - a) /
                                   (F * sin_phi * sin_phi);
                const double rad = c_t * (50.0 - 36.0 * F) + 12.0 * F * (3.0 * F - 4.0);
                const bool use_buhl = (a > kGlauertThreshold) && (rad >= 0.0);

                double a_new;
                if (use_buhl)
                {
                    a_new = (18.0 * F - 20.0 - 3.0 * std::sqrt(rad)) /
                            (36.0 * F - 50.0);
                }
                else
                {
                    const double denom = 4.0 * F * sin_phi * sin_phi + sigma * cn;
                    a_new = sigma * cn / denom;
                }
                if (a > kGlauertThreshold && rad < 0.0)
                {
                    status.warnings.push_back(
                        "Buhl discriminant < 0; low-induction branch used");
                }

                // Tangential induction.
                const double denom_ap = 4.0 * F * sin_phi * cos_phi - sigma * ct;
                double a_prime_new;
                if (denom_ap <= kTangentialDenomGuard)
                {
                    status.warnings.push_back("tangential denom <= 0");
                    a_prime_new = (ct > 0.0) ? 1.0 : 0.0;
                }
                else
                {
                    a_prime_new = sigma * ct / denom_ap;
                }

                // Clamps (same as induction.cpp).
                if (a_new < 0.0)
                {
                    status.warnings.push_back("axial induction clamped to 0");
                    a_new = 0.0;
                }
                if (a_new >= 1.0)
                {
                    status.warnings.push_back("axial induction clamped to <1");
                    a_new = 0.999999;
                }
                if (a_prime_new < 0.0)
                {
                    status.warnings.push_back("tangential induction clamped to 0");
                    a_prime_new = 0.0;
                }
                if (a_prime_new > 1.0)
                {
                    status.warnings.push_back("tangential induction clamped to 1");
                    a_prime_new = 1.0;
                }

                // Under-relaxation.
                const double da = a_new - a;
                const double dap = a_prime_new - a_prime;
                a += config_.under_relaxation * da;
                a_prime += config_.under_relaxation * dap;

                // Convergence check (tolerance relative to current induction).
                const double ref_a = std::max(a, 1e-2);
                const double ref_ap = std::max(a_prime, 1e-2);
                if (std::fabs(da) < config_.tolerance * ref_a &&
                    std::fabs(dap) < config_.tolerance * ref_ap)
                {
                    converged = true;
                    break;
                }
            }
            if (!converged)
            {
                status.warnings.push_back("momentum balance: induction not "
                                          "converged at r = " + std::to_string(st.r));
            }

            // ── Element forces (mirror bem_solver.cpp) ──
            const double sin_phi = std::sin(phi);
            const double cos_phi = std::cos(phi);
            const double dL = 0.5 * flow.density * W * W * st.chord * cl * st.dr;
            const double dD = 0.5 * flow.density * W * W * st.chord * cd * st.dr;
            const double dQ = static_cast<double>(blade_count) * r *
                              (dL * sin_phi - dD * cos_phi);
            const double dT = static_cast<double>(blade_count) *
                              (dL * cos_phi + dD * sin_phi);

            // Per-blade 3D loads distributed over the rotating blades.
            const double per_blade_axial = dT / static_cast<double>(blade_count);
            const double per_blade_tangential = dQ / (r * static_cast<double>(blade_count));
            for (int b = 0; b < blade_count; ++b)
            {
                const double theta = static_cast<double>(b) * 2.0 * M_PI /
                                     static_cast<double>(blade_count);
                const RotorFrame frame = make_rotor_frame(theta, axis, status);

                mechanics::ElementForce3D out;
                out.r = st.r;
                out.ref = {
                    axis.origin[0] + blade.z_rotor * frame.e_z[0] + st.r * frame.e_r[0],
                    axis.origin[1] + blade.z_rotor * frame.e_z[1] + st.r * frame.e_r[1],
                    axis.origin[2] + blade.z_rotor * frame.e_z[2] + st.r * frame.e_r[2],
                };
                // Thrust acts into the inflow (downwind): -dT on +e_z when the
                // inflow is -e_z. Torque follows the direction of omega.
                out.force = {
                    -per_blade_axial * frame.e_z[0] + per_blade_tangential * frame.e_t[0],
                    -per_blade_axial * frame.e_z[1] + per_blade_tangential * frame.e_t[1],
                    -per_blade_axial * frame.e_z[2] + per_blade_tangential * frame.e_t[2],
                };
                out.force_pressure = out.force; // diagnostic only, not a physical split
                out.force_shear = {0, 0, 0};
                out.moment = {0, 0, 0};
                per_element.push_back(out);
            }
        }
    }

private:
    MomentumBalanceConfig config_;
    const PolarDatabase* polars_ = nullptr;
    PolarDatabase fallback_;
};

} // anonymous namespace

std::unique_ptr<IForceEvaluator> make_momentum_balance_evaluator(
    const MomentumBalanceConfig& config, const PolarDatabase* polars)
{
    return std::make_unique<MomentumBalanceEvaluator>(config, polars);
}

} // namespace exd::physics::fluid::forces