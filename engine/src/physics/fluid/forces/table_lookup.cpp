// table_lookup.cpp
// Coefficient-based blade-element force without momentum coupling:
// same local-flow geometry and force emissions as momentum_balance.cpp,
// but axial/tangential induction is fixed at zero (one-shot local flow).

#include <exd/engine/physics/fluid/forces/force_evaluator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace exd::engine::physics::fluid::forces
{

namespace
{

constexpr double kNearZero = 1e-12;
constexpr double kSinPhiDegenerate = 1e-3;
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

std::pair<double, double> resolve_coefficients(const PolarDatabase& db,
                                               const std::string& airfoil,
                                               double alpha_deg, double re,
                                               bool& warned_fallback,
                                               exd::engine::physics::rigid_body::ModelStatus& status)
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

class TableLookupEvaluator final : public IForceEvaluator
{
public:
    TableLookupEvaluator(TableLookupConfig config, const PolarDatabase* polars)
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

    std::string_view name() const override { return "table_lookup"; }

    void compute(const BladeGeometry& blade,
                 const SurfaceFlow& flow,
                 double omega,
                 const exd::engine::physics::rigid_body::RotationAxis& axis,
                 std::vector<exd::engine::physics::rigid_body::ElementForce3D>& per_element,
                 exd::engine::physics::rigid_body::ModelStatus& status) const override
    {
        per_element.clear();
        if (!flow.valid())
        {
            status.ok = false;
            status.error = "invalid SurfaceFlow";
            return;
        }

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
        const double v_inf = std::fabs(v_ax);
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

            // Induction is fixed at zero; the local flow is a one-shot
            // fixed point that converges on the first pass.
            const double a = 0.0;
            const double a_prime = 0.0;
            double phi = 0.0;
            double W = 0.0;
            double cl = 0.0;
            double cd = 0.0;

            for (int iter = 1; iter <= config_.max_iterations; ++iter)
            {
                const double va = v_inf * (1.0 - a);
                const double vt = omega * r * (1.0 + a_prime);
                W = std::sqrt(va * va + vt * vt);
                phi = std::atan2(va, vt);
                const double sin_phi = std::sin(phi);

                if (std::fabs(sin_phi) < kSinPhiDegenerate)
                {
                    cl = 0.0;
                    cd = 0.0;
                    break;
                }

                const double alpha_deg = phi * 180.0 / M_PI - st.twist_deg;
                const double re = flow.density * W * st.chord / flow.viscosity;
                const auto [clv, cdv] = resolve_coefficients(*polars_, st.airfoil,
                                                             alpha_deg, re,
                                                             warned_fallback, status);
                cl = clv;
                cd = cdv;

                // Fixed induction: the axial/tangential momentum residuals are
                // identically zero, so the loop converges on the first pass.
                const double ref_a = std::max(a, 1e-2);
                const double ref_ap = std::max(a_prime, 1e-2);
                if (std::fabs(0.0 - a) < config_.tolerance * ref_a &&
                    std::fabs(0.0 - a_prime) < config_.tolerance * ref_ap)
                {
                    break;
                }
            }

            const double sin_phi = std::sin(phi);
            const double cos_phi = std::cos(phi);
            const double dL = 0.5 * flow.density * W * W * st.chord * cl * st.dr;
            const double dD = 0.5 * flow.density * W * W * st.chord * cd * st.dr;
            const double dQ = static_cast<double>(blade_count) * r *
                              (dL * sin_phi - dD * cos_phi);
            const double dT = static_cast<double>(blade_count) *
                              (dL * cos_phi + dD * sin_phi);

            const double per_blade_axial = dT / static_cast<double>(blade_count);
            const double per_blade_tangential = dQ / (r * static_cast<double>(blade_count));
            for (int b = 0; b < blade_count; ++b)
            {
                const double theta = static_cast<double>(b) * 2.0 * M_PI /
                                     static_cast<double>(blade_count);
                const RotorFrame frame = make_rotor_frame(theta, axis, status);

                exd::engine::physics::rigid_body::ElementForce3D out;
                out.r = st.r;
                out.ref = {
                    axis.origin[0] + blade.z_rotor * frame.e_z[0] + st.r * frame.e_r[0],
                    axis.origin[1] + blade.z_rotor * frame.e_z[1] + st.r * frame.e_r[1],
                    axis.origin[2] + blade.z_rotor * frame.e_z[2] + st.r * frame.e_r[2],
                };
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
    TableLookupConfig config_;
    const PolarDatabase* polars_ = nullptr;
    PolarDatabase fallback_;
};

} // anonymous namespace

std::unique_ptr<IForceEvaluator> make_table_lookup_evaluator(
    const TableLookupConfig& config, const PolarDatabase* polars)
{
    return std::make_unique<TableLookupEvaluator>(config, polars);
}

} // namespace exd::engine::physics::fluid::forces