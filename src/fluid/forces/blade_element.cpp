// blade_element.cpp
// Local blade-element force evaluation without an induction solve.
//
// This variant is for CFD-coupled (actuator-disk / actuator-line) runs: for
// each station it uses the MEAN of the sampled flow velocities at that
// station directly (all surface points with element_index == station) to get
// the local axial velocity, then assembles the per-blade lab-frame forces
// exactly like momentum_balance.cpp / table_lookup.cpp.
//
// IMPORTANT: do NOT use this variant with an induction-solving sampler
// upstream of the disk — the sampled flow already contains the induction
// (the wake deficit), so an additional BEM induction solve would DOUBLE
// COUNT the induction.  Use MomentumBalance for standalone reduced-order
// runs, and this variant for CFD-coupled runs.

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

constexpr double kNearZero = 1e-6;      // stall guard on the local speed W
constexpr char kFallbackAirfoil[] = "naca0012";

double normalize3(std::array<double, 3>& v)
{
    const double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len < 1e-12) return 0.0;
    return 1.0 / len;
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

class BladeElementEvaluator final : public IForceEvaluator
{
public:
    BladeElementEvaluator(BladeElementConfig /*config*/, const PolarDatabase* polars)
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

    std::string_view name() const override { return "blade_element"; }

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

        std::array<double, 3> e_z = axis.direction;
        if (e_z[0] == 0.0 && e_z[1] == 0.0 && e_z[2] == 0.0)
        {
            status.ok = false;
            status.error = "axis direction is (near) zero";
            return;
        }
        normalize3(e_z);

        const int blade_count = blade.blade_count;
        bool warned_fallback = false;
        bool all_stalled = true;

        for (std::size_t st_i = 0; st_i < blade.stations.size(); ++st_i)
        {
            const BladeStation& st = blade.stations[st_i];
            const double r = std::max(st.r, 1e-9);

            // Mean sampled velocity at this station (all points that carry
            // this station's element_index — an actuator-disk ring average).
            std::array<double, 3> v_mean = {0.0, 0.0, 0.0};
            std::size_t count = 0;
            for (std::size_t i = 0; i < flow.points.size(); ++i)
            {
                if (static_cast<std::size_t>(flow.element_index[i]) == st_i)
                {
                    v_mean[0] += flow.velocity[i][0];
                    v_mean[1] += flow.velocity[i][1];
                    v_mean[2] += flow.velocity[i][2];
                    ++count;
                }
            }
            if (count > 0)
            {
                const double inv = 1.0 / static_cast<double>(count);
                v_mean[0] *= inv;
                v_mean[1] *= inv;
                v_mean[2] *= inv;
            }
            // No samples (or all out-of-bounds, giving zero velocity) => the
            // station sees no local flow; the stall guard below zeroes dL/dD.
            if (count == 0)
            {
                v_mean = {0.0, 0.0, 0.0};
            }

            // Axial speed: positive magnitude for the convention inflow in the
            // -e_z direction (like momentum_balance, thrust into the inflow).
            const double v_ax = v_mean[0] * e_z[0] +
                                v_mean[1] * e_z[1] +
                                v_mean[2] * e_z[2];
            const double u_a = -v_ax;

            // Relative velocity: axial component plus the blade speed.
            const double vt = omega * r;
            const double W = std::sqrt(u_a * u_a + vt * vt);

            double dQ = 0.0;
            double dT = 0.0;

            if (W < kNearZero)
            {
                // Local speed at (near) zero: the station produces no blade
                // force (one collective warning reported by the caller when
                // ALL stations stall).
                // (all_stalled stays true)
            }
            else
            {
                all_stalled = false;
                const double phi = std::atan2(u_a, vt);
                const double sin_phi = std::sin(phi);
                const double cos_phi = std::cos(phi);
                const double alpha_deg = phi * 180.0 / M_PI - st.twist_deg;

                const double re = (flow.viscosity > 0.0)
                    ? flow.density * W * st.chord / flow.viscosity
                    : 0.0;
                const auto [cl, cd] = resolve_coefficients(*polars_, st.airfoil,
                                                           alpha_deg, re,
                                                           warned_fallback, status);

                // Ring totals (B blades) per element, mirroring the element
                // force loop of momentum_balance.cpp exactly (no induction).
                const double dL = 0.5 * flow.density * W * W * st.chord * cl * st.dr;
                const double dD = 0.5 * flow.density * W * W * st.chord * cd * st.dr;
                dQ = static_cast<double>(blade_count) * r *
                     (dL * sin_phi - dD * cos_phi);
                dT = static_cast<double>(blade_count) *
                     (dL * cos_phi + dD * sin_phi);
            }

            // Per-blade 3D loads distributed over the rotating blades (same
            // lab-frame assembly as momentum_balance.cpp).
            const double per_blade_axial = dT / static_cast<double>(blade_count);
            const double per_blade_tangential =
                (r > 1e-9) ? dQ / (r * static_cast<double>(blade_count)) : 0.0;
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

        if (all_stalled)
        {
            status.warnings.push_back("blade_element: no local flow at any "
                                     "station; blade forces are zero");
        }
    }

private:
    const PolarDatabase* polars_ = nullptr;
    PolarDatabase fallback_;
};

} // anonymous namespace

std::unique_ptr<IForceEvaluator> make_blade_element_evaluator(
    const BladeElementConfig& config, const PolarDatabase* polars)
{
    return std::make_unique<BladeElementEvaluator>(config, polars);
}

} // namespace exd::physics::fluid::forces