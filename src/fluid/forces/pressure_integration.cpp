// pressure_integration.cpp
// CFD-coupled force evaluator: integrates sampled pressure (p - p_ref)·n and
// optional wall shear over the blade surface, grouped by blade element.

#include <exd/physics/fluid/forces/force_evaluator.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace exd::physics::fluid::forces
{

namespace
{

struct PerElementAccum
{
    bool used = false;
    double r = 0.0;
    std::array<double, 3> ref = {0, 0, 0};
    std::array<double, 3> force = {0, 0, 0};
    std::array<double, 3> force_pressure = {0, 0, 0};
    std::array<double, 3> force_shear = {0, 0, 0};
    std::array<double, 3> moment = {0, 0, 0};
};

class PressureIntegrationEvaluator final : public IForceEvaluator
{
public:
    explicit PressureIntegrationEvaluator(PressureIntegrationConfig config)
        : config_(config)
    {
    }

    std::string_view name() const override { return "pressure_integration"; }

    void compute(const BladeGeometry& blade,
                 const SurfaceFlow& flow,
                 double /*omega*/,
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

        const std::size_t station_count = blade.stations.size();
        if (station_count == 0) return;

        for (size_t i = 0; i < flow.element_index.size(); ++i)
        {
            const int32_t e = flow.element_index[i];
            if (e < 0 || static_cast<std::size_t>(e) >= station_count)
            {
                status.ok = false;
                status.error = "SurfaceFlow element_index out of range";
                return;
            }
        }

        // Reference point follows the p_mid convention at azimuth zero.
        const RotorFrame frame = make_rotor_frame(0.0, axis, status);
        if (!status.ok) return;

        std::vector<PerElementAccum> acc(station_count);
        for (std::size_t i = 0; i < flow.points.size(); ++i)
        {
            const std::size_t e = static_cast<std::size_t>(flow.element_index[i]);
            PerElementAccum& a = acc[e];
            if (!a.used)
            {
                a.used = true;
                a.r = blade.stations[e].r;
                a.ref = {
                    axis.origin[0] + blade.z_rotor * frame.e_z[0] + a.r * frame.e_r[0],
                    axis.origin[1] + blade.z_rotor * frame.e_z[1] + a.r * frame.e_r[1],
                    axis.origin[2] + blade.z_rotor * frame.e_z[2] + a.r * frame.e_r[2],
                };
            }

            const double p_rel = flow.pressure[i] - flow.p_ref;
            const std::array<double, 3> traction_p = {
                p_rel * flow.normals[i][0],
                p_rel * flow.normals[i][1],
                p_rel * flow.normals[i][2],
            };
            const std::array<double, 3> traction_s = config_.include_shear
                                                         ? flow.shear_traction[i]
                                                         : std::array<double, 3>{0, 0, 0};
            const double area_i = flow.area[i];

            const std::array<double, 3> dF = {
                (traction_p[0] + traction_s[0]) * area_i,
                (traction_p[1] + traction_s[1]) * area_i,
                (traction_p[2] + traction_s[2]) * area_i,
            };
            a.force[0] += dF[0];
            a.force[1] += dF[1];
            a.force[2] += dF[2];
            a.force_pressure[0] += traction_p[0] * area_i;
            a.force_pressure[1] += traction_p[1] * area_i;
            a.force_pressure[2] += traction_p[2] * area_i;
            a.force_shear[0] += traction_s[0] * area_i;
            a.force_shear[1] += traction_s[1] * area_i;
            a.force_shear[2] += traction_s[2] * area_i;

            const std::array<double, 3> arm = {
                flow.points[i][0] - a.ref[0],
                flow.points[i][1] - a.ref[1],
                flow.points[i][2] - a.ref[2],
            };
            a.moment[0] += arm[1] * dF[2] - arm[2] * dF[1];
            a.moment[1] += arm[2] * dF[0] - arm[0] * dF[2];
            a.moment[2] += arm[0] * dF[1] - arm[1] * dF[0];
        }

        // Emit in station-index ascending order, skipping unused elements.
        for (std::size_t e = 0; e < station_count; ++e)
        {
            const PerElementAccum& a = acc[e];
            if (!a.used) continue;

            mechanics::ElementForce3D out;
            out.r = a.r;
            out.ref = a.ref;
            out.force = a.force;
            out.force_pressure = a.force_pressure;
            out.force_shear = a.force_shear;
            out.moment = a.moment;
            per_element.push_back(out);
        }
    }

private:
    PressureIntegrationConfig config_;
};

} // anonymous namespace

std::unique_ptr<IForceEvaluator> make_pressure_integration_evaluator(
    const PressureIntegrationConfig& config)
{
    return std::make_unique<PressureIntegrationEvaluator>(config);
}

} // namespace exd::physics::fluid::forces