// blade_surface.cpp
// Rotor frame construction and flat-plate lifting-surface sampling for the
// exd::engine::physics::fluid::forces domain. Builds the mid-chord upper/lower sample points of a
// bladed rotor at one azimuth and the local orthonormal (e_r, e_t, e_z)
// frame used by all 3D force evaluators.

#include <exd/engine/physics/fluid/forces/force_evaluator.hpp>

#include <array>
#include <cmath>
#include <cstdint>

namespace exd::engine::physics::fluid::forces
{

namespace
{

constexpr double kNearZero = 1e-12;

/// Normalize a 3-vector in place. Returns false when |v| is (near) zero.
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

} // namespace

RotorFrame make_rotor_frame(double angle_rad, const exd::engine::physics::rigid_body::RotationAxis& axis,
                            exd::engine::physics::rigid_body::ModelStatus& status)
{
    RotorFrame frame;

    // Axis direction (normalized copy).
    std::array<double, 3> e_z = axis.direction;
    if (!normalize3(e_z))
    {
        status.ok = false;
        status.error = "axis direction is (near) zero";
        return frame; // default frame
    }

    // Reference vector used to seed the rotation plane; switch away from the
    // near-axis case where u is (almost) parallel to e_z.
    std::array<double, 3> u = {1.0, 0.0, 0.0};
    const double u_dot_ez = u[0] * e_z[0] + u[1] * e_z[1] + u[2] * e_z[2];
    if (std::fabs(u_dot_ez) > 0.9)
        u = {0.0, 1.0, 0.0};

    std::array<double, 3> e_r0 = {
        u[0] - u_dot_ez * e_z[0],
        u[1] - u_dot_ez * e_z[1],
        u[2] - u_dot_ez * e_z[2],
    };
    if (!normalize3(e_r0))
    {
        status.ok = false;
        status.error = "axis direction is (near) zero";
        return frame;
    }

    // e_t0 = e_z × e_r0.
    std::array<double, 3> e_t0 = {
        e_z[1] * e_r0[2] - e_z[2] * e_r0[1],
        e_z[2] * e_r0[0] - e_z[0] * e_r0[2],
        e_z[0] * e_r0[1] - e_z[1] * e_r0[0],
    };

    const double c = std::cos(angle_rad);
    const double s = std::sin(angle_rad);

    // e_r(theta) = e_r0·cos + e_t0·sin; e_t(theta) = -e_r0·sin + e_t0·cos.
    frame.e_z = e_z;
    frame.e_r = {
        e_r0[0] * c + e_t0[0] * s,
        e_r0[1] * c + e_t0[1] * s,
        e_r0[2] * c + e_t0[2] * s,
    };
    frame.e_t = {
        -e_r0[0] * s + e_t0[0] * c,
        -e_r0[1] * s + e_t0[1] * c,
        -e_r0[2] * s + e_t0[2] * c,
    };
    return frame;
}

BladeSurface build_blade_surface(const BladeGeometry& blade, double angle_rad,
                                 const exd::engine::physics::rigid_body::RotationAxis& axis)
{
    BladeSurface surface;

    const int station_count = static_cast<int>(blade.stations.size());
    const int blade_count = blade.blade_count;
    if (station_count <= 0 || blade_count <= 0) return surface;

    const std::size_t total = static_cast<std::size_t>(2) *
                              static_cast<std::size_t>(station_count) *
                              static_cast<std::size_t>(blade_count);
    surface.points.reserve(total);
    surface.normals.reserve(total);
    surface.areas.reserve(total);
    surface.element_index.reserve(total);

    exd::engine::physics::rigid_body::ModelStatus frame_status;
    for (int b = 0; b < blade_count; ++b)
    {
        const double theta = angle_rad +
                             static_cast<double>(b) * 2.0 * M_PI / static_cast<double>(blade_count);
        const RotorFrame frame = make_rotor_frame(theta, axis, frame_status);
        if (!frame_status.ok) return BladeSurface{}; // degenerate axis

        for (int i = 0; i < station_count; ++i)
        {
            const BladeStation& st = blade.stations[static_cast<std::size_t>(i)];
            const double beta = st.twist_deg * M_PI / 180.0;
            const double cb = std::cos(beta);
            const double sb = std::sin(beta);

            // Chord direction and face normal in the e_t-e_z plane:
            //   c_hat =  cos(beta)·e_t + sin(beta)·e_z
            //   n_hat = -sin(beta)·e_t + cos(beta)·e_z  (points +axial at beta = 0)
            const std::array<double, 3> n_hat = {
                -sb * frame.e_t[0] + cb * frame.e_z[0],
                -sb * frame.e_t[1] + cb * frame.e_z[1],
                -sb * frame.e_t[2] + cb * frame.e_z[2],
            };

            const std::array<double, 3> p_mid = {
                axis.origin[0] + blade.z_rotor * frame.e_z[0] + st.r * frame.e_r[0],
                axis.origin[1] + blade.z_rotor * frame.e_z[1] + st.r * frame.e_r[1],
                axis.origin[2] + blade.z_rotor * frame.e_z[2] + st.r * frame.e_r[2],
            };

            const double half_thickness = 0.5 * st.thickness_ratio * st.chord;
            const double area = 0.5 * st.dr * st.chord;

            const std::array<double, 3> p_upper = {
                p_mid[0] + half_thickness * n_hat[0],
                p_mid[1] + half_thickness * n_hat[1],
                p_mid[2] + half_thickness * n_hat[2],
            };
            const std::array<double, 3> p_lower = {
                p_mid[0] - half_thickness * n_hat[0],
                p_mid[1] - half_thickness * n_hat[1],
                p_mid[2] - half_thickness * n_hat[2],
            };

            // Upper side (normal = n_hat).
            surface.points.push_back(p_upper);
            surface.normals.push_back(n_hat);
            surface.areas.push_back(area);
            surface.element_index.push_back(static_cast<int32_t>(i));

            // Lower side (normal = -n_hat).
            surface.points.push_back(p_lower);
            surface.normals.push_back({-n_hat[0], -n_hat[1], -n_hat[2]});
            surface.areas.push_back(area);
            surface.element_index.push_back(static_cast<int32_t>(i));
        }
    }
    return surface;
}

} // namespace exd::engine::physics::fluid::forces