// 3D rotational moment integration about a fixed axis.
// Pure geometry: reduces a set of element forces/moments to the
// scalar torque/axial-force pair plus full 3D resultants.

#include <exd/physics/mechanics/rotational_state.hpp>

#include <array>
#include <cmath>
#include <span>

namespace exd::physics::mechanics {

namespace
{

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> cross(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

} // anonymous namespace

bool normalize(std::array<double, 3>& v)
{
    const double length = std::sqrt(dot(v, v));
    if (length <= 1e-12)
        return false;

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
    return true;
}

MomentResult integrate_moment(std::span<const ElementForce3D> forces,
                              const RotationAxis& axis)
{
    MomentResult result;

    // Normalize a local copy of the axis direction; a degenerate direction
    // makes the result invalid (scalars remain all-zero).
    std::array<double, 3> direction = axis.direction;
    if (!normalize(direction))
        return result;

    for (const ElementForce3D& f : forces)
    {
        const std::array<double, 3> lever = {
            f.ref[0] - axis.origin[0],
            f.ref[1] - axis.origin[1],
            f.ref[2] - axis.origin[2],
        };
        const std::array<double, 3> r_cross_f = cross(lever, f.force);

        result.torque += dot(r_cross_f, direction);
        result.axial_force += dot(f.force, direction);

        result.total_force[0] += f.force[0];
        result.total_force[1] += f.force[1];
        result.total_force[2] += f.force[2];

        // Full resultant moment about the axis origin includes the element's
        // own moment (rigid-body transfer); element moments do not enter the
        // scalar torque evaluated along the axis.
        result.total_moment[0] += r_cross_f[0] + f.moment[0];
        result.total_moment[1] += r_cross_f[1] + f.moment[1];
        result.total_moment[2] += r_cross_f[2] + f.moment[2];
    }

    result.valid = true;
    return result;
}

} // namespace exd::physics::mechanics