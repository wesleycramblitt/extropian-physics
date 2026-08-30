// crank_mechanism.cpp
// Analytic slider-crank kinematics: piston position/velocity,
// equivalent inertia J_eq(θ) and its derivative.
//
//   x(θ)  = r·cosθ + sqrt(l² − r²·sin²θ)          (from crank axis)
//   dx/dθ = −r·sinθ − (r²·sinθ·cosθ)/s
//   d²x/dθ² = −r·cosθ − (r²·cos2θ)/s − (r²·sinθ·cosθ)²/s³
//   J_eq(θ) = J_f + m_pist·(dx/dθ)²
//   dJ/dθ  = 2·m_pist·(dx/dθ)·(d²x/dθ²)

#include "engine_internal.hpp"

#include <cmath>

namespace exd::physics::engine {

namespace
{
constexpr double PI = 3.14159265358979323846;
}

CrankKinematics crank_kinematics(double theta, const EngineGeometryConfig& g)
{
    const double r = g.crank_radius;
    const double r2 = r * r;
    const double sin_t = std::sin(theta);
    const double cos_t = std::cos(theta);
    const double s = std::sqrt(g.rod_length * g.rod_length - r2 * sin_t * sin_t);

    CrankKinematics k;
    k.x = r * cos_t + s;
    k.dx_dtheta = -r * sin_t - (r2 * sin_t * cos_t) / s;
    const double u = r2 * sin_t * cos_t;
    k.d2x_dtheta2 = -r * cos_t - (r2 * std::cos(2.0 * theta)) / s - (u * u) / (s * s * s);

    const double dx = k.dx_dtheta;
    k.j_eq = g.flywheel_inertia + g.piston_mass * dx * dx;
    k.dj_dtheta = 2.0 * g.piston_mass * dx * k.d2x_dtheta2;
    return k;
}

double cylinder_volume(double theta, const EngineGeometryConfig& g)
{
    const double x = crank_kinematics(theta, g).x;
    // TDC: x = l + r → V = V_clearance. Piston moves DOWN as x decreases.
    return g.clearance_volume + piston_area(g) * ((g.rod_length + g.crank_radius) - x);
}

double piston_area(const EngineGeometryConfig& g)
{
    return 0.25 * 3.14159265358979323846 * g.bore * g.bore;
}

} // namespace exd::physics::engine
