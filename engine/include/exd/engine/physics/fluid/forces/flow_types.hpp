#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::forces {

// ─────────────────────────────────────────────────────
// Flow / geometry exchange types for force evaluation.
// Domain-agnostic: any bladed machine (turbine, propeller,
// fan, pump, helicopter rotor) in any fluid.
// ─────────────────────────────────────────────────────

/// Freestream conditions in 3D.
struct Freestream
{
    std::array<double, 3> velocity = {0.0, 0.0, 0.0}; // m/s
    double rho = 1.225;                               // kg/m³
    double mu = 1.81e-5;                              // Pa·s
    double p_ref = 101325.0;                          // Pa
};

/// One spanwise blade element.
struct BladeStation
{
    double r = 0.0;          // radius from rotation axis (m)
    double dr = 0.0;         // element width (m)
    double chord = 0.0;      // local chord (m)
    double twist_deg = 0.0;  // chord angle from rotation plane toward +axis (deg)
    double thickness_ratio = 0.12; // max t/c for surface placement (0..1)
    std::string airfoil;     // airfoil id for coefficient models
};

/// Blade geometry of a rotating blade set.
struct BladeGeometry
{
    double r_hub = 0.0;                 // hub radius (m)
    double r_tip = 0.0;                 // tip radius (m)
    int blade_count = 3;                // number of blades
    double z_rotor = 0.0;               // rotor-plane axial coordinate along the
                                        // axis frame (distance from axis.origin)
    std::vector<BladeStation> stations; // ordered hub → tip
};

/// Flow sampled at surface points in 3D, grouped by blade element.
/// Produced by a coupling sampler (CFD field → surface) or by user code.
struct SurfaceFlow
{
    std::vector<std::array<double, 3>> points;          // N×3 surface point coords (m)
    std::vector<std::array<double, 3>> normals;         // N×3 outward unit normals
    std::vector<std::array<double, 3>> velocity;        // N×3 fluid velocity (m/s)
    std::vector<std::array<double, 3>> shear_traction;  // N×3 wall shear force/area (Pa);
                                                        // zeros when unused
    std::vector<double> pressure;                       // N static pressure (Pa)
    std::vector<double> area;                           // N surface area per point (m²)
    std::vector<int32_t> element_index;                 // N → station index (0..stations-1)
    double density = 1.225;                             // kg/m³
    double viscosity = 1.81e-5;                         // Pa·s
    double p_ref = 101325.0;                            // Pa reference static pressure

    /// All parallel arrays must have equal length.
    [[nodiscard]] bool valid() const
    {
        const auto n = points.size();
        return n > 0 && normals.size() == n && velocity.size() == n &&
               shear_traction.size() == n && pressure.size() == n &&
               area.size() == n && element_index.size() == n;
    }
};

} // namespace exd::engine::physics::fluid::forces
