#pragma once

// engine_internal.hpp
// Private engine contracts: crank kinematics and cycle
// thermodynamics as PURE functions of the crank angle.
// (Tests reach this header via a relative include, the
// repo-sanctioned exception for solver internals.)

#include <exd/physics/engine/engine_config.hpp>

#include <span>

namespace exd::physics::engine {

/// All crank-derived quantities at crank angle theta.
struct CrankKinematics
{
    double x = 0.0;           // piston position from crank axis (m)
    double dx_dtheta = 0.0;   // m/rad
    double d2x_dtheta2 = 0.0; // m/rad²
    double j_eq = 0.0;        // J_flywheel + m_pist·(dx/dθ)² (kg·m²)
    double dj_dtheta = 0.0;   // 2·m_pist·(dx/dθ)·(d²x/dθ²) (kg·m²/rad)
};

CrankKinematics crank_kinematics(double theta, const EngineGeometryConfig& geometry);

/// Cylinder volume at crank angle theta (m³). Monotone with x:
/// V = V_clearance + A_pist·((l + r) − x(theta)); TDC = min.
double cylinder_volume(double theta, const EngineGeometryConfig& geometry);

/// Piston area (m²).
double piston_area(const EngineGeometryConfig& geometry);

/// Cyclic pressure (Pa) and gas temperature (K) at crank
/// angle theta; pure function of theta (no persisted cycle
/// state). Valve windows use smooth sin² ramps (no hard
/// discontinuities for RK4/adaptive integrators).
void cylinder_state(double theta,
                    const EngineGeometryConfig& geometry,
                    const EngineThermoConfig& thermo,
                    double& p_cyl,
                    double& T_cyl);

/// Load moment opposing rotation at omega (N·m):
/// friction_constant + friction_viscous·ω + generator curve.
double load_moment(const EngineLoadConfig& load, double omega);


/// The engine ODE right-hand side for state [theta, omega]:
///
///   dtheta = omega
///   domega = (M_gas(theta) − M_load(omega) − ½·(dJ/dθ)·omega²) / J_eq(theta)
///
/// `throttle` scales heat release (Otto; ignored for steam). This is the
/// canonical RHS — step_engine integrates it, and tests verify energy
/// bookkeeping against it with sub-step quadrature.
void engine_rhs(const EngineConfig& config, double throttle,
                std::span<const double> state, std::span<double> dstate);

} // namespace exd::physics::engine
