#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace exd::physics::fluid::reduced_order::bem
{

/// Freestream / operating conditions for a BEM solve.
struct OperatingConditions
{
    double v_inf = 10.0;                // m/s, > 0 (validation)
    double rho   = 1.225;               // kg/m^3, air default
    double mu    = 1.81e-5;             // Pa*s  (Re = rho*W*c/mu)
    double p_ref = 101325.0;            // Pa, freestream static pressure
    std::optional<double> rpm_override; // else BladeRow::rotational_speed
};

/// Denominator used for non-dimensional Cp/Ct.
enum class ReferenceArea { RotorDisk, Annulus };

/// Induction correction model for the high-induction regime (a > threshold).
/// - Standard:       Buhl closed-form algebraic solution to Glauert's empirical
///                   C_T = 4a(1 - 0.25(5-3a)a). One-shot, no iteration needed.
/// - GlauertIterative: Iterative application of Glauert's empirical correction.
///                   Converges to the same answer as Buhl for steady-state but
///                   allows under-relaxation between momentum and empirical
///                   regimes. Canonical form found in NREL AeroDyn.
/// - Snel:           Exponential blending between momentum and empirical.
///                   a_corrected = a_momentum * (1 - exp(-4F sin^2(phi) /
///                   (sigma Cn) * a_momentum)). Smooth transition, no hard
///                   threshold. Snel et al. (1993).
enum class InductionCorrection { Standard, GlauertIterative, Snel };

/// Tip/hub loss correction model.
/// - Prandtl:         Classic Prandtl tip/hub loss factor
///                    F = (2/pi) * acos(exp(-f)), applied multiplicatively.
/// - DuSelig:         Du-Selig (1993) modified tip loss with loading-dependent
///                    correction. Better agreement with experiment near the tip.
/// - Chaviaropoulos:  Chaviaropoulos-Hansen (2000) tip loss with
///                    loading-dependent modification.
enum class LossCorrection { Prandtl, DuSelig, Chaviaropoulos };

/// Span (0 = hub, 1 = shroud) -> airfoil id mapping.
struct AirfoilAssignment
{
    double span = 0.0;
    std::string airfoil;
};

/// Configuration for the Level-3 duct-coupled BEM solver.
struct BEMSolverConfig
{
    uint32_t element_count = 32;        // 0 disallowed; >= 4 validated; warning if
                                        //   not in {16,32,64,128} (spec's sweep set)
    double   k_duct = 0.5;              // [0,1] duct acceleration coefficient
    double   under_relaxation = 0.25;   // (0,1]
    double   induction_tolerance = 1e-5;
    uint32_t max_iterations = 100;
    double   glauert_threshold = 0.4;   // Buhl branch guard (see induction.hpp)
    uint32_t row_index = 0;
    double   hull_cd = 0.2;             // total hull drag C_D (frontal ref area)
    ReferenceArea reference_area = ReferenceArea::RotorDisk;  // A=pi*R_tip^2 per spec;
                                    // Annulus: pi*(R_tip^2 - R_hub^2) — both documented
    uint32_t field_axial_points = 64;
    uint32_t field_radial_points = 16;
    double   upstream_extent = 0.0;     // 0 -> 3*R_tip
    double   wake_length     = 0.0;     // 0 -> 8*R_tip
    double   wake_decay_length = 0.0;   // 0 -> 4*R_tip
    double   wake_expansion  = 0.05;    // k_w, dR_w/dz
    double   wake_radius_initial = 0.0; // 0 -> R_tip
    std::vector<AirfoilAssignment> airfoils;  // empty -> all "naca0012"
    bool     include_flow_field = true; // cost is negligible (~us at 64x16)

    // Correction model selection.
    InductionCorrection induction_correction = InductionCorrection::Standard;
    LossCorrection      loss_correction      = LossCorrection::Prandtl;
};

} // namespace exd::physics::fluid::reduced_order::bem
