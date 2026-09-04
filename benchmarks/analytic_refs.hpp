// ──────────────────────────────────────────────────────────────────────────
// Analytic reference math shared by the benchmark cases (harness-level; the
// analytic anchors themselves are legitimately case-local, the 2x2 linear
// machinery is not — it was duplicated between case_plenum and case_pi).
// ──────────────────────────────────────────────────────────────────────────
#pragma once

#include <algorithm>
#include <cmath>
#include <array>

namespace bench {

/// Real-eigen decomposition of a 2x2 matrix [[a11, a12], [a21, a22]] with
/// distinct real eigenvalues (asserted by the caller's parameter choice).
struct Eigen2x2 {
    double lam1 = 0.0, lam2 = 0.0;
    double v1x = 1.0, v1y = 0.0;   // eigenvector for lam1
    double v2x = 1.0, v2y = 0.0;   // eigenvector for lam2
};

inline Eigen2x2 eigen2x2(double a11, double a12, double a21, double a22) {
    const double tr = a11 + a22;
    const double det = a11 * a22 - a12 * a21;
    const double disc = tr * tr - 4.0 * det;
    Eigen2x2 e;
    e.lam1 = (tr - std::sqrt(disc)) / 2.0;
    e.lam2 = (tr + std::sqrt(disc)) / 2.0;
    e.v1x = a12; e.v1y = e.lam1 - a11;
    e.v2x = a12; e.v2y = e.lam2 - a11;
    return e;
}

/// Coefficients c1, c2 such that x(t) = c1*v1*e^(lam1 t) + c2*v2*e^(lam2 t)
/// matches x(0).  The caller guarantees a non-singular eigenbasis.
inline std::array<double, 2> expansion2x2(const Eigen2x2& e, double x0, double y0) {
    const double den = e.v1x * e.v2y - e.v2x * e.v1y;
    const double c1 = (x0 * e.v2y - y0 * e.v2x) / den;
    const double c2 = (y0 * e.v1x - x0 * e.v1y) / den;
    return {c1, c2};
}

} // namespace bench
