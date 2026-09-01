#pragma once

// ─────────────────────────────────────────────────────
// Polytropic process primitives (STAGNATION family).
//
// Relations between total-temperature ratio
//   τ = T02/T01  and  total-pressure ratio π = p02/p01
// through a polytropic efficiency η_p (0 < η_p ≤ 1):
//
//   compression:  τ = π^((γ−1)/(γ·η_p))          π = τ^(γ·η_p/(γ−1))
//   expansion:    τ = π^(η_p·(γ−1)/γ)            π = τ^(γ/((γ−1)·η_p))
//
// η_p = 1 reproduces the exact isentropic relations. Compressor and
// turbine polytropes are NOT reciprocal (π_comp·π_turb ≠ 1 at η_p < 1).
//
// Deliberately NOT wired into IEos: an equation of state is a state law
// (p(ρ,T), h, s...); a polytrope is a PROCESS law (relates two states
// through an efficiency). The engine's in-cylinder p·Vⁿ volume polytropes
// are heat-transfer stand-ins and stay local to `engine` (Phase G §10 G.2).
// ─────────────────────────────────────────────────────

#include <exd/physics/model_status.hpp>

namespace exd::physics::thermo::polytropic {

/// Total-temperature ratio τ = T02/T01 across a compressor stage/process
/// with polytropic efficiency `eta_poly` and pressure ratio `pi` (> 0).
double temp_ratio_compression(double pi, double gamma, double eta_poly,
                              exd::physics::ModelStatus& status);

/// Total-temperature ratio τ = T02/T01 across an expansion (turbine)
/// process with polytropic efficiency `eta_poly` and pressure ratio
/// `pi` (> 0, < 1 for a real expansion).
double temp_ratio_expansion(double pi, double gamma, double eta_poly,
                            exd::physics::ModelStatus& status);

/// Total-pressure ratio π = p02/p01 from a temperature ratio `tau`.
/// Inverse of temp_ratio_compression.
double pressure_ratio_compression(double tau, double gamma, double eta_poly,
                                  exd::physics::ModelStatus& status);

/// Total-pressure ratio π = p02/p01 from a temperature ratio `tau`.
/// Inverse of temp_ratio_expansion.
double pressure_ratio_expansion(double tau, double gamma, double eta_poly,
                                exd::physics::ModelStatus& status);

/// Recover the polytropic efficiency of a compression process from
/// measured/derived (π, τ). Domain: 0 < π, 0 < τ; result clamped to (0,1].
double polytropic_efficiency_compression(double pi, double tau, double gamma,
                                         exd::physics::ModelStatus& status);

/// Recover the polytropic efficiency of an expansion process from
/// measured/derived (π, τ). Domain: 0 < π, 0 < τ; result clamped to (0,1].
double polytropic_efficiency_expansion(double pi, double tau, double gamma,
                                       exd::physics::ModelStatus& status);

} // namespace exd::physics::thermo::polytropic
