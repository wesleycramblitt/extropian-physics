#pragma once

// ---------------------------------------------------------------------------
// Partitioned FSI-lite: a 6-DOF rigid body moving under gravity in an
// incompressible fdm3 flow, two-way coupled through a Gaussian-smeared
// point force (W12).
//
// Model (v1):
//   body   m·dv/dt = m·g + F_drag,   F_drag = -½·ρ·(C_d·A)·|v-v_f|·(v-v_f)
//   fluid  momentum source per cell: a_i = w_i·(-F_drag)/(ρ·V_cell)
//          with Gaussian weights w_i = exp(-r²/2ε²) DISCRETE-normalized so
//          that Σw_i = 1 over the ±4ε support (apply_smeared_point_force;
//          the exact sum is unit-pinned).  In a CLOSED box the naive
//          partitioned identity m_b·v_b + ∫ρu dV = m_b·g·t carries a wall-
//          pressure correction: once the entrainment column reaches a wall
//          the pressure returns momentum to the ledger (measured deficit ≈
//          blob-mass × probe velocity).  The integrated two-way anchor is
//          therefore the open-box (wind-tunnel) absolute terminal speed
//          v_t = sqrt(2·m·g/(ρ·C_d·A)) -- the smeared reaction momentum
//          leaves the domain and the drag probe stays clean.
//   loop   staggered: per exchange -- sample v_f AHEAD of the body (an
//          upstream probe, see below), compute F_drag on RELATIVE velocity,
//          apply the negated force to the fluid (β-under-relaxed), advance
//          the body one window (loads held constant -- the documented Phase
//          B contract), then run fluid_steps_per_exchange fluid steps under
//          the smeared source.
//
// v1 deliberately treats the body as a POINT with an UPSTREAM SENSOR: the
// fluid velocity entering the drag law is sampled `sample_lead·ε` ahead of
// the CoM along the oncoming relative flow (i.e. the attacker side of the
// body).  Sampling at the CoM itself would measure the fluid the body's own
// smeared reaction just pushed -- the classical point-force self-induction,
// which rings up the two-way loop (measured: drag sign flips once the local
// v_f catches the body speed).  The upstream sensor excludes that
// self-induced part, so the drag acts against the flow the body is about to
// enter; wake feedback AT the body surface is W13 territory (true
// body-fitted sampling).  No torque, no orientation coupling in v1.
// Body-frame inertia is carried for interface symmetry with mechanics.
//
// Determinism: fixed loops, no atomics -- identical configs give
// bit-identical results.  No exceptions; failures surface through ModelStatus.
// ---------------------------------------------------------------------------

#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/physics/mechanics/rigid_body.hpp>
#include <exd/physics/model_status.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::coupling {

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

struct DragBodyConfig
{
    // ── Body ────────────────────────────────────────────────
    double mass = 0.0;                                 // kg, > 0
    std::array<double, 3> inertia_principal = {1.0, 1.0, 1.0}; // kg·m², each > 0
    double drag_area = 0.0;                            // C_d·A (m²), >= 0
    std::array<double, 3> initial_position = {0.5, 0.5, 0.5}; // m (CoM)
    std::array<double, 3> initial_velocity = {0.0, 0.0, 0.0}; // m/s
    std::array<double, 3> gravity = {0.0, 0.0, -9.81}; // m/s²

    // ── Coupling ───────────────────────────────────────────
    int fluid_steps_per_exchange = 4;   // fluid steps per body step, >= 1
    double force_relaxation = 0.5;      // β in (0,1]; the fluid force applied
                                        // is β·F_new + (1-β)·F_previous.
                                        // Use 1.0 when exact impulse
                                        // bookkeeping matters (acceptance).
    double smear_cells = 1.5;           // ε = smear_cells · min(dx,dy,dz), >= 0.5
    double sample_lead = 5.0;           // upstream probe distance in units of
                                        // ε, >= 0 (0 = sample at the CoM --
                                        // then the point-force self-induction
                                        // applies and the loop rings up)
    uint64_t flow_precondition_steps = 0; // run the flow to a steady state
                                        // WITHOUT the body force before the
                                        // first exchange (wind-tunnel setups:
                                        // the duct must be established or its
                                        // fill transients contaminate the
                                        // drag probe), >= 0

    // ── Flow (fdm3). MUST use fixed dt (adaptive_dt == false): the body
    //    window is fluid_steps_per_exchange · flow.dt. ─────────
    fluid::fdm3::FDM3Config flow;

    uint64_t max_steps = 2000;          // body steps (exchanges), > 0
};

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------

struct DragBodyResult
{
    bool ok = false;
    ModelStatus status;

    std::vector<double> probe_time;             // per exchange, after the body step (s)
    std::vector<std::array<double, 3>> probe_position;
    std::vector<std::array<double, 3>> probe_velocity;
    std::vector<std::array<double, 3>> probe_fluid_momentum; // Σ ρ·u·dV per exchange

    double terminal_velocity = 0.0;             // |v| after the final exchange (m/s)
    double drag_at_end = 0.0;                   // |F_drag| at the final exchange (N)
    std::array<double, 3> sampled_fluid_velocity_final = {0.0, 0.0, 0.0}; // v_f at CoM (m/s)

    std::array<double, 3> fluid_momentum = {0.0, 0.0, 0.0}; // Σ ρ·u·dV at the end (kg·m/s)
    uint64_t exchanges = 0;
};

// ---------------------------------------------------------------------------
// Drag law + smearing (exported for unit tests)
// ---------------------------------------------------------------------------

/// Quadratic drag on the RELATIVE velocity:
///   F = -½·ρ·(C_d·A)·|v_rel|·v_rel
/// Zero at rest; aligned against v_rel.
std::array<double, 3> body_drag_force(const std::array<double, 3>& relative_velocity,
                                      double rho, double drag_area);

/// Build the per-cell momentum source (acceleration, m/s²) for a point force
/// `force` applied at `position`: Gaussian weights over the ±4ε support,
/// DISCRETE-normalized so that Σcells ρ·a·dV ≡ force (exactly, up to the
/// arithmetic of the sum).  Cells outside the grid are clipped, which breaks
/// the identity -- the caller must keep the blob inside the domain.  The
/// arrays are sized nx·ny·nz and reset to zero; `force` is the action on
/// the FLUID (the negated body drag).
void apply_smeared_point_force(std::vector<double>& fx, std::vector<double>& fy,
                               std::vector<double>& fz,
                               const fluid::fdm3::FDM3Config& flow,
                               const std::array<double, 3>& position,
                               const std::array<double, 3>& force, double eps,
                               double rho);

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Run the partitioned FSI-lite simulation.  Validation failures return
/// `result.ok == false` with `status.error` describing the problem; the
/// body leaving the box or approaching the walls degrades conservation and
/// is reported as a warning.  Deterministic.
DragBodyResult simulate_drag_body(const DragBodyConfig& config, ModelStatus& status);

} // namespace exd::physics::coupling
