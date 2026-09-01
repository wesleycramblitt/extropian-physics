# Capability Matrix — Real Multiphysics Runs

What each system can do TODAY, what is genuinely simulated (vs engineering
model), and the exact configuration to run it. The library is a static C++
library; a "run" is a caller program linking `exd-physics` and calling the
entry points below (pattern + recipes: `docs/real_run_guide.md`).

## 1. Wind turbine / water turbine — real multiphysics over a mesh

| Piece | What runs | Status |
|---|---|---|
| Fluid | 3D incompressible FDM over a **structured mesh** (staggered-role collocated SIMPLE): velocity/pressure fields, 6-face BCs (Inlet/Outlet/Wall/Symmetry/Periodic/FixedPressure), SOR pressure solve, adaptive dt | real solver (verified: uniform flow exact, Poiseuille 0.13%, Taylor–Green 4.7%) |
| Rotor coupling | blade-element forces from the **local CFD field** (no induction double-count), distributed as smeared body force (actuator-disk), under-relaxed + ramped | real coupling (soak: wake deficit ~20%, cp ratio 0.63 vs reduced-order, energy balance closed) |
| Rotor dynamics | 1-DOF spin-up against generator load curve | real |
| Generator/electric | DC motor / curve loads (moment models) | real |
| Control | PI speed governor (load-fraction) regulates ω | real (setpoint tracking < 0.1%) |
| Output | velocity+pressure field stamps (`.fld` + timeline), rotor state CSV | real (streamed per step) |

**Wind** — any `V_inf`, `rho = 1.225`, `mu = 1.81e-5` via
`turbine::default_grid_config` + `run_coupled_turbine`.

**Water** — same code path, working fluid only: `rho = 1025`, `mu = 1.08e-3`
(seawater); scale the generator torque curve and rotor inertia by
`rho_w/rho_a ≈ 837` (aero torque scales with ρ). Verified by the hydro test
(same settle, Betz-bounded cp, finite states).

Boundaries of the model (documented): structured grid only (unstructured
meshes + FVM/FEM are Phase J); collocated scheme is engineering-grade;
no blade structural response yet (per-element forces ARE exported for a
future structural solver); no free-surface/cavitation.

## 2. Combustion engine — real multiphysics (0D/1D cycle)

| Piece | What runs | Status |
|---|---|---|
| Mechanism | analytic slider-crank: piston x(θ), v(θ), J_eq(θ) with the ½(dJ/dθ)ω² inertia torque | real (TDC/BDC, energy-conservation tests) |
| Cylinder | 4-stroke cycle: polytropic compression/expansion, Wiebe heat release anchored at TDC volume, sin² valve windows, per-phase γ | real thermodynamics (η = 0.384 vs 0.489 Otto bound) |
| Control | PI governor throttling heat release → speed regulation | real |
| Loads | friction + generator curve (pure T(ω)); DC motor via moment models | real |
| Output | per-step θ, ω, piston x/v, p_cyl, T_cyl, torque, power, throttle → CSV | real |

The 0D/1D cylinder is the industry-standard fidelity for cycle studies
(GT-Power class); a moving-piston 3D in-cylinder CFD is a Phase J-class
project, not required for the capability.

## 3. Steam engine — real multiphysics (Rankine-lite) ✅ upgraded

| Piece | What runs | Status |
|---|---|---|
| Steam EOS | `thermo::steam`: Clausius–Clapeyron saturation line (p_sat(T)/T_sat(p)), liquid/vapor/latent enthalpies, saturated vapor density — engineering model, ~5–10% vs IAPWS across 0–200 °C | real |
| Cycle | single-acting: saturated admission to cutoff, wet-steam polytrope expansion (x = x_cut·(V_cut/V)^(n−1), n=1.13), exhaust to condenser; cylinder temperature on the saturation line while wet | real |
| Energy | per-cycle boiler heat m·(h_g − h_f) (Rankine-lite bookkeeping); `efficiency_estimate = W/heat` | real |
| Output | same per-step CSV as the combustion engine (θ, ω, piston, p/T, torque, power) | real |

Validation: admission temperature = T_sat(p_boiler) within 0.5%; efficiency
in a plausible simple-engine band; steam properties unit-tested (anchor,
round-trips, monotonicity, enthalpy consistency).

Boundary of the model: single-constant latent heat (no superheat schedule,
no wetness-at-exhaust reporting); swap `SteamConstants` for IF97 tables
without touching the cycle.

## 3b. Compressor / turbomachinery — mean-line machine capability (W9) ✅ new

| Piece | What runs | Status |
|---|---|---|
| Compressor row | axial mean-line stage/stack: velocity-triangle Euler work, polytropic exit state (total-state closure), relative-Mach choke guard, documented envelope (M_rel < 0.7, hub/tip > 0.5) | real solver (geometry-emergent sense: the SAME code produces compressor AND turbine — reversal + turbocharger tests) |
| Operating map | `solve_operating_map` (ω × ṁ sweep → π/η/T surfaces, surge + choke lines) / `sample_operating_map` (bilinear; accepts rig-test maps) | real (surge line = dπ/dṁ = 0 surrogate, documented) |
| Surge dynamics | `fluid::lumped::plenum` — Greitzer 2-state cell on `IEos`, Jacobian-verified stability, surge limit cycle | real (verified both sides of the stability boundary) |
| Drive + control | Phase D DC motor + Phase C PI governor (throttle-gain modulation) | real (spin-up settles at operating point; governing 0.0002%) |
| Output | per-step CSV (t, ω, p, ṁ, π, torques, gain) | real (one row per step) |

Run: caller program linking `exd-physics`, calling
`fluid::turbomachinery::simulate_compression_system` (or `solve_stage_stack`
for pure/batchable stack solves — optimizer-ready). Turbocharger
(compressor+turbine on one shaft) is a verified acceptance test.

## 3c. Grid-first domains + real coupling (W10) ✅ new

| System | What runs | Verification |
|---|---|---|
| Thermal | steady conduction + advection (SOR), fixed/insulated faces, uniform source | linear profile exact, parabolic profile 1%, Péclet warning |
| Acoustics | scalar wave equation (leapfrog), pressure-release box | plane-wave period 1%, box eigenfrequencies 0.1–1%, energy bounded 10 periods |
| Structural | linear elasticity displacement form (SOR), thermal-strain channel, masked traction | uniaxial (ν=0) + oedometer + thermal column EXACT; bending qualitative + load-linearity |
| Particles | Lagrangian cloud over sampled flow channels (RK4) | ballistic exact, terminal velocity, channel advection |
| Chemistry | 0D mass-action + Arrhenius reactor | decay/equilibrium/second-order analytic, Arrhenius ratio 0.1% |
| Coupling | real `CouplingManager` (nearest/trilinear transfer, relaxed implicit links), `CoupledSimulation` multi-rate driver | staggered vs implicit both land on analytic; multi-rate converges |

Run: caller programs link `exd-physics`; coupling is configured in code
(domains register channels/sinks, links carry probe points) — pattern in
`docs/coupling_and_grid_first_domains.md`.

## 4. Shared capability guarantees

- **No exceptions**: `ModelStatus` error channels everywhere; invalid
  configs fail with messages, never crash.
- **Deterministic + batchable**: `simulate_engine`, `solve_fdm3`,
  `run_coupled_turbine` are pure functions of their config; identical
  inputs → identical outputs (pinned by tests) → optimizer loops run them
  with no I/O overhead.
- **Real-time output**: `OutputPolicy`/`OutputScheduler` with injected
  clock; binary field stamps + CSV streams (contract:
  `docs/output_channels.md`).
- Verified by 74/74 unit tests including analytic anchors for every
  physics model.
