# Capability Matrix — Real Multiphysics Runs

What each system can do TODAY, what is genuinely simulated (vs engineering
model), and the exact configuration to run it. The library is a static C++
library; a "run" is a caller program linking `exd-physics` and calling the
entry points below (pattern + recipes: `docs/real_run_guide.md`).

## 0. W14 FDM domain breadth + multiphysics presets

| System | Run entry | Fidelity/anchors |
|---|---|---|
| Species transport | `physics::species::solve_species` | operator-split advection/diffusion/reaction on the core runtime; exact decay e^{−kt}, A→B conservation, advective-decay steady profile, diffusion variance 2Dt |
| Elastic waves | `physics::structural::solve_elasticity(config.transient)` | velocity-Verlet + mirror ghosts; P-wave arrival within 2%, flight-energy drift < 5% |
| Porous media | `physics::porous::solve_porous` | Darcy pressure diffusion, implicit CG; direct steady solve via the affine-operator linear part; exact linear steady profile |
| Aeroacoustics | `presets::multiphysics::run_aeroacoustics` | fdm3 duct → acoustics mean flow; arrival-window anchor vs the plain-wave control |
| Thermal stress | `presets::multiphysics::run_thermal_stress` | thermal channel → thermal strain; free-bar expansion u(L) = α·g·L²/2 within 2% |
| Joule heating | `presets::multiphysics::run_joule_heating` | static-field E → thermal source channel; steady power-balance (source = outflux) |
| Species in flow | `presets::multiphysics::run_species_in_flow` | fdm3 channel → species advection-decay; outlet c matches exp(−k·L/u) within 5% |

## 0b. W15: static-fields Neumann faces, robotics, chiplet boards

| System | Run entry | Fidelity/anchors |
|---|---|---|
| Electrostatics (fixed) | `physics::electromagnetics::solve_static_field` + `FaceKind::Neumann` | the parallel-plate capacitor with Neumann side walls reproduces the EXACT discrete linear bridge (1e-7) and E = −V/L; the grounded-box sag (real physics) documented |
| Serial manipulator (CLASS) | `physics::robotics::step_manipulator` (N links) + `presets::robotics::run_manipulator_trajectory` | any serial arm from the user's link/joint/gain data: mass coupling, Coriolis, torque limits, soft stops, PD anti-windup; energy conserved, tracking < 1e-3 rad, stops bound |
| Heat conduction (CLASS) | `physics::thermal::solve_heterogeneous_conduction` + `presets::thermal::run_steady_conduction` | per-node k/q fields (the CAD/import contract) or box regions, per-face FixedValue/Insulated BCs; eliminated-Dirichlet + symmetric half-row operator, CG; source == sink flux exactly, linearity, spreader cooling, exact 1D bridge. Chiplets/heat sinks/reactor walls are fixtures, not presets |

---
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
`presets::turbine::default_grid_config` + `run_coupled_turbine`.

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
| Steam EOS | `physics::thermo::steam`: Clausius–Clapeyron saturation line (p_sat(T)/T_sat(p)), liquid/vapor/latent enthalpies, saturated vapor density — engineering model, ~5–10% vs IAPWS across 0–200 °C | real |
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
`physics::fluid::turbomachinery::simulate_compression_system` (or `solve_stage_stack`
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
| Coupling | real `CouplingManager` (nearest/trilinear transfer, relaxed links with target read-back, Jacobi-ordered exchange), `CoupledSimulation` multi-rate driver | staggered vs implicit both land on analytic; two-slab CHT converges to the exact interface profile (350 in a 400→300 rod) |
| Thermal transient (W11) | implicit time stepping on the same SOR operator; velocity-channel advection (per-node sampling) | Fourier-series match 0.2%; effective-Péclet profile 0.3% |
| Aeroacoustic-lite (W11) | convected wave equation with uniform mean flow | pulse arrival at c+u / c−u within 5%; zero-flow control |
| FSI-lite (W12) | `simulate_drag_body`: 6-DOF body ⇄ fdm3 via normalized Gaussian-smeared point force, upstream drag probe, preconditionable flow | smearing sum exact; zero-force purity; wind-tunnel terminal v_t ±10%; \|F_drag\|→mg ±5%; deterministic |
| CHT-lite demo (W12) | steady duct field drives thermal advection (W11 velocity channel) | outlet enthalpy ≈ source power ±35% (upwind diffusion); mean outlet T = 1D estimate ±5%; exact conduction control |

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
- Verified by 88/88 unit tests including analytic anchors for every
  physics model (W9: 80/80, W10: 87/87, W11: 88/88).
