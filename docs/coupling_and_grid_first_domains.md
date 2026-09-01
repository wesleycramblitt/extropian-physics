# Real Coupling Exchange (Phase H-lite) + Phase I Grid-first Domains (W10–W11)

Companion to `modular_solver_architecture.md` §11/§12. Product doctrine:
80–90% engineering-grade accuracy, simple physics coupled together.

## 1. Phase H-lite — real coupling exchange

The placeholder `CouplingManager` (headerless, `transfer(nullptr,…)`) is
replaced by a real, channel-based exchange:

- `coupling::SurfaceMapper` — `transfer_nearest` (scattered point sets,
  brute-force O(n·m), configuration-size surfaces) and `transfer_trilinear`
  (structured-grid → structured-grid, physical-coordinate trilinear,
  out-of-bounds → NaN + warning). `InterpolationMode { Nearest, Trilinear }`.
- `coupling::CouplingManager` — domains register read channels
  (`IScalarField3D`/`IVectorField3D`) and point-wise write sinks via
  callbacks (the manager never owns solver state). `CouplingLink` carries
  probe points, interval gating, under-relaxation `(0,1]`, and implicit
  sub-iterations. Exchange modes: explicit (staggered, `exchange(t)`) and
  implicit (`exchange_implicit(t, tol)` — sub-iterate until the max probe
  delta < tolerance). Relaxation bookkeeping is per-probe
  `(1−ω)·prev + ω·sampled` (no target read-back exists — documented).
- `coupling::CoupledSimulation` — multi-domain, multi-rate driver: each
  domain steps at its own dt under a common macro-step clock (sub-cycling);
  staggered or implicit exchange per macro-step; `RunReport` with exchange
  and convergence accounting.
- **Acceptance test** — a linear 2-ODE system partitioned across two fake
  domains, coupled through probe links: staggered (dt→0, first-order) AND
  implicit (relaxed fixed-point) both land on the analytic solution
  (1% / 1e-3 targets); multi-rate with 2× dt also converges (2%).

## 2. Phase I domains — grid-first, one module each

| Domain | Namespace | Solver | Verification anchors |
|---|---|---|---|
| Thermal | `exd::physics::thermal` | steady 3D conduction/advection (SOR, central 7-pt Laplacian + first-order upwind), fixed/insulated faces, uniform source | fixed-wall linear profile; parabolic source profile; advection shifts + Péclet warning; insulated-flat test; 3D z-linear sanity |
| Acoustics | `exd::physics::acoustics` | scalar wave equation, explicit leapfrog, pressure-release walls, CFL-adaptive dt (clamp + warning) | 1D plane-wave period (1%); 3D box eigenfrequencies (1,1,1) & (1,1,2) (~0.1–1%); energy bounded over 10 periods; CFL guard |
| Structural | `exd::physics::structural` | static linear-elasticity displacement form (Navier–Cauchy, 19-pt stencil, ghost-Robin free faces, component-wise SOR, thermal strain via optional temperature channel, masked +z traction) | uniaxial ν=0 exact; confined oedometer exact (ν=0.3); thermal column exact (confinement factor (3λ+2μ)/(λ+2μ)); cantilever qualitative + exact load linearity |
| Particles | `exd::physics::particles` | Lagrangian cloud over sampled `IVectorField3D` channels, shared integrator (RK4), deterministic lattice spawn | ballistic exact; drag terminal velocity; channel advection analytic; number conservation + determinism |
| Chemistry | `exd::physics::chemistry` | 0D well-stirred reactor, mass-action + Arrhenius, shared integrator, non-negativity clamp with warning | first-order decay exact; reversible equilibrium K; second-order analytic; Arrhenius scaling ratio (~0.1%); mole conservation + determinism |

Channels: thermal (temperature) and acoustics (pressure) export
`IScalarField3D` via `make_scalar_grid_field`; structural exports
displacement `IVectorField3D`. Particles **consume** a flow channel;
chemistry and particles are lumped/0D modules exempt from the per-domain
channel rule (engine/compressor precedent) — recorded in the headers.

## 3. Documented envelopes (where it is 80 %, where it is wrong)

- Thermal: first-order upwind advection is numerically diffusive
  (Péclet > 1 warning); steady-state only in v1.
- Acoustics: pressure-release (soft) walls only; leapfrog dispersion
  grows with mode count (verified ≤ 2% on the tested modes).
- Structural: free-surface ghost treatment is first-order at edges/corners
  (mirror + face corrections, tangential corrections single-face only);
  quantitative beam-bending accuracy at coarse z resolution is deferred to
  a future boundary-scheme upgrade — the suite pins the EXACT regimes
  (ν=0 uniaxial, confined compression, thermal column) and the qualitative
  + linearity properties of bending. SOR stagnation guard stops
  unrealizable 1e-10 chases with a warning.
- Particles/chemistry: integrator-tolerance-bound (RK4 default).

## 3b. W11 — transient thermal, CHT-lite, mean-flow acoustics

- `thermal` gains transient implicit stepping: `ThermalState` +
  `init_thermal_state`/`advance_thermal`/`simulate_thermal`, one
  backward-Euler-style step per advance (unconditionally stable), with the
  time level held FIXED on the RHS (an in-place RHS cancels the time term
  at the SOR fixed point — caught by the Fourier-series test).
  `set_temperature_point` pins interface nodes for link-driven coupling;
  `make_temperature_channel` exposes a trilinear read channel.
- Velocity-channel advection: `ThermalConfig.velocity_channel`
  (IVectorField3D) sampled per node (CHT-lite); falls back to
  body_velocity out of bounds.
- `acoustics` mean flow: linearized convected wave equation
  `(dt + u·grad)^2 p = c^2·lap p`; CFL keyed on c + |u|; warning at
  |u| >= c. Verdict: pulse arrival at c+u / c-u within 5%.
- CouplingManager exchange is JACOBI-ORDERED (all sources sampled +
  targets read back BEFORE any write) with optional target read-back
  relaxation.  This makes the interface fixed point independent of link
  write order — the two-slab CHT acceptance test pins it at the
  flux-matched 350 K (a sequential per-link order biased it to 377 K).

## 3c. W12 — FSI-lite drag body + CHT-lite channel demo

`exd::physics::coupling::simulate_drag_body` (include
`coupling/drag_body_solver.hpp`): a 6-DOF rigid body under gravity in fdm3,
two-way coupled through a Gaussian-smeared point force with DISCRETE
normalization Σwᵢ = 1 (so the injected fluid momentum equals the negated
body drag by construction — unit-pinned).  Quadratic drag on the RELATIVE
velocity, F = −½ρ·C_dA·|v−v_f|·(v−v_f); β-under-relaxed; optional flow
preconditioning (run the duct to steady before releasing the body).

**Acceptance (tests/unit/coupling/drag_body_fsi_test.cpp):**
- smearing normalization exact (Σ ρ·a·dV ≡ F, 1e−12);
- zero-force purity: C_dA=0 → v(t) = −g·t within 1%, the fluid never
  disturbed and the body drifts nowhere;
- wind-tunnel terminal: preconditioned laminar duct (μ=0.2 Pa·s, Re≈1,
  conditioned ~15 s), body falling ALONG the duct: |v_x − u| → v_t =
  sqrt(2mg/(ρ·C_dA)) within 10%, |F_drag| → mg within 5%, the drag probe
  sees ~zero vertical flow (centerline symmetry), deterministic across
  identical configs.

**Physics lessons recorded (each measured during the wave):**
1. POINT-FORCE SELF-INDUCTION: sampling v_f at the body CoM measures the
   fluid the body's own reaction just pushed; once the local v_f catches
   the body speed the drag flips sign and the two-way loop rings up
   (termed "exponential").  v1 senses the flow `sample_lead·ε` AHEAD of
   the CoM along the relative motion (the attacker side); wake feedback AT
   the body surface is W13 (body-fitted sampling).
2. CLOSED-BOX LEDGER: in a closed inviscid box the smeared reaction builds
   a deep entrainment column; once it reaches a wall the pressure returns
   momentum, so m_b·v_b + ∫ρu dV = m_b·g·t carries a wall-impulse
   term (measured deficit = blob-mass × probe velocity).  The clean
   integrated anchor is therefore the OPEN box (wind tunnel), where the
   reaction momentum leaves the domain.
3. EXPLICIT-DIFFUSION LIMIT: fdm3's Laplacian is explicit; dt must satisfy
   ν·dt·Σ(1/h²) < 0.5 (dt=0.02 at μ=0.05 on h=0.05 → NaN; dt=0.01 OK).
4. ADVECTION-DOMINATED SOR FLOOR: the thermal SOR (ω=1.5) stalls at
   ~1.7e−5 residual on the upwind advection matrix regardless of the
   requested tolerance — seeking 1e−6/1e−8 silently burns the full sweep
   cap (22 s); target the achievable floor (1e−4) when O(1e−2) checks
   suffice.
5. OUTLET CHOKE: near the fixed-pressure outlet the duct flow decelerates,
   the drag undershoots, and the body can reverse; keep the run clear of
   the outlet (longer duct).

## 4. Cross-domain demos (tests, not apps)

- CoupledSimulation staggered-vs-implicit linear agreement (H acceptance).
- TWO-SLAB CHT (W11 acceptance): two transient thermal domains on
  `CoupledSimulation`, interface temperature exchanged through relaxed
  read-back links → exact joined profile 400 → 350 → 300.
- Structural thermal-expansion coupling-lite: a temperature channel drives
  the displacement field directly.
- Particle channel advection: a structured velocity grid (or any
  `IVectorField3D`) carries the cloud.
- CHT-lite channel demo (W12): a steady inlet/outlet laminar duct solves
  first; the field is wrapped as a velocity channel (StructuredVectorGrid
  on the exact cell-center lattice) and a thermal domain on the same lattice
  solves with it — outlet-plane enthalpy flux within 35% of the source
  power (the first-order upwind's numerical diffusion and the hot-core ×
  fast-core covariance inflate the node-weighted integral; documented),
  the mean outlet temperature matches the 1D advective estimate within 5%,
  and a zero-flow control reproduces the exact 1D conduction quadratic.
  One-way by physics (incompressible — no Boussinesq).
- Full conjugate heat transfer (thermal ↔ fdm3 velocity) and aeroacoustics
  (acoustics ↔ fdm3) as `CoupledSimulation` demos are the next milestone.
