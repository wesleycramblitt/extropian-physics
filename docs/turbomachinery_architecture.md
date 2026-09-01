# Turbomachinery Architecture — Mean-Line Stages, Maps, Plenum, System Driver (W9)

Companion to `modular_solver_architecture.md` §10 (Phase G re-scope) and
the BEM reference doc pattern. Product doctrine: 80–90% engineering-grade
accuracy, simple physics coupled together, fast and easy to use.

## 1. What this is

The compressor use case, solved product-agnostically. **No `compressor`
app module exists.** The same axial mean-line physics produces compressor,
turbine, fan, and pump senses from blade geometry alone — the direction of
energy transfer is geometry-emergent, never a mode flag:

```
Δh₀ = u·(c_w2 − c_w1)        specific work (J/kg); sign is emergent
τ   = ṁ·r_m·(c_w2 − c_w1)    shaft torque (momentum — valid at ω = 0)
```

| Module | Namespace | Entry |
|---|---|---|
| Stage | `fluid::turbomachinery` | `solve_stage(StageConfig, StageInlet, ω, ṁ, IEos, ModelStatus)` |
| Stack | `fluid::turbomachinery` | `solve_stage_stack(StageStackConfig, …)` — N stages, propagates (p0, T0, c_θ) |
| Map | `fluid::turbomachinery` | `solve_operating_map(...)` / `sample_operating_map(...)` |
| Plenum | `fluid::lumped` | `plenum_derivative` / `step_plenum` (Greitzer, isentropic plenum) |
| Driver | `fluid::turbomachinery` | `simulate_compression_system(...)` (shaft+plenum+governor+CSV) |
| Polytrope | `thermo::polytropic` | stagnation-family relation primitives |

## 2. Stage physics (total-state bookkeeping)

All stage relations operate on **total** states (T0, p0); static states are
recovered from the resolved velocity field. Mixing static and total states
is an O(M²) error (~16 % in π at M = 0.5).

1. Inlet density-velocity fixed point: `c_a = ṁ/(ρ·A)` coupled with the
   static-state chain `T1 = T01 − c1²/(2cp)`, `p1 = p01·(T1/T01)^(γ/(γ−1))`,
   `ρ1 = eos.density(p1, T1)` — direct iteration, 1e-10 tolerance, 50-iter
   guard → `ModelStatus` failure.
2. Triangle geometry at mean radius: `u = ω·r_m`; `c_w1 = c_a·tan(α1)`;
   `w_u1 = c_w1 − u`; `w_u2 = c_a·tan(β2)`; `c_w2 = u + w_u2`.
3. Euler work and polytropic exit: `T02 = T01 + Δh₀/cp`;
   compression `π = τ^(γ·η_p/(γ−1))`, expansion `π = τ^(γ/((γ−1)·η_p))`
   via `thermo::polytropic` (compressor and turbine polytropes are NOT
   reciprocal; verified by test).
4. Static exit: `T2 = T02 − c2²/(2cp)`, `p2 = p02·(T2/T02)^(γ/(γ−1))`.
5. Choking on **relative** Mach at the rotor LE `M_rel = √(c_a² + w_u1²)/a1`
   (axial Mach is the wrong proxy — the LE relative flow is u-dominated).

**Validity envelope (documented promise — you know where it is 80 % and
where it is wrong):** `M_rel < 0.7` (warning above), hub/tip ratio > 0.5
(warning below), single-stage π < ~1.5, axial-only (radial stages need a
different loss model — future variant). Warnings ride the `ModelStatus`
channel; results stay valid.

## 3. Multi-stage and maps

- Stack: constant ṁ/ω, inter-stage propagation of (p0, T0, c_θ) — residual
  swirl IS the next stage's inlet swirl; total π = Ππᵢ; inter-stage reheat
  falls out naturally (‖π₂ < π₁ for equal work on a hotter inlet — tested).
- `solve_operating_map`: rectangular ω × ṁ sweep (dimensionless scale:
  corrected-flow/corrected-speed compatible), invalid points = NaN; surge
  line = argmax-π per speed line (classical stall-line surrogate — real
  surge lines come from rigs/CFD); choke line = first choked point per
  speed line. `sample_operating_map`: bilinear with clamp + nearest-valid
  fallback + warnings. Any map — computed OR rig-test data — plugs into the
  same sampler (the `TableLookup` variant pattern).
- Known mean-line artifact, documented: the choked region enters from LOW
  flow (u-dominated M_rel) and may clear at higher flow; the map keeps the
  first-choke semantics per speed line.

## 4. Greitzer plenum

`fluid::lumped::plenum` — the classical two-state surge cell:

```
I = L/A (duct inertance)
ṁ̇ = (Δp_c(ṁ) − (p − p_amb)) / I
ṗ = (a²/V)·(ṁ − ṁ_t(p)),   a² = γ·R·T_p,  T_p = T_amb·(p/p_amb)^((γ−1)/γ)
```

Isentropic plenum assumption documented in the header. Compressor and
throttle characteristics are injected functions — the plenum itself is a
generic gas-network node (engine manifolds, gas networks, HVAC later).

Stability verification uses the true two-state Jacobian (trace/det
conditions weighted by the Greitzer B parameter) — the naive "compressor
slope vs throttle slope" comparison is WRONG in the shallow-throttle
regime and is tested against time-march on both sides, plus a classic
surge limit-cycle case (cubic compressor characteristic, large volume:
sustained oscillation, collapsed mean flow).

## 5. Compression-system driver

`simulate_compression_system(CompressionSystemConfig, ModelStatus&)`:
thin (~460 lines) composition — DC motor (Phase D) + compressor stage
stack (this module) + Greitzer plenum (this module) + optional PI speed
governor (Phase C) modulating the throttle gain; RK4 via the shared
integrator module; CSV streaming via `io::CsvSeriesWriter`. Conventions
mirror the engine: governor and motor update exactly once per step.

- Compressor characteristic = `solve_stage_stack` output (p0_out − p_amb);
  identical physics path to the standalone plenum.
- 1e-9 flow floor recovers blocked-flow pressure rise at ṁ = 0 so a drive
  can spin up from rest.
- Energy bookkeeping: rectangle rule with step-average ω (exact for
  constant torques) → `energy_balance_closed` checks
  `|W_motor − W_compressor − ΔKE|/W_motor < 5 %`.

## 6. Verification (tests, all analytic/hand-anchored)

- Exact Euler work / torque / power on synthetic triangles (closed-form
  re-implementation in the test).
- Zero absolute-swirl-change ⇒ zero work; geometry-emergent sense reversal:
  an undo stage returns T0 exactly, Δh₀ exactly cancels, total π < 1 with
  Δs > 0 (η_p = 1 ⇒ π = 1 exactly reversible) — the product-agnosticism
  proof. Compressor and turbine polytropes are explicitly NOT reciprocal.
- Closed-form static/total/continuity/energy closure at M_axial ≈ 0.4–0.5.
- Stack: single ≡ single; reheat effect; determinism.
- Map: echoes direct solves at nodes; surge argmax; relative-Mach choke
  boundary; determinism; clamp/fallback warnings.
- Plenum: derivative signs; analytic fixed point; Jacobian stability both
  sides vs time-march; surge limit cycle; validation.
- System: motor spin-up settles at the operating point (ω, π ≈ 1.08, ṁ≈0.14,
  energy closure 6e-7); governor regulates speed 0.0002 % with clamps held;
  turbocharger balance — compressor and turbine stages from ONE code path
  self-balance on a shaft (ω* within 0.1 % of the geometry identity,
  torque balance to 2e-16, Δs > 0).

## 7. Boundaries (honest margins)

- Mean-line, single radius: poor where hub/tip is extreme or radial work
  distribution matters (envelope warnings enforce the limits).
- Losses are a single per-stage polytropic efficiency; no profile/secondary/
  tip-clearance correlations yet (calibratable later against maps or CFD —
  real surge lines and efficiency surfaces come from rigs, which the map
  sampler already accepts).
- No field channels by design: this is a lumped machine model like the
  engine. Compressible-grid coupling is an unowned future capability.
- Real-gas EOS for refrigerants etc. lands behind the existing `IEos` seam
  (steam EOS precedent) without touching the stage code.
