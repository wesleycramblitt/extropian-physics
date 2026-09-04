# Architecture Review — modularity, shared core math, and hack inventory

Status: W18 review (docs only).  Two independent passes: a repo-wide core-math
duplication audit, and an adversarial review of the W16–W18 work (immersed
solids, the Heun fix, the field()-sync, the benchmark demos).

---

## 1. Verdict

The physics claims in W16–W18 are believable and better-verified than most of
the repo.  The ARCHITECTURE claims are not: the fdm3 module became MORE
self-contained with every wave, and the three most recent weeks were the
three most self-contained of all.  There are no crash-level blockers, but
three genuine hacks and a family of duplicated core math.

What the shared core already is: `numerics` (8 ODE integrators, CG/BiCGSTAB/
GMRES, matrix-free operator algebra, TimeStepper, Newton), `discretization`
(FdmLaplacianOperator, DiffusionStepOperator, raw stencils), `mesh`
(StructuredGrid family), `coupling` (channels, samplers, mapping).  It serves
~5 of ~15 solvers; fdm3/fdm/wave/static/elasticity re-implement the same
math privately.  The plan of record even names fdm3 "the module that never
joined the shared layer."

Load-bearing fdm3 exceptions (justified): collocated padded storage (vs
node-centered StructuredGrid), the coupled u/p SIMPLE system (not
expressible through the scalar-operator pattern), Central/Hybrid advection
(the shared layer has upwind only).

Not load-bearing (debt): the SOR Poisson kernel, the dt policy, the stencil
family, the index conventions, the BC ghost-fill.

---

## 2. The three hacks to remove first

### H1. The W17 Heun fix triplicates the pressure-projection kernel
`fdm3_integration.cpp` (Heun ~119-149 and CrankNicolson ~247-277) embeds a
verbatim copy of the solver's SIMPLE block (`fdm3_solver.cpp:247-256`) plus
`correct_velocity()` (`fdm3_pressure.cpp:81-101`): each does
compute_pressure_rhs → fill(p,0) → solve_pressure_poisson → central-difference
gradient correction, differing only in variable names.  Consequences:
- the integrators silently own the Poisson knobs (sor_omega, tolerance,
  ghosts, dt semantics) — future Poisson changes must be applied in THREE
  places;
- Heun/CN now run TWO full SOR solves per step at the tight OUTER tolerance
  unconditionally (the docs name SOR the dominant phase);
- the 2D fdm solver still has the plain trapezoid — the same bug class with
  no fix.

Fix: extract ONE fractional-step operator in fdm3_pressure:
`project_velocity(g, config, dt, mode)` + `apply_pressure_gradient(...)`;
the integrators become pure stage schedulers; the solver's SIMPLE block and
correct_velocity both call it.  (~½–1 day.)

### H2. The W18 field()-sync is a convention adapter with a silent gate
`fdm3_solver.cpp:214-226`: every step copies the full adapter into the grid
(u,v,w,p), gated by a silent dimension check — a caller who resizes the
adapter loses their edits with NO error; the copy is ~134 MB/step at 128³
even when nothing touched field(); and the p-copy is overwritten by the
SIMPLE solve in the same step (accidental semantics, not intended DX).

Fix: dirty flag set by the non-const `field()` accessor + ingest only when
dirty + hard ModelStatus error on mismatch.  (~½ day.)

### H3. The dt policy lives in benchmarks, and the solver can still blow up
`clamp_dt_from_cfl` (fdm3_solver.cpp:153-161) is advection-only; the
explicit-diffusion bound ν·dt/h² ≤ 1/6 exists only inside benchmark configs
(cases_fdm3.cpp), and the shared `numerics::TimeStepper` (which already
computes BOTH compute_cfl_u and compute_cfl_nu) has ZERO production callers
while four solvers hand-roll weaker clamps (fdm3, species, wave, fdtd).

Fix: fold the diffusive bound into the fdm3 clamp when adaptive_dt, warn in
validate() on ν·dt/h² > 1/6, and port the solvers onto TimeStepper (or
delete it — shared-but-unconsumed is a trust signal).  (~½ day.)

---

## 3. Duplication inventory (repo-wide)

| Topic | Sites | Verdict / fix |
|---|---|---|
| 7-point Poisson SOR | fdm3_pressure.cpp:36-70, fdm_pressure.cpp:17-45, static_fields.cpp:197-277, thermal_solver.cpp:246-346 (+19-point elasticity) | One `numerics::sor` engine; kills the live SIGN FORK (fdm solves Lap(p)=−rhs, fdm3 flips it — fdm3_pressure.cpp:27-32 vs fdm_pressure.cpp:32, a latent correctness hazard) |
| Explicit RK schemas | integrators.cpp:165, fdm_integration.cpp:111, fdm3_integration.cpp:218, rigid_body.cpp:227, manipulator.cpp:302 | Generic `explicit_rk(span, scratch, StageFn)` beside numerics; CFD keeps the projection between stages |
| "CrankNicolson" | integrators.cpp:202 (true implicit) vs fdm/fdm3 (explicit trapezoid = Heun) | Rename the CFD entry `Trapezoid`; share the skeleton between fdm and fdm3 |
| Trilinear interpolation | field_sampler.cpp:62-90, field_channels.cpp:28-56 (CHARACTER-IDENTICAL duplicate), surface_mapping.cpp:50-109, thermal adapter | One `interp::trilinear` primitive; the duplicate is a pure delete |
| Central-difference stencils | fdm3_internal.hpp:59-140, fdm_internal.hpp:53-103, wave_solver.cpp:255-366, static_fields.cpp:36-44, elasticity.cpp:234-262 — vs discretization/fdm/operators.hpp used by species/porous only | Storage adapter (padded ghost view over StructuredGrid) → shared stencil library (highest impact, highest effort) |
| BC ghost-fill | fdm_boundary.cpp:6-202 (per-edge boilerplate) vs fdm3_boundary.cpp:59-184 (generic face loop) | Shared BC spec + ghost-fill helper; 2D becomes an instance |
| SIMPLE corrector + diagnostics | fdm_solver.cpp:143-175 vs fdm3_pressure.cpp:81-101 + duplicated max_velocity/max_divergence/windowed convergence | Span-based corrector/diagnostics helper |
| CFL/dt policy | TimeStepper (zero callers) vs fdm3/species/wave/fdtd hand-rolled | Port onto TimeStepper or delete it |
| Box containment | static_fields.cpp:48-61, heterogeneous.cpp:252-261 (+ fdm3 signed distance) | mesh-level shape predicate |
| 2×2 eigen/expansion | cases_modules.cpp twice (plenum + PI) | one `bench/util` header |
| Index conventions (flat vs padded) | fdm3_cell_index, FDM3FieldData::index, the +1-shifted grid idx, the W18 ingest, add_body_force_to_rhs, benchmark locals | one `fdm3_layout` header owning both layouts + conversions (prerequisite for a real adapter view) |

---

## 4. Smaller hacks / smells (minor)

- `velocity_under_relaxation` is semantically overloaded as the Heun stage
  relaxation: lowering SIMPLE damping silently changes the integrator scheme,
  and classical Heun (s=1) becomes unrepresentable.  The commit message
  presents the two layers without isolating them — TEST the projection-only
  variant; if stable, delete the perturbation (and measure the temporal
  order of the fixed Heun — no dt-sweep temporal-order test exists anywhere;
  B1 measures only spatial order).
- The B5/B6 drag measurement reimplements solver internals through the
  two-way adapter (order-sensitive; will silently rot under H2 refactors).
  Expose `FDM3Solver::last_solid_reaction(mask)` or a named helper in
  fdm3_sources.
- The W16 "collocated projection cancels uniform interior body forces" is
  documented but unquantified and unfixed: add a regression that MEASURES the
  bulk-acceleration deficit vs the applied force, and ticket the real fixes
  (pressure-consistent force splitting f' = f − ∇φ, Δφ = ∇·f; or
  Rhie-Chow).  The label "classic pressure-velocity decoupling" is imprecise
  — this is a zero-frequency mean-mode phenomenon, not checkerboarding.
- Ghost cells are never refreshed between integrator stages (fdm3_integration
  never calls apply_boundary_conditions mid-step) — FE/RK4 read stale ghosts
  at boundary-adjacent rows too.  The standard explicit-CFD fix (stage BC
  refresh) deserves the primary experiment before trusting the
  under-relax+project band-aid.
- fdm3_sources: geometry/fraction recomputed three times
  (mask/penalty/freeze); "first match wins" solids union is order-dependent
  (use max-union); the freeze blends the SMOOTH fraction while the
  benchmark masks use hard containment — the concrete mechanism of the 2.4×
  drag finding; the blend clamp [0,1] is silent.
- benchmarks: private-header include (engine_internal.hpp); rel_l2's silent
  min-size truncation.

---

## 5. Priority plan

| # | Item | Effort | Gate |
|---|---|---|---|
| P1 | Shared fractional-step projection operator; integrators become stage schedulers (H1) | ½–1 d | fdm3 tests + the W17 regression stay green |
| P2 | dt layer: diffusive bound in the clamp + validate() warning; TimeStepper adoption (H3) | ½ d | fdm3_solver adaptive_dt tests |
| P3 | field()-sync hardening: dirty flag + hard error + p-copy decision (H2) | ½ d | W18 sources tests (real freeze semantics) |
| P4 | numerics::sor engine (kills 4 copies + the sign fork); single trilinear; bench 2×2 helper | 1 d | all SOR-consumer tests + TGV/thermal anchors |
| P5 | W16-transparency quantifying regression + fix ticket | hrs–2 d | new test records the deficit number |
| P6 | last_solid_reaction() / de-brittle B5/B6; stage-BC-refresh experiment; order test for Heun | 1 d | bench smoke rows unchanged |
| P7 | fdm3_layout index header, then the stencil/BC consolidation onto discretization | incremental | module boundary doc updated |

Net guidance: fix P1–P3 now (they are the hacks), P4 next (pure
consolidation + a live hazard), then let P5–P7 proceed as the module matures.
Write the doc's "module boundary" section as the controlled-exception
statement: collocated storage, coupled SIMPLE, and central/hybrid advection
are load-bearing; SOR, dt policy, stencils, and interpolation are debt that
migrates to the shared core.
