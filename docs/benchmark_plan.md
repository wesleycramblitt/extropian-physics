# fdm3 Benchmark Plan — accuracy vs analytical/public references, performance & cost

Status: plan (W18, docs only).  One benchmark per solver use case, each anchored to an
analytical solution or a published reference dataset, each with explicit metrics,
acceptance tiers, and a performance protocol.  Implementation is phased (Section 5).

---

## 0. Principles

1. **Every benchmark measures two independent axes**: accuracy (against the reference) and
   cost (run time, memory, phase shares).  The two meet in an *accuracy-per-second* report.
2. **Every benchmark is reproducible**: fixed config hash + grid + dt + scheme + integrator,
   recorded as one CSV row per run with the git commit; serial, warmed-up, median-of-N timing.
3. **Tiers keep CI cheap**: smoke tier (seconds) can run in CI; the full sweeps are manual
   (`benchmarks/run_all.sh`).
4. **The collocated-grid caveats are documented per benchmark**, not hidden:
   - the pressure projection cancels uniform interior body forces in the fully developed
     state (W16) → drag must be measured via the kinematic-freeze reaction force or the
     penalty's transient response, not a developed-state uniform force;
   - Upwind is first-order (expect order-1 slopes deliberately);
   - 2D published values (cavity, cylinder, step) are compared against a thin-z slab with
     z-symmetry; the midplane is 2D-like, the third dimension adds documented deviation;
   - the SOR Poisson tolerance interacts with accuracy floors → sweep it.

Solver surface the plan relies on (verified): BC types Inlet/Outlet/Wall/Symmetry/
Periodic/FixedPressure; advection Central/Upwind/Hybrid; time integration
ForwardEuler/Heun/RK4/CrankNicolson (Heun stabilized in W17); `adaptive_dt` + `cfl_target`;
per-cell body forces (`set_body_force`) — the manufactured-source hook; immersed solids
(penalty + kinematic freeze, W16); Boussinesq buoyancy; thermal channel (`solve_thermal`).
The solver is serial (no OpenMP/TBB) — benchmark throughput as-is, report the SOR cost model.

---

## 1. Use-case coverage map

| Solver use case | Benchmark family | Reference type |
|---|---|---|
| Core momentum/transport accuracy, schemes, integrators | B1 MMS, B2 Taylor–Green | Exact (MMS source; closed-form decay) |
| Boundary treatment (walls, in/out, periodic, fixed p) | B3 Cavity, B4 Channel | Published tables; exact profiles |
| Body forces / actuator momentum sources | B1 MMS (enters via body force), B2 | Exact |
| Immersed solids, blockage, drag/FSI | B5 Stokes sphere, B6 Schäfer–Turek cylinder | Asymptotic / public CFD benchmark |
| Boussinesq buoyancy, natural convection | B7 de Vahl Davis, B8 Rayleigh–Bénard onset | Published Nusselt tables; linear stability |
| Thermal-fluid coupling (forced convection) | B11 Graetz / Nu∞ | Series solutions, classic Nu values |
| Separated flow | B9 Backward-facing step | Gartling / Eça–Hoekstra reference data |
| Wall-bounded external flow | B10 Blasius layer, B3D-6 yawed plate | Exact similarity (2D + 3D) |
| 3D geometry/spanwise physics | B3D-1…B3D-8 (Section 2A) | Exact series / published 3D tables / stability onsets |

---

## 2. Accuracy benchmarks

### B1 — Method of Manufactured Solutions (core order verification)

- Case: 3D periodic box (Periodic BCs on all faces), smooth exact solution satisfying
  continuity (e.g. trig-polynomial product), residual source terms injected through
  `set_body_force` (the acceleration API is exactly the manufactured-source hook).
- Sweep: nx = 8, 16, 32, 64 (dx halving) × {Central, Hybrid, Upwind} ×
  {FE, Heun, RK4, CN}.  dt from a fixed CFL (0.1) so spatial errors dominate; separately a
  dt-sweep at fixed dx for the temporal order.
- Metrics: global L2 error vs dx → observed spatial order (expect ≈2 Central/Hybrid,
  ≈1 Upwind); L2 vs dt per integrator (expect ≈2 Heun/RK4/CN, ≈1 FE).  This is also the
  permanent regression guard for the W17 Heun fix.
- References: Roache, *Code Verification by the Method of Manufactured Solutions* (2002);
  Salari & Knupp, SAND2000-1444 (2000).
- Tiers: smoke = order slope within ±20% at 3 refinement levels on 16/32/64;
  full = ±10% + absolute error floors at 64.

### B2 — Taylor–Green vortex (exact transient, + flagship perf case)

- Case: periodic box, u = sin x cos y cos z · e^(−νt) (cyclic), ν from config.
- Metrics: kinetic-energy decay E(t) vs the closed form (effective ν within a band);
  pointwise L2 error growth vs t (up to t = 5/ν); the dissipation-peak timing for the
  sustained variant.  The existing unit test is the smoke anchor; the benchmark refines
  grids and adds the performance sweep (see Section 3).
- Tiers: energy decay within 3% at 64³.

### B3 — Lid-driven cavity (walls + steady convergence)

- Case: 1³ box, no-slip bottom/sides, lid u = 1; Re = 100, 400, 1000 (ν = 1/Re).
- Reference: Ghia, Ghia & Shin, *High-Re solutions for incompressible flow…* (JCP 48,
  1982) — centerline u(y), v(x) tables; 3D: Albensoeder & Kuhlmann (2005).
- Metrics: L2 error of the centerline profiles vs the tables; corner-eddy extents
  (qualitative at coarse, quantitative at fine); steady-iteration count vs the Poisson
  tolerance.
- Tiers: Re=100 centerline L2 < 5%; Re=1000 < 15% at the finest grid; eddy structure
  qualitative.

### B4 — Channel flows: Poiseuille + Couette (exact, incl. mass flow and shear)

- Case: thin-z slab (z = 1 cell + Symmetry), pressure-driven (FixedPressure in/out) and
  moving-wall (shear-driven) variants; also an Inlet-driven entry-length transient.
- References (exact): u(y) = (Δp/2μL)(hy − y²), Q = h³Δp/(12μL) per unit width;
  wall shear τ_w = 6μQ/h²; Couette u = Uy/h; entry length L_ent ≈ 0.05·Re·h.
- Metrics: L2 profile error vs dx (order re-check); mass-flow relative error; τ_w relative
  error; recovered-vs-imposed Δp; entry-length estimate.
- Tiers: Q and τ_w within 2% at the finest grid; the exact-profile L2 slope confirms the
  spatial order at a second geometry.

### B5 — Stokes / Oseen drag on a sphere (immersed solid + FSI drag path)

- Case: small sphere (immersed box-approx or sphere mask) in slow uniform flow;
  Re = 0.1–10 (based on D), penalized (transient) or kinematic-frozen (steady) per the
  W16 caveat; drag measured from the immersed-solid reaction force.
- References: Stokes F = 6πμaU (Re→0); Oseen Cd = (24/Re)(1 + 3Re/16);
  Schiller–Naumann Cd for the wider range.
- Metrics: Cd(Re) vs the asymptotics; deviation vs a/D resolution (a ≥ 4 cells ok);
  cross-check the standalone `drag_body_fsi` path against the pure-fluid value.
- Tiers: Re = 0.5 within 5% at a ≥ 4 cells.

### B6 — Schäfer–Turek cylinder (flagship PUBLIC benchmark, immersed/FSI)

- Case: channel 2.2 × 0.41 m, cylinder D = 0.1 at (0.2, 0.2) as an immersed solid
  (kinematic freeze — the W16 hard-blockage path); thin-z + Symmetry; Re = 20 (steady)
  and Re = 100 (shedding).
- Reference: Schäfer & Turek, *Benchmark computations of laminar flow around a cylinder*
  (1996): Re=20 → Cd = 5.5795, Cl = 0.0106, Δp = 0.1175 (12-digit reference values);
  Re=100 → St ≈ 0.30, Cd/Cl amplitude bands from the literature.
- Metrics: Cd/Cl/Δp point values; Strouhal from the lift-spectrum FFT; grid convergence
  (6 → 12 → 24 cells across D); documented deviation for the slab-2D vs published 2D.
- Tiers: Re=20 Cd within 10% at fine grids; St within 5% of 0.30; smoke = qualitative
  shedding at Re=100.

### B7 — Natural convection: de Vahl Davis square cavity (Boussinesq)

- Case: 1×1 slab, ΔT hot/cold vertical walls, adiabatic horizontals, Pr = 0.71;
  Ra = 10³ → 10⁶ (body force via the W16 Boussinesq sampler).
- Reference: de Vahl Davis, *Natural convection of air in a square cavity: a benchmark
  numerical solution* (IJNMF 3, 1983): Nu_mean, Nu_max, u_max, v_max and positions.
- Metrics: mean-wall Nu vs the published bands; u_max / v_max / positions;
  midplane T-profile; Ra × grid sensitivity (Nu within a documented band per Ra).
- Tiers: Ra = 1e3/1e4 Nu_mean within 5%; Ra = 1e5 within 10%; smoke = qualitative rolls.

### B8 — Rayleigh–Bénard onset (linear-stability anchor)

- Case: box heated from below (no-slip top/bottom, Symmetry/Periodic sides), Pr = 0.71,
  small sinusoidal T-perturbation; observe growth/decay.
- Reference (exact): critical Rayleigh number Ra_c = 1707.762 (linear stability of the
  conduction state between rigid conducting plates).
- Metrics: the bracket [Ra_lo, Ra_hi] that flips the perturbation-growth sign, bracketing
  Ra_c within the resolution-dependent width; steady roll wavelength ≈ 2√2·h (qualitative).
- Tiers: bracket within ±15% of Ra_c at 32³; smoke = sign flip exists.

### B9 — Backward-facing step (separation/reattachment)

- Case: channel with a step (immersed frozen block at the inlet height drop, or domain
  geometry), Re = 800.
- Reference: Gartling (1990) and Eça–Hoekstra QNET-CFD data: primary reattachment
  x1/H ≈ 6.26 (literature band ±0.1), secondary recirculation extent, velocity profiles at
  x/H = 4, 6, 10.
- Metrics: x1/H value and trend vs grid; profile-L2 vs the reference sections.
- Tiers: x1 within 10% at fine grids; smoke = primary recirculation present.

### B10 — Blasius boundary layer (wall-bounded external flow)

- Case: flat plate (Wall at y=0 downstream of a Symmetry leading section), uniform inlet,
  low Re_x laminar.
- Reference (exact similarity): Blasius profile (tabulated f′), Cf = 0.664/√Re_x,
  δ*/x = 1.7208/√Re_x, θ/x = 0.664/√Re_x.
- Metrics: Cf(x) pointwise error away from the leading edge; profile vs the tabulated
  Blasius function; δ*/θ ratio.
- Tiers: Cf within 15% for x beyond the entry-corner transient; smoke = δ*/θ ≈ 2.59.

### B11 — Forced-convection thermal entry (Graetz) and Nu∞

- Case: plane channel with the thermal module: constant-wall-temperature and constant-flux
  variants; also a round-tube config if available.
- References: fully developed Nu∞ = 7.54 (channel, const T) / 8.23 (const q′);
  3.66 / 4.36 (tube); entry curve vs the Graetz series solution (public tables).
- Metrics: asymptotic Nu error; Nu(x) entry curve L2 vs the series; bulk-temperature
  recovery.
- Tiers: Nu∞ within 5%; entry curve qualitative at coarse grids.

---

## 3. Performance & cost methodology

1. **Harness**: a `benchmarks/` target (not unit tests) — one binary per family or a
   registered-case runner; CSV result rows.  Timing: dry run (warm caches), then median of
   ≥ 5 timed repetitions; report per-phase clocks.
2. **Per-phase breakdown** (using the existing step structure): RHS evaluation, Poisson SOR
   (per-iteration), velocity correction, BC pass, field extraction → the SOR share is the
   first thing to report (it dominates).
3. **The record set** (each with runtime + phase shares + memory):
   - throughput clips: 32³ → 256³ channel/cavity steady runs: cells/s, steps/s;
   - SOR iteration count vs grid size and tolerance (the cost model: iterations ~ how they
     scale with N and tolerance);
   - memory: live working set and bytes/cell (grid doubles + old-state copies + field
     adapters are all measurable);
   - accuracy-per-second Pareto: for B2/B3/B4 sweeps plot (runtime, error) per grid/scheme/
     integrator — this is the deliverable for the "schemes and integrators" decisions;
   - scheme/integrator comparison on the SAME case: Central vs Hybrid vs Upwind;
     FE vs Heun vs RK4 vs CN (W17 made Heun viable — benchmark its cost/accuracy against
     RK4 at the same dt; expect it to win on phase count).
4. **Reproducibility**: fixed machine + CPU pinning, serial, compiler/arch recorded in the
   row; `run_all.sh` writes `results/` with commit-hash metadata; a `results/README` logs
   the machine fingerprint.

---

## 4. Acceptance summary (per tier)

| Tier | Runtime budget | Runs where |
|---|---|---|
| smoke | < 60 s total | CI (ctest, tagged) |
| full | hours | manual `benchmarks/run_all.sh` |
| reference | overnight | manual, dedicated runs recorded to CSV |

---

## 5. Phased implementation roadmap

- **Phase 1 — core (B1, B2, B4)**: harness + MMS + TGV + channel.  All three run on the
  existing APIs (body force, Periodic BCs, FixedPressure); fastest to deliver and they pin
  the order/scheme/integrator verification the other phases build on.
- **Phase 2 — boundary/steady (B3 cavity, B9 step, + B3D-1, B3D-2, B3D-7)**: steady-state
  convergence, wall accuracy; the 3D family starts here with the exact duct series
  (B3D-1, the cheapest genuine-3D anchor) and the cubic cavity; needs longer runs and the
  SOR-tolerance sweep.
- **Phase 3 — immersed/FSI (B5 sphere, B6 cylinder, + B3D-3, B3D-4, B3D-5)**: reuses the
  W16 solids; the flagship public anchors (Schäfer–Turek 2D, the DFG/FeatFlow 3D band,
  and the Barkley–Henderson mode-A/B onsets); drag measurement per the W16 caveat; the 3D
  wake cases (B3D-4/5) are the most compute-hungry → ran last in the phase.
- **Phase 4 — thermal (B7, B8, B11, + B3D-6, B3D-8)**: Boussinesq + thermal channel; de
  Vahl Davis, Fusegi's cubic cavity, the Rayleigh–Bénard linear-stability anchor, and the
  exact yawed-plate 3D similarity solution (B3D-6, the cheapest 3D thermal-adjacent
  anchor); needs the T-field coupling wiring.
- **Phase 5 — performance harness + registry + CI smoke tags**: Pareto reports, CSV
  registry, `results/README` machine fingerprint, CI wiring for the smoke tier.

Each phase: implement → verify anchors → record rows → commit with the results summary.
Phase 1 is the recommended starting point (fastest complete accuracy story + the perf
harness skeleton).
