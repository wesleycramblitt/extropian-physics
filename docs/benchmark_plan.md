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

## 2A. 3D coverage — every 2D-anchored case also gets a genuine-3D path

**Rule for every case whose canonical reference is 2D**: two obligations —
(1) the **periodic-extrusion recovery check** (thin span, Periodic z BCs: the
3D solver must reproduce the published 2D value on the midplane inside the
documented slab-deviation band — this is itself a 3D validation with a public
anchor), and (2) a **genuine-3D benchmark** with its own reference, from the
table below.  Where a case's 3D reference is only a single-group dataset
(DFG-3D cylinder, 3D BFS DNS) we treat it as a **band anchor** (± published
spread), never as digit-exact like the 2D tables.

| 2D-anchored case | Extra 3D anchor(s) | Type |
|---|---|---|
| Channel B4 | **Rectangular-duct exact double-series** (Boussinesq 1868; the classic laminar-duct series solution) — genuinely 3D profile u(y,z); plus exact 3D Couette | exact |
| Cavity B3 | **Cubic-cavity tables**: Albensoeder & Kuhlmann, *Accurate three-dimensional lid-driven cavity flow*, JCP 206 (2005) — Re 100/400/1000; the 2D midplane must also match Ghia | published tables |
| Cylinder B6 | **DFG/FeatFlow 3D cylinder benchmark** (Turek group; Re = 20, published 3D Cd/Δp — fetch the exact dataset from the FeatFlow benchmark page); **Barkley & Henderson**, *3D Floquet stability of the cylinder wake*, JFM 322 (1996): mode-A onset Re ≈ 188.5, mode-B ≈ 259 — quantitative 3D stability anchors | published + exact-ish |
| Step B9 | Spanwise-periodic recovery of x1/H ≈ 6.26; finite-aspect-ratio end-wall trends vs **Armaly et al.** (expt, 1983) and published 3D-step DNS (fetch; band anchor) | published (band) |
| Blasius B10 | **Yawed (swept) flat plate**: exact 3D similarity — both boundary-layer components follow f′(η), closed-form skin friction 0.664/√Re_x per component | exact |
| de Vahl Davis B7 | **Fusegi, Hyun, Kuwahara, Farouk** (1991): differentially heated *cubical* cavity — 3D Nu_mean / u_max tables at Ra = 10³–10⁶ | published tables |
| Stokes sphere B5 | already intrinsically 3D; extend to the unsteady wake: **Tomboulides & Orszag** (2000) and **Johnson & Patel** (1999) sphere-wake DNS — Cd(Re), St, wake-structure onset bands | published (band) |
| TGV / MMS / RB onset | already 3D by construction | exact |

New 3D benchmark families (metrics + tiers):

- **B3D-1 Rectangular-duct Poiseuille**: exact double-series u(y,z) (truncation < 1e-6 at
  ~200 terms) on a 3D duct (Wall × 4, FixedPressure in/out).  Metrics: profile L2 vs dx
  sweep (spatial-order re-check in 3D); Q vs the series integral.  Tier: L2 < 2% at the
  finest grid; CI smoke at 24 × 16 × 16.
- **B3D-2 Cubic lid-driven cavity**: Albensoeder–Kuhlmann tables at Re 100/400/1000
  (midplane + spanwise structure, corner vortices).  Tier: midplane centerline-L2 < 5%
  (Re = 100), < 15% (Re = 1000); smoke = v-sign structure on the midplane.
- **B3D-3 DFG/FeatFlow 3D cylinder** (Re = 20, span 0.41 with wall z): published 3D Cd/Δp
  band.  Tier: Cd inside the published band at fine grids; the 2D-recovery check
  (thin-span Periodic z) within 10% of the Schäfer–Turek 2D value.
- **B3D-4 Cylinder wake 3D instabilities**: spanwise-periodic run (span ≥ 4D for mode A);
  amplify small random 3D noise at Re = 150/200/260/300; measure the fastest-growing
  spanwise mode.  Metrics: onset bracket (sign flip of the 3D-mode growth) within ±15% of
  the Barkley–Henderson Re_c (188.5 / 259); fastest-growing wavelength within ±20%
  (3.96D / 0.82D); 3D-wake St ≈ 0.2–0.21 band.  Tier: bracket-only at coarse grids;
  wavelength at fine; smoke at Re = 200 = a 3D mode appears.
- **B3D-5 Sphere wake** (Re = 200–350): Cd(Re) vs Schiller–Naumann + the DNS tables; the
  steady → planar-symmetric wake onset; St in the unsteady regime vs the DNS bands.
  Tier: Cd within 10% at Re ≤ 300; wake character qualitative.
- **B3D-6 Yawed flat plate**: uniform streamwise U + spanwise W over a plate; exact 3D
  self-similar profile; metrics: per-component Cf vs 0.664/√Re_x (beyond the corner
  transient), δ*/θ ≈ 2.59.  Tier: Cf within 15%; smoke = spanwise component shape-matches
  f′(η) (L1 < 10%).
- **B3D-7 3D backward-facing step**: (a) spanwise-periodic: x1/H recovers 2D (6.26 band)
  on the midplane; (b) finite-span runs at aspect ratios 2/4/8 vs Armaly's end-wall /
  spanwise-variation trends (laminar region, Re ≈ 300–800).  Tier: x1 within 10%
  (periodic span); spanwise trends qualitative.
- **B3D-8 Cubic differentially-heated cavity** (Ra = 10³–10⁶): Nu_mean vs the Fusegi
  tables; span aspect-ratio sweep 1 → 4 must show the 3D Nu approaching the de Vahl Davis
  2D value.  Tier: Nu within 10% (Ra ≤ 10⁵); the AR-sweep monotonic and matching the
  literature trend.

3D mesh/cost note for Section 3: each refinement doubles the cells per axis = ×8 cost per
step PLUS the SOR iteration-count growth — the 3D smoke tier stays at 32³-ish, full at
128³, reference at 256³ for the exact-anchor families; the band-anchor families are
benchmarked at the grid where their own published runs saturate.

---

## 2B. Full-engine coverage — every physics module and every multiphysics path gets an anchor

**Coverage criteria (what "good coverage" means here):**
1. Every physics module in the engine has at least one benchmark family anchored to an
   exact solution or a public reference (exact preferred; single-group datasets become
   band anchors).
2. Every registered multiphysics coupling path has at least one composite or
   consistency/conservation benchmark (partitioned result vs an analytic composite or an
   exact steady state).
3. Every family has a CI-smoke tier (seconds) and a full tier (sweeps).
4. The matrix below is the living coverage ledger: each family carries a status
   (planned → implemented → verified) and rolls up into the demo-list map.

| Module (dir) | Family | Anchor | Type |
|---|---|---|---|
| fdm3 CFD | B1–B11, B3D-1–8 | as listed in Sections 2/2A | exact/public |
| fdm 2D (legacy) | shares B1/B2/B4 anchors where applicable | — | exact |
| BEM aero (reduced_order/bem) | C1 turbine BEM | Betz 16/27, Glauert a = 1/3 (exact); NREL Phase VI measured power (band) | exact + band |
| fluid/forces | inside C1 (blade element, momentum balance, pressure integration) | force/moment integral checks vs the analytic loads | exact |
| lumped plenum | C2 plenum/stages | isothermal exponential relaxation (exact); Greitzer B-criterion | exact + semi-analytic |
| turbomachinery | C3 compression system | Greitzer surge boundary (analytic); stage-map conservation | semi-analytic + exact |
| electromagnetics (fdtd) | C4 FDTD | Fresnel reflect/transmit 1D (exact); Mie series 2D cylinder (exact); cavity modes (closed form); CFL bound | exact |
| electromagnetics (circuit) | C5 circuits | RC/RLC closed forms; AC phasor steady state | exact |
| electromagnetics (static) | C6 static fields | dipole/plate potentials closed form; Laplace series | exact |
| acoustics | C7 wave solver | d'Alembert + impedance reflection (exact); box eigenmodes (closed form) | exact |
| structural | C8 elasticity/waves | cantilever PL³/3EI + mode β_nL (closed form); c = √(E/ρ) (exact) | exact |
| thermal (CHT) | C9 conjugate HT | two-layer series resistance; fin cosh solution + efficiency; transient slab (Fourier) and semi-infinite (erfc) | exact |
| thermo | C10 EoS/polytropic/steam | polytropic relations (exact); IAPWS-IF97 steam spot values (public table) | exact + public |
| species | C11 species transport | exact decay; Gaussian release σ² = 2Dt; conservation | exact |
| reaction | C12 reactor | CSTR algebraic steady state; PFR exponential profile | exact |
| porous | C13 porous | 1D Darcy linear profile; Forchheimer closed form | exact |
| rigid_body / multibody | C14 rigid dynamics | slider-crank closed forms; pendulum (small-angle exact, elliptic tables); assembly inertia | exact + tables |
| robotics | C15 manipulator | FK vs closed-form DH; trajectory vs analytic profiles | exact |
| particles | C16 particle track | Stokes settling v_t = Δρgd²/18μ; paths in exact flow fields | exact |
| control | C17 controller | PI on 1st-order plant vs closed-loop closed form | exact |
| coupling (manager) | C18 partitioned coupling | the same coupled problem vs an analytic composite; energy/mass conservation across exchanges; fixed-point contraction | exact |
| coupling (turbine system) | C19 coupled turbine | BEM+generator+drive steady operating point vs the analytic power balance | exact |
| EM+thermal+fluid chain | C20 composite multiphysics | 1D Beer–Lambert absorption → Joule heat → steady conduction (quadratic profile) | exact composite |
| species+fluid | C21 composite transport | release in uniform fdm3 flow vs the advected Gaussian (σ² = 2Dt) | exact composite |

Family definitions (metrics + tiers; all get a CI smoke tier):

- **C1 BEM turbine**: Betz C_P = 16/27 at a = 1/3 (smoke, exact); Glauert optimum; NREL
  Phase VI power curve band (full; the definitive public BEM dataset); the water-turbine
  duct augmentation sanity.  Tier: Betz within 1%; Phase VI within the measured spread.
- **C2 Plenum/stages**: isothermal plenum pressure relaxation vs exp(−t/τ) exact; stage
  stack conservation (mass/energy balance < 1%).  Tier: < 1% at 3 time constants.
- **C3 Compression system**: Greitzer B-parameter surge boundary (analytic criterion —
  on/off surge prediction vs the criterion); operating-map consistency.  Tier: boundary
  within ±10%.
- **C4 FDTD**: 1D plane wave at a dielectric interface vs Fresnel coefficients (exact
  R,T); 2D TMz cylinder scattering vs the Mie series coefficients (public); cavity
  eigenfrequencies (closed form); CFL-bound adherence.  Tier: Fresnel within 1% at 20
  cells/λ; Mie within 5% (band); smoke = reflection sign and magnitude.
- **C5 Circuits**: RC step (exact e^{−t/RC}); RLC ringing (closed-form damping/ω_d);
  AC steady state vs phasors.  Tier: < 1% at 3 τ.
- **C6 Static fields**: point-dipole potential vs the closed form; parallel-plate field;
  Laplace box vs the series solution.  Tier: < 2% L2.
- **C7 Acoustics**: 1D pulse propagation vs d'Alembert (exact shift); partial reflection
  at an impedance step (exact amplitudes); box-mode frequencies.  Tier: < 1% (modes),
  < 5% (reflection).
- **C8 Structural**: cantilever tip deflection vs PL³/3EI and mode resonant freqs
  (β₁L = 1.8751…); axial wave speed √(E/ρ).  Tier: < 2%.
- **C9 CHT**: two-layer conduction (exact series resistance); straight fin vs the cosh
  solution and published fin efficiency; transient slab vs the Fourier series and the
  semi-infinite erfc regime.  Tier: < 2% steady, < 5% transient.
- **C10 Thermo**: polytropic pV^n relations exact; steam spot-values vs IAPWS-IF97
  (saturation + superheated table points, < 0.1% band); ideal-vs-real EoS consistency.
- **C11 Species**: exact decay; Gaussian release σ² = 2Dt (L2 < 3%); global conservation
  < 1%.
- **C12 Reaction**: CSTR steady state (algebraic, < 1%); PFR exponential profile
  (< 2%); Arrhenius temperature sensitivity vs the exact ODE solution.
- **C13 Porous**: 1D Darcy u = −(k/μ)∇p linear profile (< 1%); Forchheimer closed form
  (< 3%).
- **C14 Rigid dynamics**: slider-crank position/velocity/acceleration closed forms
  (< 1% at 3 crank speeds); pendulum period vs small-angle (exact) and elliptic-integral
  tables; assembly inertia (exact).
- **C15 Robotics**: forward kinematics vs closed-form DH trig identities (< 1e-6
  relative — pure algebra); trajectory vs analytic profiles (< 1%).
- **C16 Particles**: Stokes settling v_t (< 1%); path in uniform/shear flows vs closed
  forms (< 2%).
- **C17 Control**: PI on a 1st-order plant vs the closed-form closed-loop step response
  (< 1%); gain/phase margins vs the analytic values.
- **C18 Partitioned coupling**: same composite problem solved partitioned vs monolithic/
  analytic steady state (consistency < 2%); energy/mass conservation across the exchange
  boundaries (< 1%); fixed-point convergence rate vs the analytic contraction factor
  (qualitative).
- **C19 Coupled turbine**: the matched-speed operating point: BEM torque = generator
  torque = drive losses at the analytic balance (< 2% power imbalance); speed set-point
  tracking (exact 1st-order response).
- **C20 EM+thermal+fluid composite**: 1D slab: I(x) = I₀e^{−αx}, absorbed power →
  steady conduction with the quadratic interior temperature profile (exact steady state;
  < 2%); the fluid side carries the heat out (B11 channel Nu as the boundary condition
  — links to the verified B-family).
- **C21 Species+fluid composite**: scalar release advected by the uniform fdm3 flow:
  center shift = Ut (exact), spread σ² = 2Dt (< 3%); total mass conserved (< 1%).

Explicitly OUT of scope for the anchor suite (documented, not silently skipped):
compliant-wall FSI (the engine's drag-FSI path cannot do structural-deformation
coupling — the Turek–Hron FSI family needs it), combustion/turbulence DNS-grade
statistics, and the 12-digit CFD excercise beyond the Schäfer–Turek level.

---

## 2C. Delivery record — runnable demo status (benchmarks/benchmark_suite)

Status: FIRST WAVE implemented and smoke-verified (all numbers below from the
smoke tier on the dev machine; full sweeps available via `--full`).  One
runnable binary, `benchmarks/benchmark_suite`, cases dispatched by name
(`--list`, `--case NAME`, `--full`, `--grid N`); CMake target built with the
tests.  SECOND WAVE added with it (all smoke-verified, same suite binary):
B3 cavity, B5 Stokes sphere, B6 Schäfer–Turek cylinder, C4 FDTD, C6 static
fields, C9 thermal.  Remaining (deferred): B7/B8/B11 coupled thermal (needs
the T-field coupling wiring), C18–C21 composites, the B3D series, and the
Re=100 shedding branch of B6 (full tier).

Verified smoke results (metric, measured, tier target):

| Case | Metric | Smoke result | Tier |
|---|---|---|---|
| B1 MMS | window-L2 spatial order (Central) | 1.94 (12/24³) | ~2 |
| B2 TGV | decay-rate ratio (1.0 exact) | 1.050 @ 16³ → 0.993 @ 32³ | within 10% |
| B4 channel | profile L2 / Q rel err | 3.9e-3 / 1.5e-5 | Q < 2% |
| C2 plenum | linearized-Greitzer 2×2 rel err | 0.28% | < 1% |
| C5 circuit | i(4τ) rel err | 3.7e-5 | < 1% |
| C7 wave | f = c/2Lx rel err | 2.6e-4 | < 1% |
| C8 oedometer column | u_z max rel err | 1.2e-3 (M-modulus anchor) | < 2% |
| C10 polytrope | τ closed form + η roundtrip | 3.5e-16 | < 1e-9 |
| C11 species | exponential decay | 8.6e-15 | < 1% |
| C12 reactor | A→B + conservation | 2.6e-11 / 2.2e-16 | < 2% |
| C13 Darcy | steady linear profile L2 | 2.7e-16 | < 1% |
| C14 crank | closed-form kinematics | 1.4e-16 | < 1e-6 |
| C15 FK | 2-link trig chain | 0 | < 1e-9 |
| C16 settling | exact v,z | 7.2e-14 | < 2% |
| C17 PI | closed-loop step response | 2.2e-5 | < 1% |
| B3 cavity | centerline tables vs Ghia Re=100 | L2(u) 0.21; L2(v) coarse-tilt | full: < 5% |
| B5 sphere | Stokes/Oseen (freeze-reaction drag) | 2.4x over at a = 4 cells | full: < 5% at a >= 8 |
| B6 cylinder | Schäfer-Turek Re=20 | Cd 2.20 vs 5.58; dP 0.165 vs 0.118 | full tier |
| C4 FDTD | exact CFL bound + travel time | bound OK; travel 3.9%; drift 0.94% | < 1% |
| C6 static | exact linear plate bridge | 2.2e-8 | < 2% |
| C9 thermal | exact quadratic + Fourier series | 9.2e-6 / 8.2e-5 | < 2% / < 5% |

Recorded findings during the bring-up (all written into the case comments and
the conformance log):
1. **FDM3Solver field()-sync bug FIXED (W18 code)**: `field()` was an
   extract-only adapter — edits (immersed freeze/penalty, hand-set ICs) never
   reached the solver's grid, so the W16 external-manipulation seam was inert
   and its tests passed tautologically.  `step()` now ingests the adapter
   first; the W16 anchors still pass (109/109) with the REAL freeze physics.
2. **TGV**: the (1,1,1) mode decays with k² = 3 (e^(−2νk²t)); report the
   rate ratio, which converges quadratically (1.050/0.981/0.993 at 16/24/32³).
3. **MMS**: the continuous-residual residual drifts linearly in time (operator
   mismatch); the fixed-window L2 is the order metric (1.94 ≈ 2) and the
   drift is reported separately.  Two traps documented: Hybrid→upwind branch
   at cell Pe > 2 adds ~13ν numerical diffusion; the explicit-diffusion
   stability bound is ν·dt/h² ≤ 1/6 (unclamped by the CFL logic).
4. **C8**: the free-lateral body-force column at ν = 0.3 deviates O(1) from
   the uniaxial formula (1.7× at 11×25, ~13× at 7×24) — a limitation of the
   discrete free-surface (Robin) treatment under Poisson contraction; the
   oedometer-confined anchor (M-modulus) is exact to 0.12% and the ν = 0
   free column to 6e-11.  The full-plane-pin mask also stagnates the SOR
   (residual stuck at 1.0); the roller+corner-pins convention converges.
5. **B5/B6 (wave 2)**: the wall-plane differential momentum balance is
   contaminated by the channel's own pressure drop at low Re (10x at the
   smoke horizon); the correct low-Re drag measure is the kinematic-freeze
   REACTION force (sum of rho*V*(u_post - u_pre)/dt over the mask,
   steady-averaged).  At the coarse masks (a = 4 cells, D = 6.4 cells) both
   bodies show a consistent ~2.4x over-prediction - the collocated frozen
   shell's effective-radius blur; the fine-grid convergence is the full-tier
   gate (recorded, not hidden).
6. **B3 (wave 2)**: at 40^2 the u-centerline matches the Ghia tables well
   (core height exact, mid values ~3%) while the v-centerline shows the
   classic coarse-grid tilt (the Ghia 129^2 reference needs 80^2+, the full
   tier).
7. **B6 dP is a direct field measurement** and lands 40% high at the coarse
   mask (0.165 vs 0.118) - the mask-blurred pressure distortion; the same
   fine-grid gate applies.
8. **FDTD**: the uniform-material v1 cannot do the Fresnel/Mie interface
   (recorded); the exact CFL bound is enforced by the config validation
   (courant > 1 rejected), the pulse travel time is reproduced to 3.9%
   with the source-delay-corrected reference, and the PEC-box energy drift
   is 0.94%.

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
- **Phase 6 — module-wide exact anchors (C1–C17)**: cheapest-exact first — C13 Darcy,
  C5 circuits, C11 species, C10 polytrope, C8 cantilever, C7 acoustics, C14 crank,
  C15 FK, C16 settling, C12 reactor, C17 PI, C2 plenum, C1 Betz-level BEM, C6 static,
  C3 Greitzer, C4 Fresnel/Mie (the compute-heavy FDTD cases last in the phase).
- **Phase 7 — multiphysics composites (C9, C18–C21)**: C9 CHT fin (exact, fast),
  C20 EM+thermal+fluid (the demo-list composite), C21 species+fluid, C19 turbine
  operating point, C18 partitioned-coupling consistency/conservation.  Each composite
  family reuses already-verified single-physics anchors as its boundary conditions —
  the multiphysics suite is only as believable as its parts.

Each phase: implement → verify anchors → record rows → commit with the results summary.
Phase 1 is the recommended starting point (fastest complete accuracy story + the perf
harness skeleton).  Coverage criterion: the Section 2B matrix is complete when every
module row has status = verified.
