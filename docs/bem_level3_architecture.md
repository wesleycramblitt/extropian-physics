# Level-3 Duct/Hull-Coupled BEM Solver — Architecture

Status: implementation contract (2026-08-26, rev 2 — post-review)
Owner: extropian-physics
Namespace head: `exd::physics::fluid::reduced_order::bem`

## 0. Scope

Level-3 turbine solver: given the parametric `exd::geometry::TurbineDefinition`
(hull/duct meridional profile + rotor geometry) and operating conditions, estimate
rotor torque/thrust/power, Cp/Ct, radial loading, and an *engineering-approximate*
axial velocity/pressure field. Not a CFD solve; no Navier–Stokes, no mesh solve,
no turbulence model, no arbitrary-axis support, no transient behavior.

Single rotor per solve. Row-resolution truth table (see §5).

## 1. Coordinate system (user decision: geometry frame directly)

The solver consumes `exd::geometry::TurbineDefinition` verbatim. All math happens
in the machine-local meridional frame the geometry module already uses:

- `z` — rotor axis coordinate (geometry "forward"/nose is +z; flow travels −z)
- `r` — radial coordinate
- `θ` — circumferential (only used for chord/solidity conventions)

Conventions (also documented in geometry's turbine.hpp):
- rotor plane: mid-chord plane of the blade row, `z_r = 0.5·(z_le(0.5) + z_te(0.5))`
- upstream = +z, downstream = −z; `V_inf` is a positive magnitude along the axis
- positive rotation about +z is CCW; turbine sense = +stagger + CCW spin
- lengths in meters, angles in degrees (from `TurbineParam.unit`), rpm → rad/s
  internally (`Ω = rpm·2π/60`)

No origin/axis transform layer. World placement is the geometry module's mesh
transform concern.

## 2. Category taxonomy (user decision)

Rearchitect extropian-physics into domain categories. Top-level dirs under
`include/exd/physics/` and `src/` today:

```
physics/
├── fluid/                  ← NEW (domain category)
│   ├── cfd/                ← future: fvm/, fem/, lbm/
│   └── reduced_order/      ← analytical/engineering solvers
│       ├── bem/            ← THIS SOLVER (first occupant)
│       ├── actuator_disk/  ← future
│       ├── vortex/         ← future
│       └── potential_flow/ ← future
├── thermal/                ← future
├── structural/             ← future
├── electromagnetics/       ← future
└── (shared infrastructure, stays at top level: mesh/, field/, bc/,
    material/, coupling/, solver/)
```

- No `PhysicsDomain::TurbineAerodynamics`: BEM is `PhysicsDomain::FluidFlow`,
  fidelity tier expressed by the tree, not the enum.
- All BEM types in namespace `exd::physics::fluid::reduced_order::bem`.
- Shared infra stays where it is; zero churn.
- README rewritten with the full intended taxonomy (Phase 2).

## 3. File layout

```
include/exd/physics/
├── physics.hpp                          # umbrella; includes BEM headers in
│                                        #   dependency order: bem_config.hpp,
│                                        #   airfoil.hpp, bem_result.hpp,
│                                        #   bem_solver.hpp
└── fluid/reduced_order/bem/
    ├── bem_config.hpp                   # BEMSolverConfig, OperatingConditions
    ├── airfoil.hpp                      # AirfoilPolar, AirfoilAssignment,
    │                                    #   PolarDatabase (member evaluate)
    ├── bem_result.hpp                   # TurbineResult + sub-results
    └── bem_solver.hpp                   # solve_turbine() entry point

src/fluid/reduced_order/bem/
    ├── bem_solver.cpp                   # orchestration, element loop, integration
    ├── blade_geometry.cpp               # TurbineDefinition → blade stations
    ├── induction.cpp                    # InductionModel (low/high induction)
    ├── losses.cpp                       # Prandtl tip/hub loss (LossModel hook)
    ├── polar_db.cpp                     # builtin tables + CSV loading
    ├── duct_flow.cpp                    # K_duct acceleration model
    ├── hull_drag.cpp                    # hull forces (CD split)
    └── flow_field.cpp                   # Gaussian wake + Bernoulli pressure field

data/airfoils/{naca0012,naca4412}_{re}.csv   # sample polars (synthetic, labeled)
tests/unit/bem/*.cpp                     # doctest suites (EXT_PHYSICS_BUILD_TESTS)
docs/bem_level3_architecture.md          # this document
```

## 4. Public API (namespace exd::physics::fluid::reduced_order::bem)

### 4.1 airfoil.hpp

```cpp
struct AirfoilPolar {
    std::string name;                 // e.g. "naca0012"
    double re = 0.0;                  // design Re; 0 = Re-independent fallback
    std::vector<double> alpha_deg;    // strictly increasing
    std::vector<double> cl;
    std::vector<double> cd;
    [[nodiscard]] bool valid() const;
    [[nodiscard]] std::pair<double, double> evaluate(double alpha_deg) const;
    // linear interp; clamp flat beyond table ends
};

struct AirfoilAssignment {            // span (0=hub,1=shroud) → airfoil id
    double span = 0.0;
    std::string airfoil;
};

class PolarDatabase {
public:
    void add(AirfoilPolar p);                              // replace same name+re
    void add_builtin_polars();                             // synthetic NACA0012/4412
    bool load_csv(const std::string& path);                // cols: name,re,alpha_deg,cl,cd
    bool load_directory(const std::string& path);          // all *.csv
    const AirfoilPolar* find(const std::string& name, double re) const;
    bool has(const std::string& name) const;
    std::vector<std::string> airfoil_names() const;
    bool empty() const;
    explicit operator bool() const { return !empty(); }
};
```

Polar Re-selection rule (single, unambiguous): among polars with the requested
`name`, pick the one minimizing |re_target − re|; a polar with `re == 0`
(Re-independent) has distance +inf and is used only when no re-tagged polar
exists for that name at all. `find` returns `nullptr` for unknown names.
(Callers must null-check; the solver converts a missing airfoil into a clear
invalid-result error, never a crash.)
`add_builtin_polars()` tables are the single source of truth; the vendored
CSVs in `data/airfoils/` are exported copies of the same data (labeled
"synthetic — replace with XFoil data").

### 4.2 bem_config.hpp

```cpp
struct OperatingConditions {
    double v_inf = 10.0;                // m/s, > 0 (validation)
    double rho   = 1.225;               // kg/m³, air default
    double mu    = 1.81e-5;             // Pa·s  (Re = ρ·W·c/μ)
    double p_ref = 101325.0;            // Pa, freestream static pressure
    std::optional<double> rpm_override; // else BladeRow::rotational_speed
};

enum class ReferenceArea { RotorDisk, Annulus };   // Cp/Ct denominator

struct BEMSolverConfig {
    uint32_t element_count = 32;        // 0 disallowed; ≥ 4 validated; warning if
                                        //   not in {16,32,64,128} (spec's sweep set)
    double   k_duct = 0.5;              // [0,1] duct acceleration coefficient
    double   under_relaxation = 0.25;   // (0,1]
    double   induction_tolerance = 1e-5;
    uint32_t max_iterations = 100;
    double   glauert_threshold = 0.4;   // Buhl branch guard (see §6.2)
    uint32_t row_index = 0;
    double   hull_cd = 0.2;             // total hull drag C_D (frontal ref area)
    ReferenceArea reference_area = ReferenceArea::RotorDisk;  // A=πR_tip² per spec;
                                    // Annulus: π(R_tip²−R_hub²) — both documented
    uint32_t field_axial_points = 64;
    uint32_t field_radial_points = 16;
    double   upstream_extent = 0.0;     // 0 → 3·R_tip
    double   wake_length     = 0.0;     // 0 → 8·R_tip
    double   wake_decay_length = 0.0;   // 0 → 4·R_tip
    double   wake_expansion  = 0.05;    // k_w, dR_w/dz
    double   wake_radius_initial = 0.0; // 0 → R_tip
    std::vector<AirfoilAssignment> airfoils;  // empty → all "naca0012"
    bool     include_flow_field = true; // cost is negligible (~µs at 64×16)
};
```

Validation errors (element_count < 4 or 0-sentinels where invalid, k_duct∉[0,1],
under_relaxation∉(0,1], glauert_threshold∉(0,1), row_index OOB, negative
lengths) → invalid result with message.

### 4.3 bem_result.hpp

```cpp
struct RadialStation {
    double radius_m, span, axial_velocity, relative_velocity,
           induction_axial, induction_tangential,
           inflow_angle_deg, angle_of_attack_deg, reynolds,
           cl, cd, lift_per_m, drag_per_m,
           torque_per_m, thrust_per_m,   // ring totals (B blades)
           pressure_jump;
    bool   converged;
    uint32_t iterations;
};

struct RotorResults {
    double torque, thrust, power;       // SI
    double cp, ct;                      // ref area per ReferenceArea
    double efficiency;                  // = Cp / (16/27) (Betz-relative)
};

struct HullForces { double drag, pressure_drag, viscous_drag, cd, reference_area; };
struct DuctState  { double k_duct, area_ratio, m_duct, v_rotor; };

struct SystemResults {
    double net_thrust;      // T_rotor − D_hull
    double net_power;       // = rotor power (hull drag is a support load, not shaft loss)
    double efficiency;      // = rotor.efficiency · (T_net/T_rotor) — engineering FOM,
                            //   documented as such in the header
};

struct FlowFieldGrid {
    std::vector<double> z;              // axial stations, ascending (upstream high z)
    std::vector<double> r;              // radial stations, 0 → R_max;
                                        //   R_max = max(R_shroud(z)) over the grid
                                        //   domain (≥ R_tip)
    std::vector<double> velocity;       // axial flow magnitude, row-major [iz*nr+ir]
    std::vector<double> pressure;       // gage vs p_ref, same layout
};

struct TurbineResult {
    bool valid = false;
    std::string error;
    bool converged = false;
    std::vector<std::string> warnings;
    RotorResults rotor;
    HullForces hull;
    DuctState duct;
    SystemResults system;
    std::vector<RadialStation> radial;
    FlowFieldGrid flow_field;
};
```

### 4.4 bem_solver.hpp

```cpp
TurbineResult solve_turbine(const exd::geometry::TurbineDefinition& turbine,
                            const OperatingConditions& conditions,
                            const PolarDatabase& polars,
                            const BEMSolverConfig& config = {});
```

## 5. Geometry ingestion (blade_geometry.cpp)

Mimics exd-geometry's turbine.cpp semantics exactly:

1. Rotor-row resolution (truth table):
   - 0 Rotor rows → invalid result, "no Rotor row in TurbineDefinition".
   - 1 Rotor row → solve it. If `config.row_index` does not match that row's
     index → warning (row_index ignored).
   - >1 Rotor rows → invalid result, "N Rotor rows; single-rotor solver".
   Non-rotor rows are always ignored (warning emitted once).
2. Per span fraction f∈[0,1], straight lines (as geometry does):
   `le(f) = (1−f)·LE_hub + f·LE_shroud`, `te(f)` likewise (Vec2f in (z,r)).
3. Rotor plane `z_r = 0.5·(le_x(0.5) + te_x(0.5))`.
4. `chord(f) = |te(f) − le(f)|` — Euclidean in the (z,r) meridional plane,
   identical rule to turbine.cpp. Aerodynamic chord is approximated by this
   meridional 2-D magnitude (documented; slightly overestimates
   perpendicular-to-span chord for swept blades). `BladeSection::chord` is
   NOT used (matches geometry).
5. Stagger β(f) in deg: linear interpolation over `sections` by span. If
   `sections` is empty → mirror geometry synthesis: stagger = 0 (zero-camber
   flat plate) + warning "BladeRow.sections empty; zero-stagger defaults".
   Out-of-range spans → clamp + warning if extrapolated.
6. Airfoil assignment: `config.airfoils` bracket-selected (nearest assignment
   with span ≤ element span; none → take lowest-span assignment; empty list →
   all "naca0012").
7. Radial stations: N midpoints `r_i = R_hub + (i+0.5)·dr`,
   `dr = (R_tip−R_hub)/N`, with `R_hub = 0.5·(le_y(0)+te_y(0))`,
   `R_tip = 0.5·(le_y(1)+te_y(1))`.
   Tip-clearance clamp: if `row.tip_feature == TipFeature::Clearance`, clamp
   `R_tip ≤ R_shroud(z_r) − tip_clearance` (mirrors the mesh generator; no
   clamp otherwise). Validate R_tip > R_hub.
8. Shroud/hub splines: `MonotoneCubicSpline` from `flow_path.shroud_points` /
   `hub_points` (public exd::geometry type; float-precision, upcast to double —
   documented; guards in the duct model use a small epsilon margin).
9. `HubDefinition` is intentionally ignored: the BEM treats the hub as a body
   of revolution defined by the flow-path hub spline only.
10. Errors: no rotor row (per table), rpm ≤ 0, chord ≤ 0, R_tip ≤ R_hub,
    empty flow path, element_count < 4.

Element inputs: {r, dr, chord, β_deg, span f, airfoil id}.

## 6. Flow model

### 6.1 Duct acceleration (duct_flow.cpp) — labeled engineering estimate

- `A(z) = π·R_shroud(z)²` (outer wall), `A_r = A(z_r)`, upstream reference
  `A_u = A(z_u)` with `z_u` = front end of the shroud spline (or
  `inlet_station` front z if provided; straight/absent shroud → M = 1).
- `M_duct = 1 + K_duct·(A_u/A_r − 1)`; error if `M_duct ≤ 0.05`.
- `V_rotor = M_duct·V_inf` (uniform at the rotor plane).

### 6.2 Induction (induction.cpp)

`InductionModel` with `solveLowInduction()` / `solveHighInduction()` (spec §14;
isolated so the correction is replaceable later — e.g. iterative Glauert or
Pitt–Peters).

Per element, per iteration:

```
V_a = V_rotor·(1−a)      V_t = Ω·r·(1+a′)
W = √(V_a² + V_t²)       φ = atan2(V_a, V_t)
α = φ − β
Cl,Cd = polar.evaluate(α);  Re = ρ·W·c/μ (polar re-selected per iteration:
        find(name, Re) each time — nearest-Re polar, see §4.1)
C_n = Cl·cosφ + Cd·sinφ   C_t = Cl·sinφ − Cd·cosφ   (spec §12; Hansen-consistent)
F = F_tip·F_hub (Prandtl), clamped ≥ 1e-3
σ = B·c/(2πr)
```

Axial — low induction (branch chosen from the PREVIOUS iterate's `a` for
stability):

```
use_buhl = (a_prev > glauert_threshold) && (rad_discriminant >= 0)
rad = C_T·(50 − 36F) + 12F·(3F − 4),  C_T = σ·C_n·(1−a)²/(F·sin²φ)
a_calc = use_buhl
       ? (18F − 20 − 3·√rad) / (36F − 50)      // Buhl closed form
       : σ·C_n / (4·F·sin²φ + σ·C_n)           // momentum closed form
if (a_prev > glauert_threshold && rad < 0) warning("Buhl discriminant < 0; "
    "low-induction branch used")
```

Glauert threshold: default 0.4, per user spec. This is the safe default for
the *algebraic* Buhl closed form (at F=1 the discriminant is positive only for
C_T > 0.857, i.e. a > 0.311; at a = 0.2 the algebraic form is undefined —
Buhl's empirical NREL fit uses 0.2 with a different model). Document this in
the header; value is configurable.

Tangential: `a′_calc = σ·C_t / (4·F·sinφ·cosφ − σ·C_t)` with an explicit guard:

```
denom = 4·F·sinφ·cosφ − σ·C_t
if (denom <= 1e-12) { a′ = (C_t > 0) ? 1.0 : 0.0; warning("tangential denom ≤ 0"); }
else a′ = clamp(σ·C_t/denom, 0.0, 1.0);   // clamp + warning if engaged
```

Same for `a`: clamp to [0, 1) with warning if clamped (C_n < 0 deep stall).

Under-relax λ (default 0.25) both factors. Relative convergence:
`|Δa| < tol·max(a, 1e-2)` and same for a′ (tol = induction_tolerance; absolute
1e-5 floor behavior preserved near a→0). Max 100 iterations; element not
converged → result.converged = false, per-station flag, warning.

Degenerate station guard: where `sinφ < 1e-3` (hub of very high-TSR rotors),
set a = a′ = 0, converged = true, iterations = 0 (station contributes no
induction; documented).

### 6.3 Losses (losses.cpp)

Prandtl tip/hub (spec §11): `f_tip = (B/2)·(R−r)/(r·sinφ)`,
`f_hub = (B/2)·(r−R_hub)/(R_hub·sinφ)`, `F_tip|F_hub = (2/π)·acos(exp(−f))`,
`F = F_tip·F_hub`, clamp F ≥ 1e-3. Hosted behind a `LossModel` hook so Snel /
Du–Selig corrections can be added later. Behavior at φ→0 documented (§6.2).

## 7. Integration

```
Q = Σ B·r·(dL·sinφ − dD·cosφ)          T = Σ B·(dL·cosφ + dD·sinφ)
P = Q·Ω
λ_tsr = Ω·R_tip/V_inf
Cp = P / (½ρ·A_ref·V_inf³)   Ct = T / (½ρ·A_ref·V_inf²)
A_ref = π·R_tip²             (RotorDisk, default; per user spec)
      | π·(R_tip²−R_hub²)    (Annulus; documented non-standard convention:
                               standard for a hub-occupied disk, but references
                               usually quote πR² for open rotors)
```

## 8. Hull forces (hull_drag.cpp) — engineering estimate

- `D_hull = ½ρV_inf²·C_D·A_f`, `A_f = π·R_shroud_max²`, C_D = config.hull_cd.
- Viscous split: turbulent flat-plate `C_f = 0.074·Re_L^(−0.2)`,
  `Re_L = ρV_inf·L/μ`, L = shroud axial length; wetted area
  `A_wet = 2π·∫ r·√(1+r′²) dz` over shroud spline.
- `D_viscous = ½ρV_inf²·C_f·A_wet`; `D_pressure = D_hull − D_viscous`
  (clamped ≥ 0, warning if clamp engaged).
- `T_net = T_rotor − D_hull`.

## 9. Flow field (flow_field.cpp) — explicitly engineering estimate

Grid: `z ∈ [z_r + upstream_extent, z_r − wake_length]` (ascending, upstream
high z), `r ∈ [0, R_max]` (R_max = max shroud radius over the grid domain).
Interpolate element induction `a(r)` across stations; `a(r) = 0` outside
[R_hub, R_tip] (no double-count on the outboard region).

Velocity (axial magnitude):
- upstream (z > z_r): `V(z) = V_inf·(1 + K_duct·(A(z)/A_r − 1))` inside the
  shroud (r ≤ R_shroud(z)); `V_inf` outside the duct or beyond the duct ends.
  At z = z_r⁺ this equals `V_rotor` by construction — the duct model is
  continuous with the disk.
- disk (z = z_r): `V = V_rotor·(1 − a(r))` — continuity with the BEM.
- downstream (z ≤ z_r): deficit relative to `V_rotor`, growing toward the
  classical far-wake value `2a` and spreading radially:

```
V(z,r) = V_rotor·(1 − a(r)·[1 + (1 − g(z))·exp(−(r/R_w(z))²)]),   z ≤ z_r
g(z) = exp(−(z_r − z)/L_w)      L_w = wake_decay_length (default 4·R_tip)
R_w(z) = R_w0 + k_w·(z_r − z)   R_w0 = wake_radius_initial (default R_tip)
```

  Check anchors:
  • z = z_r (g = 1): V = V_rotor·(1 − a(r)) EXACTLY at every radius —
    the disk plane matches the BEM velocity field by construction
    (the radial Gaussian modulates only the downstream development).
  • z → −∞ centerline (r = 0): V → V_rotor·(1 − 2a) ✓ classical far wake.
  • r ≫ R_w: deficit → a(r) (no false recovery inside streamtubes).
  • no overshoot: deficit factor ∈ [a, 2a] everywhere downstream.
  No target-plane closure: the far wake does not return to V_inf within the
  grid — single-streamtube momentum result, documented.
- solid hub body: `r < r_hub(z)` (hub spline; 0 if absent) → V = 0.
- V clamped ≥ 0.01·V_inf to keep pressure finite.

Pressure:
- Bernoulli everywhere: `p = p_ref + ½ρ(V_inf² − V²)`.
- Pressure drop across the rotor at the disk plane only:
  `Δp(r) = dT(r)/(2πr·dr)` (per element), applied for z ≤ z_r:
  `p(z,r) −= Δp(r)`.
- Explicit caveat (in code + README): total pressure is not conserved
  downstream and the field does not recover to freestream; this is an
  analytical/engineering pressure estimate, NOT a CFD pressure solution.

## 10. CMake

- FetchContent `exd-geometry` (mirror exd-core pattern from
  extropian-geometry/CMakeLists.txt, which itself supports `EXD_CORE_DIR`;
  plumb the same override here):
  ```
  if(DEFINED EXD_GEOMETRY_DIR)
      FetchContent_Declare(exd-geometry SOURCE_DIR "${EXD_GEOMETRY_DIR}")
  else()
      FetchContent_Declare(exd-geometry
          GIT_REPOSITORY https://github.com/wesleycramblitt/extropian-geometry.git
          GIT_TAG main)
  endif()
  if(DEFINED EXD_CORE_DIR)   # mirrored from exd-geometry's own pattern
      FetchContent_Declare(exd-core SOURCE_DIR "${EXD_CORE_DIR}")
  endif()
  set(ENABLE_TEXT  OFF CACHE BOOL "" FORCE)   # geometry default needs FreeType
  set(BUILD_TESTS  OFF CACHE BOOL "" FORCE)   # /HarfBuzz; force BEFORE
  set(BUILD_DEMO   OFF CACHE BOOL "" FORCE)   # FetchContent_MakeAvailable
  set(EXT_CORE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(exd-geometry)
  ```
- Link `exd::geometry` PUBLIC (its types appear in the public BEM API).
- Add the 8 BEM sources to the `exd-physics` target.
- Tests: (exists) `EXT_PHYSICS_BUILD_TESTS` option; create `tests/CMakeLists.txt`
  with doctest v2.4.11 FetchContent (mirror geometry) and per-file test
  executables; the test target gets
  `target_compile_definitions(<test> PRIVATE EXT_PHYSICS_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")`
  and tests build the airfoil dir path as `EXT_PHYSICS_DATA_DIR + "/airfoils"`
  (no relative paths — ctest runs from the build dir).

## 11. Tests (doctest, tests/unit/bem/)

1. `polar_test.cpp` — add/find nearest-Re; re==0 fallback only when no re-tagged
   polar; `find("unknown", re) == nullptr` and solver converts to a clean error;
   evaluate() linear interp + flat clamp; CSV load (valid + malformed);
   builtins present.
2. `blade_geometry_test.cpp` — annulus TurbineDefinition (geometry-test style):
   R_hub/R_tip/chord/β extraction, z_r; tip-clearance gating on feature flags;
   row-resolution truth table (0/1/many rotors); error paths (chord ≤ 0,
   rpm = 0, element_count < 4); empty-sections warning; HubDefinition ignored.
3. `bem_betz_test.cpp` — synthetic 2π-slope Cd=0 polar, cylindrical shroud:
   positive T/Q/P; Cp slightly below 16/27 (tip/hub loss keeps it under);
   a ∈ (0, 0.5); converged. NOTE: Betz test never exercises the Buhl branch
   (a = 1/3 < 0.4) — do not assert Buhl there. Golden values pinned after
   first run (doctest::Approx).
4. `duct_test.cpp` — K_duct=0 → M=1; K_duct=1 → M=A_u/A_r; converging shroud
   raises Cp vs cylindrical; divergence lowers it; M_duct ≤ 0.05 → error.
5. `induction_test.cpp` — high-loading case (low TSR): Buhl branch engages
   (a > 0.4, discriminant ≥ 0), converged a < 1, no NaN; discriminant < 0
   fallback + warning; tangential denom guard; a/a′ clamps; station converged
   flag/iteration bounds; relative-convergence behavior at light load.
6. `flow_field_test.cpp` — upstream duct model continuity at z_r: V(z_r⁺) =
   V_rotor; disk plane matches BEM exactly at every radius: V(z_r) = V_rotor(1−a(r)); disk pressure
   jump: ΣΔp·dA ≈ rotor thrust; wake: centerline deficit grows from a toward
   2a downstream, radial profile spreads (R_w grows), V → V_rotor far
   radially; axisymmetry (r=0 sample = centerline).
7. `hull_test.cpp` — hull_cd=0 → zero drag; viscous ≤ total; T_net = T − D.
8. `integration_test.cpp` — full solve with builtin polars + vendored CSV DB:
   valid, converged, warnings empty; element_count ∈ {16,32,64,128}: Cp
   converges (|ΔCp| < 1% between 64 and 128).

## 12. Phasing

- Phase 1: CMake + taxonomy dirs + headers + blade_geometry + airfoil/polars
  + BEM core (losses, induction, integration) + tests 1–3, 5.
- Phase 2: duct_flow, hull_drag, flow_field + tests 4, 6–8 + README/docs
  rewrite (taxonomy + BEM usage section).

## 13. Known limitations (documented)

- K_duct duct model, Gaussian wake, Bernoulli pressure: engineering estimates
  (upgrade path: panel/actuator-disck-coupling → Level 4).
- Total pressure not conserved downstream; field does not recover to V_inf.
- Chord = meridional |TE−LE| (mirrors geometry; approximates aero chord).
- Airfoil id lives in solver config, not BladeSection (no field exists in
  exd-geometry; additive `airfoil` field there is a follow-up).
- Single rotor, axial-only; hub body viscous drag folded into hull_cd.
- Synthetic built-in polars must be replaced with XFoil data for real studies.
- Float-precision splines from exd::geometry upcast to double.

## 14. References

- Burton, N. Jenkins, D. Sharpe, E. Bossanyi, "Wind Energy Handbook", 2nd ed.,
  Wiley, 2011 — BEM formulation (ch. 3), Prandtl losses, Glauert correction.
- Moriarty & Hansen, "AeroDyn Theory Manual", NREL/TP-500-36881, 2005 — Buhl
  high-induction closed form (a = (18F−20−3√(CT(50−36F)+12F(3F−4)))/(36F−50)).
- Hansen, "Aerodynamics of Wind Turbines", 3rd ed., Earthscan, 2015 — classical
  induction closed forms a = σCn/(4F sin²φ + σCn), a′ = σCt/(4F sinφcosφ − σCt).
- User specification (Level-3 duct/hull-coupled BEM, 2026-08): formulas
  §6–§19 as cited inline in this document; deviations noted where they occur.

## 15. Correction models (Phase 3)

Phase 3 adds pluggable correction models for the induction and loss computations.
Models are selected via `BEMSolverConfig` enums and dispatched by factory functions.

### 15.1 Induction correction models

Enum: `InductionCorrection { Standard, GlauertIterative, Snel }`

**Standard** (existing, default):
- Buhl closed-form algebraic solution to Glauert's empirical CT curve
  `CT = 4a(1 - 0.25(5-3a)a)` for a > glauert_threshold.
- One-shot, no iteration needed beyond the momentum/BE coupling.

**GlauertIterative**:
- Same empirical CT curve as Buhl, but solved via fixed-point iteration
  with under-relaxation instead of the algebraic closed form.
- For CT_blade < 0.889: `a_new = CT / (4F)`.
- For CT_blade >= 0.889: `a_new = 0.5*(1 + sqrt(1 - CT/F))`.
- Converges to the same answer as Buhl for steady-state; the iterative
  approach allows smooth blending between momentum and empirical regimes.
- Reference: Glauert (1935), per NREL AeroDyn convention.

**Snel**:
- Exponential blending: `a = a_mom * (1 - exp(-4F sin^2(phi) / (sigma Cn) * a_mom))`.
- Smooth transition, no hard threshold.
- For small a_mom: a ≈ a_mom (momentum regime).
- For large a_mom: growth is sub-linear, naturally capping induction.
- Reference: Snel et al. (1993).

### 15.2 Loss correction models

Enum: `LossCorrection { Prandtl, DuSelig, Chaviaropoulos }`

**Prandtl** (existing, default):
- Classic: `F = (2/pi) * acos(exp(-f))` with `f = (B/2)*(R-r)/(r*sin(phi))`.
- Applied multiplicatively to tip and hub.

**DuSelig**:
- Du & Selig (1993) modified tip loss with Lorentzian g-factor.
- For eta = r/R_tip > 0.95: `f = f_prandtl * (1 + g)` where
  `g = (1.386B - 1.964) * (eta-1) / sqrt(1 + ((eta-1)/0.464)^2)`.
- Smoother transition near the tip than plain Prandtl.
- Hub loss: unmodified Prandtl.

**Chaviaropoulos**:
- Chaviaropoulos & Hansen (2000) loading-dependent tip loss.
- `f = f_prandtl * (1 + alpha * sigma)` where alpha = 1.5 and
  sigma is the local solidity `B*c/(2*pi*r)`.
- Loss grows with increasing local blade loading.
- Hub loss: unmodified Prandtl.

### 15.3 Model selection

```cpp
BEMSolverConfig config;
config.induction_correction = InductionCorrection::GlauertIterative;
config.loss_correction = LossCorrection::DuSelig;
TurbineResult result = solve_turbine(turbine, conditions, polars, config);
```

3 × 3 = 9 possible model combinations. All combinations are valid; the solver
does not enforce any pairing constraints.

### 15.4 File layout additions (Phase 3)

```
src/fluid/reduced_order/bem/
    induction_corrections.cpp   # GlauertIterative, Snel induction models
    loss_corrections.cpp        # DuSelig, Chaviaropoulos loss models

tests/unit/bem/
    induction_corrections_test.cpp   # 5 tests for induction models
    loss_corrections_test.cpp        # 6 tests for loss models
```

### 15.5 Tests 9–14 (Phase 3)

9. `induction_corrections_test.cpp` — GlauertIterative engages for high loading;
   GlauertIterative Cp within 5% of Standard; Snel smooth induction; Snel
   induction bounded below Standard; all three models selectable and valid.
10. `loss_corrections_test.cpp` — DuSelig valid with Cp in (0.2, 16/27);
    DuSelig vs Prandtl comparison; Chaviaropoulos valid; loss factor bounded;
    all models produce finite results; all models selectable.
