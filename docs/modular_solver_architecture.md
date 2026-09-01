# Modular Solver System — Architecture & Implementation Plan

Authoritative roadmap for turning extropian-physics into a modular, multi-domain,
multi-discretization solver system. Companion to `agent_guide.md` (conventions)
and `bem_level3_architecture.md` (reference implementation pattern).

---

## 1. Vision

Build physics solvers like lego: independent **domain modules** that plug
together through stable interfaces ("seams"), where each domain may be solved
by one or more **discretization methods**, all orchestrated by a **coupling
framework**.

```
 APPLICATION LAYER     engine · compressor · turbine · generator apps
 MACHINE LAYER         RotatingAssembly — generic rotating machine
 DOMAIN LAYER          mechanics · fluid::forces · thermo · electrical
                       · control · structural · thermal · acoustics ...
 DISCRETIZATION LAYER  fdm (active) · fdtd (planned) · fvm (planned)
                       · fem (planned) · immersed/particle (later)
 DATA/CONTRACT LAYER   fields & field channels · mesh · bc · time stepping
                       · integrators · units · material
```

**Rules that make this work (from agent_guide, applied here):**

- Namespace mirrors directory: `exd::physics::<domain>`.
- Strategy + Factory for pluggable models; concrete models in anonymous
  namespaces in `.cpp`; dispatcher factories on config enums.
- Value-type results, no exceptions, `ModelStatus` (ok/error/warnings).
- Public headers only for cross-module contracts; `*_internal.hpp` in `src/`.
- Coupling stubs from day one: every solver exposes its exchange channels.
- **Extract a seam only when a second variant exists** (repo doctrine): do not
  invent speculative interfaces. E.g. `IFluidSolver` waits for FVM; EM field
  channels generalize the fluid sampler only when EM lands.

---

## 2. Current state (inventory)

### 2.1 Modules

| Module | Namespace | Status |
|---|---|---|
| FDM 2D fluid | `fluid::fdm` | ✅ active (2D **collocated** cell-centered, SIMPLE, RK/Euler/Heun/CN — docs previously claimed staggered; grid is collocated, one ghost layer, SOR 5-point Poisson) |
| BEM turbine (reduced-order) | `fluid::reduced_order::bem` | ✅ active (standalone) |
| Generic force evaluation | `fluid::forces` | ✅ active (pressure/momentum/table evaluators, blade surface) |
| Generic rotating machine | `mechanics` | ✅ active (axis rotation, moment models, assembly) |
| Turbine application | `turbine` | ✅ active (step/simulate over generic stack) |
| Coupling field sampling | `coupling` | ✅ active (`IFlowField3D`, uniform/structured/FDM adapters, `sample_flow`) |
| Mesh types + I/O | `mesh` | ⚠️ types only in src (no public headers); all gmsh/vtk/cgns/stl/obj/openfoam read/write are EMPTY STUBS |
| Fields | `field` | ⚠️ infrastructure (not yet consumed by solvers) |
| BC framework | `bc` | 🔴 headerless stub (types in framework.cpp only, zero consumers, serialization empty) |
| Coupling manager | `coupling` | 🔴 placeholder (`transfer(nullptr,...)`; no real data movement; no public header) |
| Time stepping | `solver` | 🟡 private (`TimeStepper`/`ConvergenceMonitor` live in `.cpp`, not exported) |
| FDM 3D fluid | `fluid::fdm3` | ✅ active (collocated SIMPLE, 7-point SOR, 6-face BCs, body forces, persistent `FDM3Solver`) |
| Engine app | `engine` | ✅ active (Otto Wiebe + steam Rankine-lite, governor, CSV motion output) |
| Output channels | `io` | ✅ active (exd-fld v1 stamps + timeline, CSV series, cadence policy) |
| Turbine app coupled | `turbine::run_coupled_turbine` | ✅ active (actuator-disk in fdm3, BladeElement forces, soak-verified) |
| Steam EOS | `thermo` | ✅ active (Clausius–Clapeyron saturation; swap to IF97 tables later) |
| Turbomachinery | `fluid::turbomachinery` | ✅ W9: axial mean-line stage/stack, operating map, map lookup — §10 G.1 |
| Plenum / 0D volumes | `fluid::lumped` | ✅ W9: Greitzer plenum, generic gas-network node — §10 G.3 |
| Polytropic processes | `thermo` | ✅ W9: stagnation-polytrope primitives; engine volume-polytropes stay local — §10 G.2 |
| Coupling exchange | `coupling` | ✅ W10: real `CouplingManager`/`SurfaceMapper`/`CoupledSimulation` (was placeholder) |
| Thermal (grid-first) | `thermal` | ✅ W10: steady conduction/advection, `IScalarField3D` |
| Acoustics (grid-first) | `acoustics` | ✅ W10: leapfrog wave, `IScalarField3D` |
| Structural (grid-first) | `structural` | ✅ W10: Navier–Cauchy displacement SOR, `IVectorField3D`, thermal-strain coupling |
| Particles (Lagrangian) | `particles` | ✅ W10: tracks over sampled channels |
| Chemistry (0D) | `chemistry` | ✅ W10: mass-action + Arrhenius reactor |

### 2.2 Seams (the lego connectors)

| Seam | Producer | Consumer | Status |
|---|---|---|---|
| `IFlowField3D::sample(p)→{v,p}` | any fluid method | force evaluators, apps, coupling | ✅ |
| `IForceEvaluator::compute(blade, flow, ω, axis)` | fluid::forces | `RotatingAssembly`, apps | ✅ |
| `IMomentModel::moment(ω)` | mechanics (curve/linear/const) | `RotatingAssembly` | ✅ — **natural home for generator/motor/brake/friction models** |
| `IRotationalDynamics::advance(dt, M)` | mechanics | `RotatingAssembly` | ✅ |
| `ISolverPlugin` (initialize/step/get_field/get_coupling_surface) | external solver repos | `SolverManager`, `CouplingManager` | ✅ |
| Field channels (vector/scalar field at points) | any domain (EM, thermal…) | cross-domain exchange | 🔲 generalize `IFlowField3D` when EM lands |
| `IRigidBodyDynamics` (6-DOF) | mechanics | apps, contacts, FSI | 🔲 Phase B |
| `IController` (PI/PID, pitch, governor) | control | turbine/engine apps, compression systems | 🔲 Phase C |
| `IEos` (equation of state) | thermo | engine/compressor/thermal | 🔲 Phase C |
| `IFluidSolver` (method-switchable fluid step) | fdm + fvm | coupled runs | 🔲 Phase J (trigger: FVM exists) |

---

## 3. Phasing principle

Priority per project direction:

1. **Foundations + machine domains + EM with coupling** (multiple time
   integration methods, general rigid bodies, engines, turbines, compressors,
   EMF/EF/MF) — everything here runs on the FDM/grid machinery + generic
   mechanics that already exist.
2. **Rest of the major domains** (thermal transfer, structural/soft bodies,
   acoustics, particles/chemistry).
3. **Discretization breadth** (FVM, FEM, unstructured/boundary-fitted/immersed
   meshes, method-switchable fluid interface, cross-method verification).

Grid-first: FDM (and FDTD, which is structurally identical on grids) carries
fields as far as it can; FVM/FEM arrive when the mesh/ machinery demands them.

---

## 4. Phase A — Shared infrastructure: multiple time integration methods ✅ DONE

Status: `solver::integrate_step` (8 methods: ForwardEuler, BackwardEuler,
Heun, RK4, CrankNicolson, SymplecticEuler, Verlet, AdaptiveRK45) live in
`solver/integrators.*`; `TimeStepper`/`ConvergenceMonitor` promoted to public
`solver/time_stepping.hpp`; canonical `ModelStatus` moved to the umbrella level
(`exd::physics::ModelStatus`, `mechanics::ModelStatus` kept as an alias).

**Design decisions recorded during implementation:**
- FDM keeps its field-optimized integrators (they operate on the staggered
  grid layout directly; packing/unpacking a flat state per step would be pure
  overhead in the CFD hot loop). The shared module is the canonical API for
  system-level ODEs (machines, circuits, control, engine crank).
- AdaptiveRK45 reports the accepted dt via `dt_used`; with relative error
  control the tolerance-limited step is constant for pure exponentials (local
  error scales with |y|) — step *growth* is observed on low-error solutions.
  Both behaviors are pinned by tests.

Everything later (rigid bodies, engine crank dynamics, circuits, control
loops) consumes the same integrators. Today the FDM solver owns its RK4/Heun
implementations privately and `TimeStepper` is not exported.

### Deliverables

1. **Promote time stepping to public API**
   - New public header `include/exd/physics/solver/time_stepping.hpp`
     declaring `TimeStepper` (CFL-adaptive, dt clamps, `adapt()`, `advance()`)
     and `ConvergenceMonitor` (tolerance, residual, windowed checks).
   - Move the implementations from `src/solver/time_stepping.cpp` (keep the
     same file as implementation home; header now public).
   - FDM config gains `use_shared_stepper` default true; `step_fdm` switches
     to the shared stepper (behavior unchanged, `cfl_target` etc. preserved).

2. **Shared integrator module** — `include/exd/physics/solver/integrators.hpp`,
   `src/solver/integrators.cpp`, namespace `exd::physics::solver`:

   ```cpp
   enum class IntegrationMethod : uint8_t {
       ForwardEuler, BackwardEuler, Heun, RK4, CrankNicolson,
       SymplecticEuler,   // Hamiltonian systems: rigid bodies, springs
       Verlet,
       AdaptiveRK45,      // embedded RK45 with error control
   };
   struct IntegratorConfig { IntegrationMethod method; double rel_tol; double abs_tol; };
   // Integrate one step over a differentiable state vector.
   // state is std::span<double>; deriv(state, dstate, t) fills dstate.
   bool integrate_step(const IntegratorConfig&, double t, double dt,
                       std::span<double> state,
                       const std::function<void(std::span<const double>, std::span<double>, double)>& deriv,
                       ModelStatus& status);
   ```

   - Implicit methods (BackwardEuler, CrankNicolson) via fixed-point
     iteration with under-relaxation (Newton-lite later); no exceptions.
   - `AdaptiveRK45` adjusts dt internally, reports the accepted dt.
   - FDM integration, engine crank, circuits, and rigid bodies all route
     through this module (FDM first; others adopt on arrival).
   - Straight-line migration: `fdm_integration.cpp` thin wrapper around the
     shared module.

3. **Tests** — `tests/unit/solver/integrator_test.cpp`: order-of-accuracy
   checks (Euler O(h), Heun/RK4 O(h⁴) on y′=λy and y″=−y), stability of
   implicit methods at large dt, symplectic energy conservation on a harmonic
   oscillator (bounded drift), adaptive tolerance adherence; convergence
   monitor windowed residual test.

**Definition of done:** `TimeStepper`/`ConvergenceMonitor` public; all 4+ FDM
integration methods route through the shared module; integrator_test green;
`agent_guide.md` §11.5 note updated (promoted).

---

## 5. Phase B — General rigid bodies (mechanics 6-DOF) ✅ DONE

Status: `mechanics::rigid_body.*` (6-DOF state, quaternion utilities,
SymplecticEuler + RK4 advance, gravity/spring helpers) with quaternion and
rigid-body test suites. Free fall, torque-free spin, and energy-conservation
cases verified for both methods.

**Design decision recorded during implementation:** `advance(dt, loads, …)`
holds loads constant within a step (explicit in state-dependent forces).
That is the right contract for external loads (gravity/springs evaluated
per step) and SymplecticEuler behaves as intended. State-dependent forces
(gas pressure vs crank angle, aerodynamic loads) belong on the shared
integrator module's `DerivativeFn` (engine phase), whose RHS sees the full
state.

Extend `mechanics` from 1-DOF axis rotation to full rigid bodies, kept as a
*separate* module from the machine-axis stack (the 1-DOF rotation stays the
fast path for turbines/engines/compressors).

### Deliverables

1. **Types** — extend `include/exd/physics/mechanics/`:

   ```cpp
   // rigid_body.hpp
   struct RigidBodyState {
       std::array<double,3> position;      // center of mass (m)
       std::array<double,4> orientation;   // quaternion (w,x,y,z), unit
       std::array<double,3> linear_velocity;   // (m/s)
       std::array<double,3> angular_velocity;  // (rad/s, body frame)
   };
   struct RigidBodyForces {
       std::array<double,3> force;    // world frame (N)
       std::array<double,3> torque;   // body frame (N·m)
   };
   struct RigidBodyConfig {
       std::array<double,3> inertia_principal;  // principal moments (kg·m²)
       std::array<double,3> mass_center_offset; // if not at CoM
       double mass = 0.0;
   };
   class IRigidBodyDynamics {
       virtual std::string_view name() const = 0;
       virtual RigidBodyState advance(double dt, const RigidBodyForces& loads,
                                      const RigidBodyState& state,
                                      ModelStatus& status) const = 0;
   };
   ```

   - Quaternion utilities (normalize, integrate ω→q, rotate vector) as free
     functions in `src/mechanics/quaternion.cpp` (no external lib).
   - Symplectic Euler and RK4 configurations via
     `make_rigid_body_dynamics(RigidBodyIntegration, config)`.

2. **Coupling-ready surface**: rigid bodies expose a
   `RigidBodyState` getter; FDM/immersed-boundary and constraint joints
   (Phase J, engine) consume it. Motions are *driven* by forces from any
   domain (gravity, springs, fluid, contact-lite).

3. **Body assembly convenience**: constant gravity + spring force helpers
   (`mechanics::force_gravity`, `force_linear_spring`) for tests and apps.

4. **Tests** — `tests/unit/mechanics/rigid_body_test.cpp`: free fall
   (exact z(t)), torque-free rotation (|ω| conserved), quaternion
   normalization drift, symplectic orbit stability, gravity+spring
   oscillator comparison vs analytic.

**Definition of done:** 6-DOF rigid body advances with ≥2 integrators;
quaternion utilities unit-tested; energy behavior correct on test cases.

---

## 6. Phase C — Thermo + Control (enablers for engines/compressors/turbines) ✅ DONE

Status: `thermo/eos.*` (IEos + IdealGasEos, Sutherland viscosity),
`control/controller.*` (IController + PI with clamp anti-windup). EOS and PI
tested against sea-level air relations and a first-order plant (setpoint
tracking within 2%, anti-windup recovery).

### C.1 `thermo` domain — `exd::physics::thermo`

- **EOS seam** (first variant: ideal gas):

  ```cpp
  class IEos {
      virtual std::string_view name() const = 0;
      // p(ρ,T); T(p,ρ); internal energy / enthalpy; cp, cv, gamma
      virtual double pressure(double rho, double T, ModelStatus&) const = 0;
      virtual double temperature(double p, double rho, ModelStatus&) const = 0;
      virtual double gamma() const = 0;
  };
  struct IdealGasConfig { double R; double gamma; };
  std::unique_ptr<IEos> make_ideal_gas(const IdealGasConfig&);
  // Second variant later: PolytropicModel, real-gas tables (REFPROP-style).
  ```

- Files: `include/exd/physics/thermo/eos.hpp`, `src/thermo/eos.cpp`,
  `tests/unit/thermo/eos_test.cpp` (ideal gas relations, isentropic
  relations pV^γ, T_02/T_01 shock-free checks).

- Also `thermo::property` helpers: cp/cv/R lookups, Sutherland viscosity
  for pe.

### C.2 `control` domain — `exd::physics::control`

- **Controller seam** (first variant: PI):

  ```cpp
  class IController {
      virtual std::string_view name() const = 0;
      virtual double update(double setpoint, double measurement, double dt,
                            ModelStatus&) = 0;   // returns control effort
      virtual void reset() = 0;
  };
  struct PiControllerConfig { double kp; double ki; double clamp_min; double clamp_max; /* anti-windup */ };
  std::unique_ptr<IController> make_pi_controller(const PiControllerConfig&);
  // Later: PID, gain scheduling (pitch scheduling), feedforward.
  ```

- Files: `include/exd/physics/control/controller.hpp`,
  `src/control/pi_controller.cpp`, `tests/unit/control/pi_controller_test.cpp`
  (setpoint tracking on a first-order plant, anti-windup on clamp, reset).

**Definition of done:** EOS + PI controller tested; both immediately
consumable by turbine (pitch/torque control) and later engine (governor).

---

## 7. Phase D — Electrical domain (EMF / EF / MF) with coupling ✅ DONE

Status: `coupling/field_channels.*` (IVectorField3D/IScalarField3D +
structured grid channels + fluid adapters), `electrical/static_fields.*`
(3D SOR Poisson: electrostatic φ→E and magnetostatic A_z→B, verified vs
parallel-plate capacitor and Ampere's wire law), `electrical/fdtd.*` (3D
Yee FDTD with PEC box, gaussian plane source, auto-Courant dt; verified vs
c=1/√(ε₀μ₀) propagation speed, reflection, energy conservation),
`electrical/circuit.*` (DcMotorModel with exact RL armature solution +
quasi-steady `make_dc_motor_moment` as IMomentModel).

**Design decisions recorded during implementation:**
- Faraday's law requires H −= (dt/μ)·∇×E (spec typo initially; the wrong
  sign is numerically unstable — caught during implementation).
- The IMomentModel motor law is the OPPOSING convention:
  T(ω) = kt·(ke·ω − v_supply)/R. Positive = generator braking, negative =
  motoring assist. The coupled acceptance test is the physical one: a
  turbine-style constant aero torque spins a V=0 generator up from rest to
  the attracting equilibrium ω* = aero·R/(kt·ke).
- FDTD boundary is "out-of-range reads = 0" (bounded, lossless; behaves
  closer to PEC in the interior; PML/true PEC is future work).

`exd::physics::electrical`. Three nested sub-systems, all grid-friendly:

### D.1 Field channels (coupling generalization)

Generalize the sampler contracts so non-fluid domains can exchange fields:

```cpp
// coupling/field_channels.hpp (new)
class IVectorField3D {   // E, B, displacement, force density…
    virtual bool sample(const std::array<double,3>& p,
                        std::array<double,3>& value_out) const = 0;
};
class IScalarField3D {   // temperature, concentration, potential φ…
    virtual bool sample(const std::array<double,3>& p, double& value_out) const = 0;
};
// IFlowField3D stays for fluids (vector + scalar together).
```

FDTD/static EM fields implement `IVectorField3D`/`IScalarField3D`; FDM
thermal implements `IScalarField3D` (Phase I). Samplers (`sample_flow`,
`sample_scalar`) extend the pattern.

### D.2 Static EF/MF — fast win on the FDM grid

- Laplace/Poisson solvers for electric potential φ and magnetostatic
  vector/scalar potential: **reuse the existing pressure-Poisson machinery**
  (SOR solvers in `fdm_pressure.cpp`) on a regular grid.
- Files: `src/electrical/static_fields.cpp` (config/result types), tests vs
  analytic solutions (capacitor plate φ, dipole, infinite wire B).

### D.3 FDTD time-domain — `src/electrical/fdtd.cpp`

- Classic Yee lattice (E staggered, H staggered), PML/absorbing or
  conducting BC, sources (plane wave, dipole), stability on Courant limit.
- Exposes `IVectorField3D` for E and H; config/result structs per repo
  pattern; 3D grid machinery (indexing) shared with FDM.
- Tests: plane-wave propagation in vacuum (analytic), reflection at
  conducting wall, Courant stability bound, energy conservation.

### D.4 Lumped circuits + machine models — the coupling payoff

- Lumped ODE circuit state (R/L/C, sources) integrated by the Phase A
  integrator module: `src/electrical/circuit.cpp`.
- **Generator/motor models implement `mechanics::IMomentModel`** — the
  existing seam:

  ```cpp
  struct DcMotorConfig { double kt; double R; double L; double ke; double v_supply; };
  std::unique_ptr<mechanics::IMomentModel> make_dc_motor_moment(const DcMotorConfig&);
  ```

  Internal circuit state integrates i(t):  V = R·i + L·di/dt + ke·ω;
  torque = kt·i (opposing when generating). The *moment law* is the
  contract — no new coupling machinery needed.
- Tests: RL circuit analytic decay; motor spin-up from rest against a
  constant load (ω settles where kt·i(ω) = T_load); generator supplying
  torque to a turbine assembly via `RotatingAssembly`.

**Definition of done:** static EF/MF Poisson + FDTD on grids with analytic
verification; `IVectorField3D`/`IScalarField3D` channels live; DC motor as
`IMomentModel`; coupled turbine-charging-motored-rotor demo test green.

---

## 8. Phase E — Engines

`exd::physics::engine` app over the generic stack. ICE four-stroke single-
cylinder first (multi-cylinder = N cylinders summed, Phase E2).

### E.1 Slider-crank mechanism (1-DOF constraint)

- Piston translation coupled analytically to crank angle θ:
  x(θ), v(θ), a(θ) via rod length + crank radius (exact geometry, no
  integration error).
- Equivalent inertia J_eq(θ) = J_flywheel + m_piston·(dx/dθ)² + … →
  time-varying inertia handled by the Phase A integrators (ImplicitEuler
  or RK4 on θ, ω with J(θ)).
- Files: `src/engine/crank_mechanism.cpp` + tests (TDC/BDC positions,
  v(θ)=0 at TDC/BDC, energy of free-spinning crank conserved).

### E.2 In-cylinder gas force

- `thermo::IEos` (ideal gas) + heat-release model
  (`IGasForceModel` seam: `gas_force(state, crank_angle, …)` variants:
  polytropic compression/expansion, Wiebe heat release for combustion —
  second variant lands with combustion phase).
- Files: `src/engine/gas_force.cpp`, `src/engine/cycle.cpp`
  (intake/compression/power/exhaust valves modeled as flow-area/phase
  windows), tests: compression polytrope pV^γ, indicated power ≈ analytic
  Otto-cycle estimate, net cycle work sign check.

### E.3 Governor + friction

- Speed governor = Phase C PI controller (`control`) driving throttle/load.
- Friction = `mechanics::IMomentModel` (linear/bearing already available).
- Files: `src/engine/engine_simulator.cpp`: `step_engine`, `simulate_engine`
  mirroring the turbine app entry points; `tests/unit/engine/`:
  idle spin-up to governed speed, load step transient, no NaN, energy
  balance (indicated − friction − load = ΔKE).

**Definition of done:** single-cylinder engine simulates a stable governed
cycle; physics verified vs polytropic/Otto analytics; couples to load
(generator via Phase D `IMomentModel`); multi-cylinder = configuration
(Phase E2, no architecture change).

---

## 9. Phase F — Turbines: control + 3D + coupling

The turbine app already runs on the generic stack. Extend:

1. **Control integration**: pitch governor and torque controller from
   Phase C wired into `simulate_turbine` (`TurbineConfig.control`):
   rated-power regulation test (Cp held near peak, ω governed).
2. **3D FDM fluid solver** (`fluid::fdm` 3D): the current 2D solver becomes
   the 2D specialization; add 3D staggered grid, 3D pressure Poisson,
   SOR/CG; produces a `coupling::StructuredGridField` natively.
   - This is the *fluid* enabler for the coupled demo; cross-check cp/ct of
     the momentum-balance evaluator against a momentum source in FDM
     (actuator-disk-in-grid) — Phase F3.
3. **Coupled demo (reference for the framework)**: FDM 3D (or axisymmetric
   structured) field + `step_turbine` in a loop, exchanging every N steps;
   rotor state feeds back as momentum-source/moving-wall BC. This is the
   `real_solver_for_turbine.md` loop, made real.

**Definition of done:** turbine governed under control; 3D FDM fluid solver
with analytic/2D-limit verification; coupled turbine-in-grid demo runs and
matches reduced-order cp within engineering tolerance.

---

## 10. Phase G — Turbomachinery generalization (compressor capability, product-agnostic)

Re-scoped (2026-08-31): **no `compressor` app module.** The compressor use
case is met by general capabilities — the same mean-line turbomachinery
physics covers axial compressors, axial turbines, fans, and (thermodynamic
side of) rotors. Direction of work is **geometry-emergent, not a mode
enum**: `Δh₀ = u·(c_w2 − c_w1)`; swirl added → compression (work into gas),
swirl extracted → expansion (work out). Config validation warns on a
"compressor" stage that produces expansion (inlet-swirl vs blade-angle
mismatch).

**Product doctrine (whole library):** target is *80–90% engineering-grade
accuracy, simple physics coupled together, fast and easy to use*. Every
module documents its validity envelope; no speculative unification; the
measure of a model is that a caller knows where it is 80% and where it is
wrong.

### G.1 `fluid::turbomachinery` — axial mean-line stage (namespace
`exd::physics::fluid::turbomachinery`, `src/fluid/turbomachinery/`)

- `stage.hpp/.cpp` — one axial stage at mean radius. **All stage relations
  operate on TOTAL states (T0, p0)** (review C1: Euler work and the
  polytropic exponent forms are total-state laws; static recovery via
  absolute velocity incl. residual swirl; a static-path chain is an O(M²)
  error, ~16% in π at M = 0.5):
  `T02 = T01 + Δh₀/cp;` compression `π = (T02/T01)^(γ·η_p/(γ−1))`,
  expansion `π = (T02/T01)^(γ/((γ−1)·η_p))`; then static T2, p2 from the
  velocity field and `IEos`.
  - `c_a = ṁ/(ρ·A)` solved by explicit density-velocity fixed-point
    (3–8 scalar iterations, convergence guard → `ModelStatus` failure).
  - **Choking keyed on relative Mach** at the rotor LE
    `M_rel = √(c_a² + (c_θ1 − u)²)/a` (review M3): map sweep marks choked
    points invalid rather than clipping interior behavior.
  - Validity envelope (documented in header + doc): `M_rel < 0.7`,
    hub/tip ratio > 0.5, single-stage π < ~1.5, axial-first (radial stages
    are a later config'd variant, different loss correlations).
  - Entry: `StageResult solve_stage(const StageConfig&, const StageInlet&,
    double omega, double mdot, const thermo::IEos&, ModelStatus&)`.
    `StageInlet{p0, T0, c_theta}`; **`StageResult{p0, T0, c_theta, static
    p/T/rho, delta_h0, torque, power, pi, efficiencies, mach}`** — torque
    `= ṁ·Δh₀/ω` and power are first-class (review M5: system energy
    balances close on them).
- `stage_stack.hpp/.cpp` — N stages, inter-stage propagation of
  `(T0, p0, c_theta)` (residual swirl IS the next stage's inlet swirl —
  review M4); total π = Ππᵢ; reheat effect physical. `solve_stage_stack`
  pure entry, batchable.
- `operating_map.hpp/.cpp` — `solve_operating_map(stack, inlet, sweep,
  IEos, status)` sweeps ω × ṁ (dimensionless: corrected flow vs corrected
  speed, so rig-test-supplied maps are directly compatible) → `OperatingMap`
  (π/η/T surfaces, surge curve at dπ/dṁ = 0 per speed line — documented
  surrogate, real surge lines come from rigs/CFD — and relative-Mach choke
  curve). Vector grids are deliberate (config-size data, not hot path;
  agent_guide §11.1 exception noted).
- `map_lookup.hpp/.cpp` — `sample_operating_map(map, omega, mdot, status)`,
  bilinear; consumes any `OperatingMap` — computed OR rig-test data (the
  `TableLookup` variant pattern).

### G.2 `thermo::polytropic` — stagnation-polytrope math primitives

Public module `exd::physics::thermo` for the **stagnation** polytrope
family only (`polytropic_exponent`, `temp_ratio`, `press_ratio`,
efficiency relations). The engine's in-cylinder volume polytropes are
**NOT migrated** (review M2): they are documented heat-transfer stand-ins
on static p·Vⁿ and share only the word "polytropic" with total-state
turbomachinery relations — unifying them would degrade both. Two families
(volume vs stagnation) get a shared header only if a second consumer
appears; both exported as separate named relation functions.

### G.3 `fluid::lumped::plenum` — Greitzer 0D volume (namespace
`exd::physics::fluid::lumped`, `src/fluid/lumped/plenum.cpp`)

State `(p_plenum, ṁ_duct)`; plenum continuity + duct inertance
`dṁ/dt = Δp/I`; `IEos` for density (adiabatic a² = γ·R·T stated in the
header). Compressor/throttle characteristics injected as functions — the
plenum is a generic gas-network node (engine manifolds, gas networks,
HVAC later). **Stability verification uses the full two-state Jacobian**
(review C2 — the naive slope comparison is wrong in the shallow-throttle
regime): linearize, derive trace/det conditions weighted by the Greitzer
B parameter, verify against an analytic `B_critical` case AND time-march
the stable→limit-cycle transition with the Phase A integrator.

### G.4 Thin system driver (DX, not an app)

`simulate_compression_system(CompressionSystemConfig, ModelStatus&)` —
~100 lines composing `RotatingAssembly` (compressor-sense provider:
negative torque wrt +ω) + Phase D DC motor drive (motoring assist,
negative = assisting) + plenum + optional PI governor (surge margin /
speed), streaming CSV via `io::CsvSeriesWriter`. This is "easy to use +
quick output products" for the capability; the turbocharger (turbine +
compressor stages on one shaft) remains an **acceptance test** until
Phase H formalizes multi-machine coupling.

### G.5 Files, tests, verification

- Files: `include/exd/physics/fluid/turbomachinery/{stage,stage_stack,
  operating_map,map_lookup}.hpp`, `include/exd/physics/fluid/lumped/
  plenum.hpp`, `include/exd/physics/thermo/polytropic.hpp` + matching
  `src/` + `validate_*_config` per config (m1).
- `IEos` gains `density(p, T)` and a `specific_entropy`/Δs helper
  (review M1 — the seam is currently dormant infra; stage and entropy
  checks need it); consumed as non-owning const ref.
- Tests (analytic anchors, repo regime):
  - polytropic: η_p = 1 → exact isentropic `p0₂/p0₁ = (T0₂/T0₁)^(γ/(γ−1))`;
    Δs = cp·ln(T0₂/T0₁) − R·ln(p0₂/p0₁) > 0 for compression; Δh = cp·ΔT0.
  - stage: exact Euler work/torque on synthetic triangles; zero-swirl →
    Δh₀ = 0, π = 1; **full-closure reversal test** — reversed blade angles
    give identical |Δh₀|, expansion exponent, T0₂ < T0₁, π < 1, Δs > 0 in
    both senses (a reciprocal π_comp·π_turb = 1 check is NOT a test —
    compressor and turbine polytropes are not reciprocals, review C3);
    closed-form energy/kinetic/momentum closure at M_axial ≈ 0.4–0.5 vs a
    hand computation.
  - stack: single ≡ single; 2-stage reheat effect.
  - map: echoes `solve_stage_stack` at nodes; surge line at dπ/dṁ = 0;
    relative-Mach choke boundary; determinism pin.
  - plenum: steady fixed point; Jacobian stability boundary both sides;
    surge limit cycle; energy bookkeeping.
  - system: motor-driven spin-up to operating point matching the steady
    characteristic; turbocharger energy balance closed
    (`Σ(τ_motor + τ_compressor)·ω·dt = ΔKE`; steady state
    `motor_power ≈ compressor_power`); PI surge-margin setpoint tracking.

**Field channels:** deliberately none (review M6). The compressor is a
lumped machine model like the engine — no grid field to exchange;
**compressible-grid coupling is an unowned future capability**, filed as
its own roadmap item, not "Phase J".

**Definition of done:** G.1–G.3 public with analytic verification;
`simulate_compression_system` runs a motor-driven pressure rise with
EOS-consistent thermodynamics and streams CSV; turbocharger + surge
acceptance tests green; `docs/turbomachinery_architecture.md` written;
roadmap/capability matrix updated; full suite green (74 + new).

---

## 11. Phase H — Coupling framework (real exchange)

The `coupling` domain becomes the orchestrator for multi-domain runs.

1. **Rework `CouplingManager`** (keep API shape, implement the body):
   - `Link { source, target, source_channel, target_channel, interpolation
     (nearest → RBF/MLS later), interval, sub_iterations, relaxation }`.
   - `SurfaceMapper::transfer` real implementation (nearest-neighbor first,
     bilinear/trilinear where grids match).
   - Field channels (`IFlowField3D`, `IVectorField3D`, `IScalarField3D`)
     as the data currency between links.
2. **`CoupledSimulation` driver** (`src/coupling/coupled_simulation.cpp`):
   - Declare domains with per-domain dt (multi-rate), links, convergence
     criteria; explicit (staggered) mode or sub-iterated implicit mode
     with relaxation; shared `TimeStepper` per domain.
   - Solver-agnostic: each domain supplies `step(dt)` + channels (internal
     solvers adapt via thin wrappers; external solvers already expose
     `ISolverPlugin`).
3. **Reference coupled demos** (the acceptance tests):
   - Turbine ↔ FDM field (Phase F3) — fluid↔machine.
   - Turbine ↔ DC generator (Phase D) — machine↔electrical.
   - Engine ↔ generator (Phase E×D) — machine↔electrical, transient.
   - Rigid body ↔ FDM (immersed) — later, Phase J.
4. Convergence/verification: energy bookkeeping across links (checks in
   tests), staggered-vs-implicit agreement on a linear problem.

**Definition of done:** `CouplingManager` exchanges real data; demo 1 and 2
green; multi-rate staggering tested; coupling stubs removed.

---

## 12. Phase I — Rest of the major domains (grid-first)

Each follows the phase template (types → solver/config/result → channels →
tests vs analytics). FDM methods only; discretization swaps come later.

| Domain | Namespace | First solver (FDM-style) | Verification target |
|---|---|---|---|
| Thermal transfer | `thermal` | conduction (implicit CN) + advection; conjugate via bc | 1D heat equation analytic, Biot-limit check |
| Structural / soft bodies | `structural` | linear elasticity, displacement form (small strain) | beam deflection, hydrostatic patch |
| Acoustics | `acoustics` | scalar wave equation on grid | plane wave, mode frequencies of a box |
| Particles / sprays | `particles` | Lagrangian advection over sampled fields (uses channels + integrators) | ballistic motion, drag terminal velocity |
| Chemistry (combustion support) | `chemistry` | species transport + reaction ODEs (integrator module) | first-order decay, equilibrium constant |

Soft bodies additionally: mass-spring and FDM elasticity both re-use the
`IRigidBodyDynamics` integrators; FEM arrives in Phase J for the
unstructured path.

**Definition of done:** every domain has a working grid-first solver, a
public config/result pair, at least one field channel, analytic tests,
and a coupling example where adjacent domains interact (CHT, aeroacoustic,
aero-elastic via immersed force).

---

## 13. Phase J — Discretization breadth: FVM, FEM, meshes

The "prove plug-and-play" phase. Trigger rule applies: seams are extracted
here because a second variant now exists.

1. **FVM incompressible solver** on structured grids
   (`fluid::fvm`, `src/fluid/fvm/`): cell-centered velocities, face fluxes,
   Rhie-Chow or MAC staggering, reuse grid/indexing + Pressure machinery;
   config/result structs symmetric to FDM's.
2. **Method seam** — `IFluidSolver` in `coupling` or `solver`:

   ```cpp
   class IFluidSolver {        // implemented by fdm, fvm (and plugins later)
       virtual bool step(double dt) = 0;
       virtual coupling::IFlowField3D& field() = 0;
   };
   ```

   Coupled runs select the method by config
   (`fluid_method = FDM | FVM`) — first proof of plug-and-play.
3. **Cross-method verification**: lid-driven cavity and channel flow with
   both methods; grid-convergence rates; conservation checks.
4. **Unstructured meshes**: activate `mesh/` (unstructured volume,
   partition, quality) + FVM on tetra/hex; field channels adapt via
   interpolation (SurfaceMapper RBF path).
5. **FEM** (`structural`, `thermal`): element assembly on mesh/ element
   types; sparse solve (Pardiso/Eigen-lite or CG); boundary-fitted
   capability.
6. **Immersed boundary** for rigid bodies in FDM/FVM: rigid body state
   (Phase B) sources body forces into the fluid grid — the FSI demo.
7. Later: AMR, boundary-fitted moving meshes (rotor grids), spectral and
   particle methods behind the same channels.

**Definition of done:** FVM + FDM both behind `IFluidSolver` with method
switching by config; cross-method test suite green; FEM structural/thermal
with analytic patch tests; rigid body immersed coupling demo runs.

---

## 14. Cross-cutting engineering

- **Verification regime**: every solver gets (a) an analytic/limit test,
  (b) a non-dimensional sanity bound (Betz for turbines, Otto limit for
  engines, polytropic for compressors, Courant for FDTD), (c) conservation
  or energy-bookkeeping checks in coupled tests.
- **Optimizer integration**: keep `solve_x(definition, conditions, config)`
  pure and batchable (agent_guide §11.4) — new apps must expose the same
  shape.
- **Visualization**: generic field channels double as output contracts
  (VTK writers consume `IField3D` channels + mesh types).
- **Units**: exd-core has units; adopt typed quantities for new public
  config/result structs where convenient, keep doubles in hot paths.
- **Docs per module**: each new domain adds `docs/<domain>_architecture.md`
  following `bem_level3_architecture.md` as the template.
- **Checklists**: `agent_guide.md` §12 checklist applies to every new
  solver (config/result/entry/internal/factories/coupling stub/tests/CMake/
  docs).

---

## 15. Sequencing & dependencies

```
Phase A  integrators + time stepping
   │
   ├─→ Phase B  rigid bodies (6-DOF)        (uses A)
   ├─→ Phase C  thermo + control            (uses A)
   │      ├─→ Phase D  electrical/EM/FDTD   (uses A, C; channels generalize the fluid sampler)
   │      ├─→ Phase E  engines              (uses A, B, C, D)
   │      ├─→ Phase F  turbines + 3D FDM    (uses A, C, D)
   │      └─→ Phase G  turbomachinery generalization (uses A, C, D; re-scoped — §10)
   └─→ Phase H  coupling framework          (consumes D channels; hosts F/G demos)
            │
            ├─→ Phase I  thermal/structural/acoustics/particles/chemistry (grid-first)
            └─→ Phase J  FVM → method seam → FEM → meshes → immersed-boundary FSI
```

Phase A is the only hard prerequisite; B–G are largely parallel after A+C;
H formalizes what B–G already need ad hoc; J lands when the second
discretization exists.

## 15b. Wave program status (2026-08-29)

Executed: real-time output → engines → 3D FDM → coupled turbine-in-grid → capability assurance (74/74).
Next: W9 — turbomachinery generalization (Phase G re-scope).

| Wave | Deliverable | Status |
|---|---|---|
| W0 | doc corrections (collocated FDM, stub BC/mesh-IO/coupling-manager truth table) | ✅ done |
| W1 | `io` module: `IFieldWriter` (exd-fld v1 binary + timeline manifest), `CsvSeriesWriter`, `OutputPolicy`/`OutputScheduler` (injected clock) — spec in `docs/output_channels.md` | ✅ done |
| W2 | Phase E engine: `engine` module — slider-crank kinematics (analytic x,v,dx/dθ,d²x/dθ²), J_eq(θ) + ½(dJ/dθ)ω² inertia torque, polytropic + Wiebe Otto cycle, steam placeholder (admission + n=1.13 polytrope), PI governor, pure T(ω) loads, `step_engine`/`simulate_engine` + CSV machine-state streaming | ✅ done (66/66 tests) |
| W3 | Phase F.2: `fluid::fdm3` — 3D collocated SIMPLE, 6-face BCs, 7-point SOR, body-force source, persistent `FDM3Solver`, IFlowField3D adapter. Validated: uniform-flow exact, 3D Poiseuille 0.13%, Taylor–Green 4.7%, fixed-CT disk 22% | ✅ done (70+ tests) |
| W4 | Phase F.3/H-lite: coupled turbine-in-grid — `force::BladeElement` local evaluator (no induction double-count), smeared negated body force (Gaussian, β-relaxed, ramped), `turbine::run_coupled_turbine` + `default_grid_config`. Soak: ω* settles, wake decel ~20%, cp vs reduced-order ratio 0.63, energy balance closed | ✅ done (71 tests) |
| W5 | Phase F.1: PI speed governor wired into `TurbineConfig.governor` (load-fraction control, one update per step); regulation to setpoint <0.1% with mid-range throttle; batchable `simulate_engine`/`solve_fdm3` entry points stay pure | ✅ done (71 tests) |
| W7 | Real-run DX: parametric `make_turbine_definition`, rotor-state CSV streaming in `run_coupled_turbine`, `run_fdm3_simulation` stamping driver, `docs/real_run_guide.md`. Also: CSV writer creates parent dirs; Otto heat release re-anchored at TDC volume (was over-producing past the Otto bound); γ-mismatch warning | ✅ done (73 tests) |
| W8 | Capability assurance: hydro (water-turbine) soak with seawater properties; `thermo::steam` saturation EOS (Clausius–Clapeyron, enthalpies, vapor density); steam engine cycle upgraded to Rankine-lite (saturation-line temperatures, wet-steam quality, boiler-heat efficiency); `docs/capability_matrix.md` | ✅ done (74 tests) |
| W9 | Phase G re-scope (§10): product-agnostic turbomachinery — `fluid::turbomachinery` axial mean-line stage + stack (total-state closure, geometry-emergent work sign, relative-Mach choking, documented envelope), `thermo::polytropic` (stagnation family; engine polytropes NOT migrated — heat-transfer stand-ins), `fluid::lumped::plenum` (Greitzer, Jacobian-verified stability), `solve_operating_map`/`sample_operating_map`, thin `simulate_compression_system` driver (motor + plenum + PI, CSV). `IEos` + `density(p,T)` + Δs. No compressor app module | ✅ done (80/80 tests) — see docs/turbomachinery_architecture.md |
| W10 | Phase H-lite + Phase I: real `CouplingManager` exchange (`SurfaceMapper` nearest/trilinear, relaxed + implicit sub-iterated links, `CoupledSimulation` multi-rate driver — staggered vs implicit verified on a linear system) + five grid-first domains: `thermal` (steady conduction/advection, SOR), `acoustics` (leapfrog wave, box-mode frequencies), `structural` (Navier–Cauchy displacement SOR, thermal-strain channel), `particles` (Lagrangian over `IVectorField3D`), `chemistry` (mass-action + Arrhenius reactor). All analytic-verified; docs/coupling_and_grid_first_domains.md | ✅ done (87/87 tests) |
| (future) | Generic `bc` framework promotion (when a mesh-based consumer exists); `CouplingManager` real exchange; field writers for FDM3 stamps; compressible-grid coupling (unowned — lumped machine models carry no field channels by design, §10 G.5) | deferred |

## 16. Immediate next step

**W9 + W10 shipped (87/87 tests).** Next per the product doctrine
(§10): apply the now-real coupling machinery to the real domains —
conjugate heat transfer (thermal ↔ fdm3), aeroacoustics (acoustics ↔
fdm3), and the engine/turbine/compression demos migrated off ad hoc
wiring onto `CoupledSimulation`; then the remaining Phase I breadth.
Deep FVM/FEM work (Phase J) stays deliberately last.
