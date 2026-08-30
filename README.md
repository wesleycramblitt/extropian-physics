# extropian-physics

**Simulation physics library: data model, interfaces, and engineering solvers.**

Contains the mesh data structures, field types, boundary condition framework, material database, solver plugin interface, multiphysics coupling, and a growing set of reduced-order physics models. High-fidelity solver implementations live in separate `extropian-solver-*` repos.

Depends on `extropian-core` and `extropian-geometry`.

## Category taxonomy

```
include/exd/physics/ and src/
├── fluid/                      ← fluid domain
│   ├── cfd/                    ← future: fvm/, fem/, lbm/
│   └── reduced_order/          ← analytical / engineering solvers
│       ├── bem/                ← Level-3 duct/hull-coupled BEM turbine solver
│       ├── actuator_disk/      ← future
│       ├── vortex/             ← future
│       └── potential_flow/     ← future
├── thermal/                    ← future
├── structural/                 ← future
├── electromagnetics/           ← future
└── shared/                     ← mesh/, field/, bc/, material/, solver/, coupling/
```

Shared infrastructure (`mesh`, `field`, `bc`, `material`, `solver`, `coupling`) stays at the top level of the library and is consumed by every domain.

## Level-3 BEM turbine solver

A single-rotor, axisymmetric Blade-Element/Momentum solver in `exd::physics::fluid::reduced_order::bem`.

```cpp
#include <exd/geometry/turbine.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>

using namespace exd::physics::fluid::reduced_order::bem;

PolarDatabase polars;
polars.add_builtin_polars();          // synthetic NACA0012 / NACA4412 tables

BEMSolverConfig config;
config.element_count = 32;
config.k_duct = 0.5;                  // duct acceleration coefficient
config.hull_cd = 0.2;                 // total hull drag coefficient

OperatingConditions cond;
cond.v_inf = 10.0;
cond.rho   = 1.225;

TurbineResult result = solve_turbine(turbine, cond, polars, config);
```

The solver returns rotor torque/thrust/power, Cp/Ct, radial loading, hull drag with a viscous/pressure split, and an engineering-approximate axial velocity/pressure field.

### Correction models

The solver supports pluggable induction and loss correction models via `BEMSolverConfig`:

```cpp
BEMSolverConfig config;
config.induction_correction = InductionCorrection::GlauertIterative;  // or Standard, Snel
config.loss_correction = LossCorrection::DuSelig;                     // or Prandtl, Chaviaropoulos
```

| Model | Type | Description |
|-------|------|-------------|
| `Standard` | Induction | Buhl closed-form (default) |
| `GlauertIterative` | Induction | Iterative Glauert empirical correction |
| `Snel` | Induction | Exponential blending, smooth transition |
| `Prandtl` | Loss | Classic Prandtl tip/hub loss (default) |
| `DuSelig` | Loss | Du-Selig (1993) modified tip loss |
| `Chaviaropoulos` | Loss | Chaviaropoulos-Hansen (2000) loading-dependent loss |

3×3 = 9 possible model combinations. See `docs/bem_level3_architecture.md` §15.

### Optimizer integration

The BEM solver integrates with [extropian-optimization](https://github.com/wesleycramblitt/extropian-optimization) via its stepper API:

```cpp
#include <exd/opt/opt.hpp>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>

using namespace exd::opt;
using namespace exd::physics::fluid::reduced_order::bem;

Problem p;
p.variables = {
    Variable{0.01, 0.5},    // chord (m)
    Variable{-5.0, 15.0},   // twist offset (deg)
    Variable{100.0, 2000.0} // RPM
};
p.directions = {Direction::Maximize};

OptimizeOptions o;
o.max_evaluations = 5000;
o.seed = 42;

Optimizer opt(p, Algo::CMAES, o);
while (opt.running()) {
    auto batch = opt.request_batch();
    if (batch.empty()) break;
    std::vector<Evaluation> evals;
    for (const auto& x : batch) {
        auto eng = to_engineering(p, x);
        // Modify TurbineDefinition with eng[0], eng[1], eng[2]
        TurbineResult result = solve_turbine(turbine, cond, polars, config);
        evals.push_back({{-result.rotor.cp}, {}, result.valid});
    }
    opt.submit_results(std::move(evals));
}
auto r = opt.result();  // r.best_x has optimal turbine parameters
```

### Engineering-estimate caveats

- **Duct model**: `M_duct = 1 + K_duct·(A_u/A_r − 1)` is a one-dimensional area-ratio acceleration estimate, not a panel or CFD solve.
- **Wake model**: Gaussian wake deficit spreading downstream; the far wake does not recover to freestream inside the grid.
- **Pressure field**: Bernoulli pressure plus the per-element momentum pressure jump; total pressure is not conserved downstream. This is **not a CFD pressure solution**.
- **Chord**: approximated by the meridional distance `|TE − LE|`, mirroring `exd-geometry`.
- **Polars**: the built-in NACA tables are synthetic engineering approximations and must be replaced with XFoil or wind-tunnel data for real studies.

### Current status

- BEM phases 1–3 (induction/loss correction models) implemented.
- Modular stack (shared integrators, rigid bodies, thermo, control, electrical, turbine app, coupling field channels) implemented — see `docs/modular_solver_architecture.md`.
- **`fluid::fdm3`** — 3D collocated SIMPLE solver (6-face BCs, SOR Poisson, body-force sources, persistent `FDM3Solver`, IFlowField3D adapter). Verified: uniform-flow preservation, 3D Poiseuille, Taylor–Green decay, fixed actuator disk vs momentum theory.
- **`turbine::run_coupled_turbine`** — actuator-disk turbine-in-grid coupling (local blade-element forces, smeared negated body force, under-relaxation + ramp, spin-up soak matching reduced-order Cp within engineering tolerance).
- **`engine`** — single-cylinder slider-crank: analytic kinematics with J_eq(θ) and the ½(dJ/dθ)ω² inertia torque, polytropic + Wiebe Otto cycle, PI governor, CSV machine-state output (crank angle, piston x/v, p/T, torque, power over time). Steam variant is Rankine-lite: Clausius–Clapeyron saturation EOS (`thermo::steam`), saturated admission with cutoff, wet-steam polytrope expansion, boiler-heat efficiency accounting.
- **Capability matrix** — wind/water turbine, combustion engine, steam engine: what runs, at what fidelity, and the exact configuration — `docs/capability_matrix.md`.
- **`io`** — real-time output channels: binary `exd-fld v1` field stamps + timeline manifest (`docs/output_channels.md`), CSV time series, wall-clock-throttled `OutputPolicy` for "real-time if specified". Contract ready for the animation/visualization repo.
- Speed governor in `TurbineConfig` (PI load-fraction control).

## Architecture

```
ext::physics::
├── mesh/          Unstructured volume, surface, structured meshes + I/O
├── field/         Scalar, vector, tensor fields + interpolation
├── bc/            Boundary condition types + serialization
├── material/      Material property database
├── solver/        ISolverPlugin interface + manager + time stepping
└── coupling/      Surface mapping, coupling orchestration
```

## ISolverPlugin

The central interface. Every solver (FluidX3D, OpenFOAM, CalculiX, Elmer, etc.) implements this:

```cpp
class ISolverPlugin {
    virtual void initialize(mesh, bcs, materials, params) = 0;
    virtual bool step(double dt) = 0;
    virtual void finalize() = 0;
    virtual unique_ptr<FieldAccessor> get_field(name) = 0;
    virtual unique_ptr<CouplingSurface> get_coupling_surface(name) = 0;
};
```

Solver plugins are separate repos (`extropian-solver-fluidx3d`, `extropian-solver-openfoam`, etc.).

## Quick usage sketches

**Engine with CSV motion output (for animating the piston/crank in another repo):**

```cpp
#include <exd/physics/engine/engine_simulator.hpp>
using namespace exd::physics::engine;
exd::physics::ModelStatus status;
EngineConfig cfg;                       // defaults: 4-stroke Otto, Wiebe heat release
cfg.thermo.q_in_cycle = 1500.0;         // J per cycle
cfg.initial_omega = 50.0;               // starter momentum
cfg.csv_path = "engine_state.csv";      // time,theta,omega,piston_x,piston_v,p_cyl,...
simulate_engine(cfg, status);
```

**Coupled turbine-in-grid (rotor inside the 3D FDM, fields streamable at a cadence):**

```cpp
#include <exd/physics/turbine/coupled_turbine.hpp>
using namespace exd::physics::turbine;
CoupledTurbineConfig c;
c.grid = default_grid_config(10.0);     // Inlet at +z with w=−vInf, Outlet at −z
c.turbine = my_turbine;                 // exd-geometry TurbineDefinition
c.rotor_origin = {1.2, 1.2, 2.5};
c.rotor_inertia = 100.0;
auto r = run_coupled_turbine(c, status); // r.history: t, ω, θ, torque, power, exchanges
```

## Building

```bash
cmake -S . -B build -G Ninja
cmake --build build
```

To build tests with local dependency overrides:

```bash
cmake -S . -B build -DEXT_PHYSICS_BUILD_TESTS=ON \
    -DEXD_GEOMETRY_DIR=/path/to/extropian-geometry \
    -DEXD_CORE_DIR=/path/to/extropian-core
cmake --build build
ctest --test-dir build
```

Requires: `extropian-core`, `extropian-geometry`.

## License

Business Source License 1.1 — see [LICENSE](LICENSE).
Converts to Apache 2.0 on 2029-05-26.
