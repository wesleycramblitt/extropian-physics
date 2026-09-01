# extropian-physics

**Simulation physics library: data model, interfaces, and real solvers.**

A C++23 static library implementing real multiphysics solvers with a
pluggable coupling architecture: a 3D incompressible fluid solver, a
coupled turbine-in-grid (wind *and* water), slider-crank engines
(combustion and steam), reduced-order blade-element models, machine
mechanics, thermo/EOS, control, electrical, and real-time output channels.
High-fidelity external solvers (FluidX3D, OpenFOAM, CalculiX…) implement
`ISolverPlugin` in separate `extropian-solver-*` repos.

Depends on `extropian-core` and `extropian-geometry`.

**Where to start**

- Capability overview per system (what runs today, fidelity, exact config): [`docs/capability_matrix.md`](docs/capability_matrix.md)
- Run recipes (turbine, engine, standalone CFD): [`docs/real_run_guide.md`](docs/real_run_guide.md)
- Output contract for the animation/visualization repo: [`docs/output_channels.md`](docs/output_channels.md)
- Whole-system roadmap (phases A–J): [`docs/modular_solver_architecture.md`](docs/modular_solver_architecture.md)

## Module map

```
include/exd/physics/ and src/
├── fluid/
│   ├── fdm/                   2D incompressible FDM (legacy, untouched)
│   ├── fdm3/                  3D incompressible FDM — SIMPLE, SOR, 6-face BCs,
│   │                          body-force sources, persistent FDM3Solver
│   │                          [docs/fdm3_architecture.md]
│   ├── forces/                blade surfaces + force evaluators:
│   │                          PressureIntegration, MomentumBalance (BEM),
│   │                          TableLookup, BladeElement (local, coupled)
│   └── reduced_order/bem/     Level-3 duct/hull-coupled BEM turbine solver
│                              [docs/bem_level3_architecture.md]
├── turbine/                   turbine app: step/simulate over the generic
│   │                          stack, parametric builder, coupled turbine-in-grid
│   │                          driver [docs/turbine_coupling_architecture.md]
├── turbomachinery/            product-agnostic axial mean-line: stage/stack
│   │                          (geometry-emergent compressor/turbine sense),
│   │                          operating-map solver + sampler, compression-system
│   │                          driver [docs/turbomachinery_architecture.md]
├── lumped/                    Greitzer plenum — generic 0D gas-volume node
│   │                          (2-state surge cell, injectable characteristics)
├── engine/                    slider-crank engines: Otto (Wiebe) + steam
│   │                          (Rankine-lite), governor, CSV motion output
│   │                          [docs/engine_architecture.md]
├── mechanics/                 1-DOF rotors, 6-DOF rigid bodies, moment models,
│   │                          quaternions, rotating assemblies
├── thermo/                    EOS (ideal gas), Sutherland viscosity, steam
│   │                          saturation model
├── control/                   PI controller (anti-windup)
├── electrical/                static EF/MF (SOR Poisson), 3D FDTD, DC motor
├── coupling/                  field channels (IFlowField3D, vector/scalar),
│   │                          samplers, structured-grid adapters, REAL
│   │                          CouplingManager (Jacobi-ordered exchange,
│   │                          target read-back) + SurfaceMapper +
│   │                          CoupledSimulation multi-rate driver
│   │                          [docs/coupling_and_grid_first_domains.md]
├── thermal/                  transient implicit conduction/advection with
│   │                          velocity-channel transport (CHT-lite),
│   │                          temperature/read channels — IScalarField3D
├── acoustics/                 leapfrog wave + mean-flow (convected) wave —
│   │                          IScalarField3D channels
├── structural/                linear elasticity (Navier–Cauchy displacement
│   │                          SOR, thermal-strain channel) — IVectorField3D
├── particles/ chemistry/      Lagrangian tracking over channels; 0D mass-action
│   │                          + Arrhenius reactor (lumped — no channels)
├── io/                        output channels: exd-fld v1 field stamps,
│   │                          CSV time series, cadence policy
│   │                          [docs/io_architecture.md]
├── solver/                    plugin interface, SolverManager, shared
│   │                          integrators (8 methods), time stepping
├── mesh/ field/ bc/ material/ data-model infrastructure (mesh I/O stubs —
│                              real structured/FVM/FEM meshing is Phase J)
├── tools/                     (removed — no demo executables; runs are
│                              caller programs, see real_run_guide)
```

Every domain follows the repo doctrine: typed config/result structs,
`ModelStatus` error channels (no exceptions), strategy+factory for pluggable
models, value-type results, deterministic and optimizer-batchable entry
points. Conventions: [`docs/agent_guide.md`](docs/agent_guide.md).

## Capability summary (details: capability_matrix.md)

| System | Run entry | Fidelity |
|---|---|---|
| Wind turbine | `turbine::run_coupled_turbine` (3D FDM + blade-element + generator + governor) | real incompressible CFD coupling |
| Water turbine | same, `rho=1025, mu=1.08e-3` (hydro, proven) | same |
| Combustion engine | `engine::simulate_engine` (Otto, Wiebe) | real 0D/1D cycle |
| Steam engine | `engine::simulate_engine` (Rankine-lite) | real 0D/1D cycle + saturation EOS |
| Standalone 3D CFD | `fluid::fdm3::solve_fdm3` / `run_fdm3_simulation` | incompressible FDM |
| Reduced-order turbine | `bem::solve_turbine` | BEM + corrections |
| Compressor system | `fluid::turbomachinery::simulate_compression_system` (stage stack + plenum + DC motor + governor) | mean-line + Greitzer |
| Turbocharger (balance) | acceptance test: compressor + turbine stages on one shaft | mean-line, one code path |
| Coupled field domains | `coupling::CoupledSimulation` (manager + surface mapper + multi-rate driver), two-slab CHT on transient `thermal` | analytic fixed point (350 K interface, exact joined profile) |
| Transient thermal (CHT-lite) | `thermal::simulate_thermal` (implicit, velocity-channel advection, read channel) | Fourier-series match 0.2%, effective-Péclet 0.3% |
| Aeroacoustic-lite | `acoustics::simulate_wave` with uniform mean flow | c±u pulse arrival within 5% |
| Grid-first domains | `thermal`/`acoustics`/`structural`/`particles`/`chemistry` — see rows above + capability_matrix | analytic-verified (W10, W11) |
| FSI-lite (W12) | `coupling::simulate_drag_body`: 6-DOF body ⇄ fdm3, smeared point force + upstream drag probe, preconditionable flow | terminal v_t ±10%, drag→mg ±5%, exact normalization, deterministic |
| CHT-lite demo (W12) | steady duct flow advects thermal energy via the W11 velocity channel | mean outlet T = 1D estimate ±5%; conduction control exact |

## Level-3 BEM turbine solver (reduced-order)

A single-rotor, axisymmetric Blade-Element/Momentum solver in
`exd::physics::fluid::reduced_order::bem` — the standalone fast model used
for design-space sweeps, initialization, and reduced-order Cp reference.

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

Returns rotor torque/thrust/power, Cp/Ct, radial loading, hull drag with a
viscous/pressure split, and an engineering-approximate axial
velocity/pressure field. Correction models: induction
(Standard/Buhl, GlauertIterative, Snel) × loss (Prandtl, DuSelig,
Chaviaropoulos). Engineering-estimate caveats are documented in
`docs/bem_level3_architecture.md`.

## ISolverPlugin

The central interface for external high-fidelity solvers:

```cpp
class ISolverPlugin {
    virtual void initialize(mesh, bcs, materials, params) = 0;
    virtual bool step(double dt) = 0;
    virtual void finalize() = 0;
    virtual unique_ptr<FieldAccessor> get_field(name) = 0;
    virtual unique_ptr<CouplingSurface> get_coupling_surface(name) = 0;
};
```

Solver plugins are separate repos (`extropian-solver-fluidx3d`, etc.).

## Building

```bash
cmake -S . -B build -DEXT_PHYSICS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Local dependency overrides:

```bash
cmake -S . -B build -DEXT_PHYSICS_BUILD_TESTS=ON \
    -DEXD_GEOMETRY_DIR=/path/to/extropian-geometry \
    -DEXD_CORE_DIR=/path/to/extropian-core
```

Requires: CMake 3.21+, C++23, `extropian-core`, `extropian-geometry`.
89 unit tests (doctest) — the full suite runs in ~1 min.

## License

Business Source License 1.1 — see [LICENSE](LICENSE).
Converts to Apache 2.0 on 2029-05-26.
