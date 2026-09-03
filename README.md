# extropian-physics

**One composable numerical engine.** A C++23 static library implementing the
architecture of [`docs/implementation_spec.md`](docs/implementation_spec.md)
v0.3 — a GPU-portable, composable multiphysics runtime:

```text
DATA + OPERATORS + PHYSICS + DISCRETIZATION + COUPLING + EXECUTION + HARDWARE
```

Real solvers live in `engine/physics/`; the coupling engine composes them
(fluid-structure, conjugate heat transfer); presets assemble validated
configurations (turbine, slider-crank engine, incompressible CFD, CHT) —
never duplicate solvers. CPU backend today; CUDA per spec Phase 11.
Conformance record: [`docs/implementation_spec_conformance.md`](docs/implementation_spec_conformance.md).

## Repo layout (spec §57)

```
engine/
├── include/exd/engine/
│   ├── core/            State · Field · EntitySet · Units · Memory · Operators · ExecutionGraph
│   ├── mesh/            structured topology, generation, validation, boundary conditions
│   ├── discretization/  FDM operators (gradient/divergence/curl/Laplacian/advection) + matrix-free ops
│   ├── numerics/        ODE integrators (8) · time stepping · CG/GMRES/BiCGSTAB · Newton/fixed-point
│   ├── physics/         fluid (fdm, fdm3, forces, BEM, lumped, turbomachinery) · thermal ·
│   │                    structural (static + transient) · electromagnetics · reaction ·
│   │                    particles · rigid_body · species transport · porous media ·
│   │                    thermo · control · acoustics
│   ├── coupling/        field channels · samplers · mapping · CouplingManager · contracts ·
│   │                    compatibility rules · Simulation pipeline · solver plugins
│   ├── fidelity/        FidelityProfile (REALTIME … HIGH_FIDELITY)
│   ├── presets/         turbine · engine · incompressible CFD · conjugate heat transfer
│   ├── backends/        CPU backend (parallel_for/reduce/stencil primitives)
│   ├── output/          field stamps (exd-fld) · CSV series · output policy
│   └── diagnostics/     CFL · conservation · residual history · stopwatch
├── src/                 private implementation (mirrors include/exd/engine)
└── tests/               doctest suites per layer
```

Namespace mirrors directory: `exd::engine::<layer>::<module>`. One umbrella:
`<exd/engine/engine.hpp>`. `geometry/` (spec §57) is the external
`extropian-geometry` dependency.

## Configuration pipeline (spec §54)

```text
Parse → Construct → Validate → Resolve Defaults → Build Coupling Graph →
Build Execution Graph → Allocate State → Initialize Backend → Execute
```

`coupling::Simulation` implements the pipeline; `presets::multiphysics::ConjugateHeatTransfer`
is the reference preset — two thermal modules + explicit coupling contracts
(units, rank, conservation, temporal behavior) resolved into one state, one
coupling graph and one execution graph. Usage:

```cpp
#include <exd/engine/engine.hpp>
#include <exd/engine/presets/multiphysics/conjugate_heat_transfer.hpp>

exd::engine::core::ModelStatus st;
exd::engine::presets::multiphysics::ConjugateHeatTransfer cht;
exd::engine::coupling::Simulation sim;
cht.configure({}, sim, st);      // validate → coupling graph → allocate → execute
sim.run(st);                     // converges to the 400→350→300 joined profile
```

## Capability summary (details: docs/capability_matrix.md)

| System | Run entry | Fidelity |
|---|---|---|
| Wind/water turbine | `presets::turbine::run_coupled_turbine` (3D FDM + blade-element + generator + governor) | real incompressible CFD coupling |
| Combustion/steam engine | `presets::engine::simulate_engine` (Otto, Wiebe / Rankine-lite) | real 0D/1D cycle |
| Standalone 3D CFD | `physics::fluid::fdm3::solve_fdm3` / `run_fdm3_simulation` | incompressible FDM |
| Reduced-order turbine | `physics::fluid::reduced_order::bem::solve_turbine` | BEM + corrections |
| Compressor system | `physics::fluid::turbomachinery::simulate_compression_system` | mean-line + Greitzer |
| Coupled field domains | `coupling::CoupledSimulation` + `presets::multiphysics::ConjugateHeatTransfer` | two-slab CHT, exact joined profile |
| Transient thermal (CHT-lite) | `physics::thermal::simulate_thermal` | Fourier-series match 0.2% |
| Aeroacoustic-lite | `physics::acoustics::simulate_wave` | c±u pulse arrival within 5% |
| FSI-lite (W12) | `coupling::simulate_drag_body` | terminal v_t ±10%, drag→mg ±5% |
| Matrix-free Poisson | `numerics::solve_cg/gmres/bicgstab` on `discretization::fdm::FdmLaplacianOperator` | analytic sin³ verification |
| Species transport | `physics::species::solve_species` (advection–diffusion–reaction, operator-composed) | exact decay, A→B conservation, advective-decay profile, variance 2Dt |
| Elastic waves | `physics::structural::solve_elasticity` (transient) | P-wave arrival within 2%, flight energy conserved |
| Porous media | `physics::porous::solve_porous` (Darcy, implicit CG + direct steady) | exact linear steady, variance 2Kt, linear mass growth |
| Aeroacoustics | `presets::multiphysics::run_aeroacoustics` (fdm3 → acoustics mean flow) | pulse arrives on schedule, mean flow injected |
| Thermal stress | `presets::multiphysics::run_thermal_stress` | free-bar expansion u(L) = α·g·L²/2 |
| Joule heating | `presets::multiphysics::run_joule_heating` (EM → thermal q-channel) | source/outflux energy balance, heating confirmed |
| Species in flow | `presets::multiphysics::run_species_in_flow` | outlet c = c_in·exp(−k·L/u) within 5% |
| Robotic arm | `physics::robotics` (2R, revolute joints, torque limits, stops, PD) + `presets::robotics::run_arm_trajectory` | energy conserved, PD tracking < 1e-3 rad, stops bound motion |
| Chiplet board | `presets::electronics::solve_chiplet_board` (heterogeneous-k conduction, chips, spreaders) | chip power = sink flux, peak scales with power, spreader cools |
| Electrostatics | `physics::electromagnetics` Neumann faces | exact linear discrete bridge (1e-7) |
| Fidelity profiles | `fidelity::profile(FidelityLevel::…)` | REALTIME→HIGH_FIDELITY defaults |

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
107 unit tests (doctest) — the full suite runs in ~60 s.

## Docs

- `docs/implementation_spec.md` — the authoritative architecture spec (v0.3).
- `docs/implementation_spec_conformance.md` — refactor contract: layout/namespace mapping, spec§→code table, phase status.
- `docs/capability_matrix.md` — per-system fidelity and verification anchors.
- `docs/real_run_guide.md` — run recipes.
- `docs/output_channels.md` — output contract for the animation/visualization repo.
- `docs/modular_solver_architecture.md` — historical roadmap (phases A–J, waves W0–W12), superseded by the spec.

## License

Business Source License 1.1 — see [LICENSE](LICENSE).
Converts to Apache 2.0 on 2029-05-26.
