# Implementation-Spec Conformance — Refactor Plan & Record

**Status: COMPLETE — repo conforms to implementation_spec.md v0.3 (CPU today; GPU/FVM/FEM = roadmap)**
**Spec: `docs/implementation_spec.md` v0.3 (authoritative)**
**This document is the contract for the refactor: layout, namespaces, decisions,
migration record, and the spec→code conformance map.**

---

## 1. Objective

Restructure and extend `extropian-physics` so the repository **matches
`docs/implementation_spec.md` exactly** in architecture, terminology, and
layering, while preserving every existing capability (89 tests) and the repo
doctrine (no exceptions, `ModelStatus` error channels, typed config/result
structs, deterministic batchable entry points).

The spec's phased roadmap (Phases 1–13) is the completion definition:

| Spec phase | Content | Status |
|---|---|---|
| 1 — Core Runtime | State, Field, EntitySet, Memory, Operators, Execution Graph, CPU Backend | **this refactor** |
| 2 — Structured Mesh | Cartesian mesh, topology, geometry, boundary identification, mesh generation | **this refactor** |
| 3 — FDM | Stencil operators, BCs, heat/diffusion/advection-diffusion | **this refactor** |
| 4 — Basic Physics | Thermal, Fluid, Reaction-Diffusion | mostly exists (W10–W12); re-homed |
| 5 — Numerical Solvers | ODE, time integrators, CG, GMRES, BiCGSTAB, nonlinear | ODE exists; linear/nonlinear **this refactor** |
| 6 — FVM | face topology, fluxes, pressure-velocity coupling | roadmap (Phase J of old plan) |
| 7 — FEM | elements, basis, quadrature, matrix-free | roadmap |
| 8 — Bodies & Particles | rigid bodies, particles, neighbor search, constraints | rigid bodies exist; SoA sets + neighbor search roadmap |
| 9 — Coupling Engine | coupling graph, contracts, mapping, interpolation, conservative transfer, temporal coupling, iteration, validation | contracts + validation **this refactor**; conservative transfer roadmap |
| 10 — Presets | CFD, FEA, Thermal, EM, Multiphysics from modules | **this refactor** (turbine/engine demoted) |
| 11 — GPU | CUDA backend, GPU-resident state | roadmap |
| 12 — GPU-aware external interfaces | GpuContext/GpuBuffer/GpuFieldView | roadmap |
| 13 — Fidelity & Accuracy | FidelityProfile, diagnostics, error estimation | profiles + diagnostics **this refactor**; estimation roadmap |

---

## 2. Repository layout (spec §57, adapted to this repo's include/src split)

```
extropian-physics/
├── engine/                                  ← the engine (spec §57 root)
│   ├── include/exd/engine/                  ← PUBLIC headers (namespace == path)
│   │   ├── engine.hpp                       umbrella (was include/exd/physics/physics.hpp)
│   │   ├── core/                            §57 core/
│   │   │   ├── model_status.hpp             ModelStatus (was include/exd/physics/model_status.hpp)
│   │   │   ├── units.hpp                    dimensional units + validation (§25; exd-core units too thin)
│   │   │   ├── field.hpp                    Field + FieldMetadata (§6)
│   │   │   ├── entity_set.hpp               EntitySet: cells/nodes/faces/particles/bodies (§7)
│   │   │   ├── state.hpp                    State container (§5), CPU/GPU residency flags
│   │   │   ├── memory.hpp                   MemorySpace/Buffer (§59)
│   │   │   ├── operator.hpp                 operator requirements/outputs contract (§34)
│   │   │   └── execution.hpp                ExecutionGraph (DAG) + primitives (§38)
│   │   ├── mesh/                            §57 mesh/
│   │   │   ├── structured.hpp               StructuredGrid (promoted from coupling) + topology (§9)
│   │   │   ├── generation.hpp               uniform/local refinement (§10)
│   │   │   ├── validation.hpp               mesh validation (§10)
│   │   │   └── boundary.hpp                 boundary identifiers + BC framework (§32)
│   │   ├── discretization/                  §57 discretization/
│   │   │   └── fdm/                         FDM stencil operators, ghost cells (§26)
│   │   │       ├── gradient.hpp divergence.hpp curl.hpp laplacian.hpp advection.hpp
│   │   ├── numerics/                        §57 numerics/
│   │   │   ├── integrators.hpp              ODE framework: Euler…AdaptiveRK45 (§29)
│   │   │   ├── time_stepping.hpp            TimeStepper (CFL) + ConvergenceMonitor
│   │   │   ├── linear_operator.hpp          matrix-free operator interface (§35)
│   │   │   ├── cg.hpp bicgstab.hpp gmres.hpp (§36)
│   │   │   └── nonlinear.hpp                fixed point + Newton (§37)
│   │   ├── coupling/                        §57 coupling/
│   │   │   ├── field_channels.hpp  field_sampler.hpp  surface_mapping.hpp  coupling_manager.hpp
│   │   │   ├── contracts.hpp                coupling contracts: units/rank/conservation/temporal (§17–19)
│   │   │   ├── rules.hpp                    CompatibilityRule validation (§21–22, §55)
│   │   │   ├── plugin_interface.hpp         ISolverPlugin (external solver contract)
│   │   │   └── solver_manager.hpp           SolverManager (plugin registry)
│   │   ├── physics/                         §57 physics/
│   │   │   ├── fluid/                       fdm, fdm3, forces, lumped, reduced_order/bem, turbomachinery
│   │   │   ├── thermal/  structural/  acoustics/  particles/
│   │   │   ├── electromagnetics/            (was electrical/)
│   │   │   ├── reaction/                    (was chemistry/)
│   │   │   ├── rigid_body/                  (was mechanics/)
│   │   │   ├── thermo/  control/
│   │   ├── fidelity/                        §57 fidelity/
│   │   │   └── profiles.hpp                 FidelityProfile + REALTIME…HIGH_FIDELITY (§45–46)
│   │   ├── presets/                         §57 presets/ (assemblies, NOT solvers — §14, §69)
│   │   │   ├── cfd/                         incompressible CFD assembly (default grid + stamping)
│   │   │   ├── multiphysics/                conjugate heat transfer preset; FSI wrappers
│   │   │   ├── turbine/                     (was physics/turbine app) — demoted to preset
│   │   │   └── engine/                      (was physics/engine app) — demoted to preset
│   │   ├── backends/                        §57 backends/
│   │   │   └── cpu.hpp                      CPU backend: parallel_for/reduce/stencil/… (§40)
│   │   ├── output/                          §57 output/ (was physics/io)
│   │   │   ├── field_writer.hpp  series_writer.hpp  output_policy.hpp
│   │   └── diagnostics/                     §57 diagnostics/
│   │       └── diagnostics.hpp              residual/CFL/conservation/performance (§50)
│   ├── src/                                 private implementation (mirrors include/exd/engine)
│   ├── tests/                               §57 tests/ (doctest; unit suites per layer)
│   └── CMakeLists.txt                       library target exd::physics
├── docs/  data/  build/                     (unchanged)
```

`geometry/` (spec §57) is satisfied by the **external dependency
`extropian-geometry`** (`exd::geometry`): primitives, CAD-like interfaces,
implicit geometry. It stays a separate repository — documented mapping, no
code lives here.

---

## 3. Namespace mapping (namespace mirrors directory)

Root namespace changes from `exd::physics` to **`exd::engine`** (the spec's
tree root is `engine/`; the dependency library already owns `exd::core`,
`exd::math` — a collision is avoided).

| Old namespace | New namespace | Directory |
|---|---|---|
| `exd::physics` (ModelStatus, umbrella) | `exd::engine::core` / `exd::engine` (umbrella) | core/, engine.hpp |
| `exd::physics::thermal` | `exd::engine::physics::thermal` | physics/thermal/ |
| `exd::physics::fluid` (+ fdm/fdm3/forces/lumped/turbomachinery/reduced_order::bem) | `exd::engine::physics::fluid` (…) | physics/fluid/… |
| `exd::physics::acoustics`, `structural`, `particles`, `thermo`, `control` | `exd::engine::physics::<same>` | physics/<same>/ |
| `exd::physics::electrical` | `exd::engine::physics::electromagnetics` | physics/electromagnetics/ |
| `exd::physics::chemistry` | `exd::engine::physics::reaction` | physics/reaction/ |
| `exd::physics::mechanics` | `exd::engine::physics::rigid_body` | physics/rigid_body/ |
| `exd::physics::coupling` (+ drag_body_solver) | `exd::engine::coupling` | coupling/ |
| `exd::physics::solver` (integrators, time_stepping) | `exd::engine::numerics` | numerics/ |
| `exd::physics::solver` (plugin_interface, solver_manager) | `exd::engine::coupling` | coupling/ |
| `exd::physics::io` | `exd::engine::output` | output/ |
| `exd::physics::turbine` (app) | `exd::engine::presets::turbine` | presets/turbine/ |
| `exd::physics::engine` (app) | `exd::engine::presets::engine` | presets/engine/ |
| (new) | `exd::engine::mesh`, `exd::engine::core`, `exd::engine::discretization::fdm`, `exd::engine::fidelity`, `exd::engine::backends`, `exd::engine::diagnostics` | mesh/, core/, discretization/fdm/, fidelity/, backends/, diagnostics/ |

The old public umbrella `exd/physics/physics.hpp` is replaced by
`exd/engine/engine.hpp`. No alias headers are kept (single clean tree; all
consumers are in-repo and are updated).

---

## 4. Key decisions

1. **Turbine and engine "apps" are NOT necessary as solver modules — they become
   presets.** Rationale (spec §14, §69, §74): the engine must never contain
   monolithic per-application solvers; the same capability must be reachable by
   composing modules. The turbine case = `mechanics`(rigid_body)
   + `fluid::forces`/BEM + `control` + `electromagnetics` (generator)
   + `fluid::fdm3` (field) composed by `coupling`; the engine case = rigid-body
   (crank) + `thermo` + `control` + `mechanics` moments composed by the ODE
   framework. The physical building blocks stay in `physics/`; the drivers
   (`simulate_engine`, `run_coupled_turbine`, …) move to `presets/` unchanged
   in behavior — *every* previously solvable case remains solvable, and the
   existing acceptance tests follow the drivers (now preset tests).
2. **Namespace root `exd::engine`** — matches the spec tree root, avoids the
   `exd::core`/`exd::math` collision with the dependency.
3. **Dead headerless stubs are deleted, not migrated.** `src/field/*`,
   `src/mesh/*`, `src/bc/*`, `src/material/*` have zero consumers (inventory in
   `docs/modular_solver_architecture.md` §2). They are replaced by the real
   `core`/`mesh` modules built in this refactor (Field, structured mesh
   topology, boundary framework, material properties).
4. **Units** (spec §25): `exd::core::units` is angle-conversion only, so a
   dimensional unit type (SI base exponents + quantity constants + consistency
   validation) is implemented in `engine/core/units.hpp`.
5. **Physics stays hardware-independent** (§3.3): the CPU backend implements
   the execution primitives; physics modules never name a backend. Kernel
   entry points take backend + fields.
6. **Coupling contracts** (§17–21): `CouplingLink` is upgraded to carry a
   contract (quantity, units, rank, association, conservation) and the
   validation system (`rules.hpp`) enforces compatibility before execution.
7. **The old roadmap** (`docs/modular_solver_architecture.md`, phases A–J,
   waves W0–W12) remains as historical record. The spec's phased plan
   supersedes it; this document is the conformance record.

---

## 5. Spec-section → code map (completeness checklist)

| Spec § | Requirement | Where (after refactor) |
|---|---|---|
| 3.1 | data-oriented arrays | `core/field.hpp` (storage = contiguous buffer), all existing solvers |
| 3.2 | operators over fields | `core/operator.hpp`, `discretization/fdm/*` |
| 3.3 | physics independent of hardware | all `physics/` modules (no backend includes) |
| 3.4, 12, 13 | physics separate from discretization, per-module compatibility | module metadata: `physics/<mod>/module.hpp` declares supported/preferred discretizations |
| 5 | State | `core/state.hpp` (fields, residency, versioning) |
| 6 | Field metadata | `core/field.hpp` (name/rank/components/units/location/association/precision/domain/mesh) |
| 7 | EntitySet | `core/entity_set.hpp` (cells/nodes/faces/particles/rigid bodies, contiguous) |
| 8 | Geometry | external `extropian-geometry` (§2 mapping) |
| 9–10 | Mesh + generation | `mesh/structured.hpp`, `mesh/generation.hpp`, `mesh/validation.hpp` |
| 25 | Units | `core/units.hpp` |
| 26 | FDM operators | `discretization/fdm/*` (gradient/divergence/curl/Laplacian/advection, ghost cells) |
| 29 | ODE framework | `numerics/integrators.hpp` (8 methods ≥ spec's 4) |
| 30 | Rigid bodies | `physics/rigid_body/` (existing) + SoA set (roadmap) |
| 31 | Particles | `physics/particles/` (existing tracks; neighbor search roadmap) |
| 32 | BCs first-class | `mesh/boundary.hpp` |
| 33 | Materials | `core/material.hpp` (promoted from stub, made public) |
| 34–35 | Operators + matrix-free | `core/operator.hpp`, `numerics/linear_operator.hpp` (apply/apply_transpose/diagonal/jvp) |
| 36–37 | Linear/nonlinear solvers | `numerics/cg.hpp`, `bicgstab.hpp`, `gmres.hpp`, `nonlinear.hpp` |
| 38 | Execution graph | `core/execution.hpp` (+ `CoupledSimulation` macro-step driver) |
| 39–40 | GPU-first / backend abstraction | `backends/cpu.hpp`; CUDA = roadmap (Phase 11) |
| 41–42, 44 | GPU context, external consumers | interface reserved (`core/state.hpp` residency flags, coupling channels); roadmap |
| 43 | Persistent output | `output/` (was io) |
| 45–46 | Fidelity profiles | `fidelity/profiles.hpp` |
| 48–50 | Accuracy/diagnostics | `diagnostics/diagnostics.hpp` (residual/CFL/conservation); estimation roadmap |
| 14, 51–53, 68–70 | Presets | `presets/*` (turbine, engine, incompressible CFD, CHT, FSI) |
| 54–55 | Configuration pipeline + validation | `core/execution.hpp` (Simulation: validate→coupling graph→execution graph→allocate→init backend→execute), `coupling/rules.hpp` |
| 56 | Dependency direction | enforced by layout: core→numerics→discretization→physics→coupling→presets; backends orthogonal |
| 63–64 | Verification | doctest suites per module (89 existing + new infra tests) |

---

## 6. Migration record

| Step | Content | Status |
|---|---|---|
| M1 | Write this conformance doc | done |
| M2 | git-mv existing code into `engine/` tree; rewrite includes/namespaces/CMake; restore green build + 89 tests | ✅ done |
| M3 | Core runtime (Phase 1): Field/State/EntitySet/Memory/Units/Operators/ExecutionGraph + CPU backend | ✅ done |
| M4 | Mesh layer: structured topology, generation, validation, boundary | ✅ done |
| M5 | FDM operators + verification (manufactured solutions / convergence order) | ✅ done |
| M6 | Numerics: move ODE/time; add CG/BiCGSTAB/GMRES + nonlinear solvers | ✅ done |
| M7 | Coupling contracts, compatibility rules, configuration pipeline | ✅ done |
| M8 | Fidelity profiles, diagnostics, output rename | ✅ done |
| M9 | Presets: turbine/engine demotion, incompressible CFD, CHT | ✅ done |
| M10 | Docs/README/agent_guide/capability matrix update | ✅ done |
| M11 | Final: clean configure, full build, 100 tests green, conformance walkthrough | ✅ done |

