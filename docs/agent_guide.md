# Extropian Physics — Agent Guide

This document is the authoritative reference for anyone (human or AI) working on
this codebase. Follow these rules strictly. When in doubt, the code in
`src/fluid/reduced_order/bem/` is the reference implementation.

---

## 1. Project Identity

**extropian-physics** is a C++23 static library providing the physics
infrastructure layer for the Extropian simulation ecosystem. It owns:

- Mesh representation (nodes, cells, boundary faces)
- Fields (scalar, vector, tensor) on meshes
- Boundary condition framework
- Material database
- Solver plugin interface (lifecycle, stepping, coupling surfaces)
- Solver management (registration, creation, lifecycle)
- Coupling between solvers (surface mapping, data exchange)
- Domain-specific solvers that implement the plugin interface

**Not owned by this repo:** solver implementations. High-fidelity solvers
(CFD, FEA, etc.) live in separate `extropian-solver-*` repos and implement
`ISolverPlugin`. The BEM turbine solver is the one exception — it lives here
because it's a reduced-order model tightly coupled to the physics data model.

### External dependencies

| Dependency | Purpose | Source |
|---|---|---|
| `extropian-core` | Pool allocators, handles, hash, serialization, config, logging, ECS | FetchContent |
| `extropian-geometry` | `TurbineDefinition`, `BladeRow`, `FlowPath`, splines | FetchContent |
| `doctest` | Test framework (tests only) | FetchContent |
| `nlohmann_json` | JSON parsing (transitive via core) | Transitive |
| `Threads` | Solver threading (linked, not yet used in BEM) | System |

---

## 2. Repository Layout

```
extropian-physics/
├── include/exd/physics/          PUBLIC API HEADERS — contracts only
│   ├── physics.hpp               Umbrella header
│   ├── solver/                   Plugin interface + manager
│   └── fluid/reduced_order/bem/  BEM public types (config, result, airfoil, solver)
├── src/                          PRIVATE IMPLEMENTATION
│   ├── mesh/                     Mesh types + I/O stubs
│   ├── field/                    Field types + interpolation
│   ├── bc/                       Boundary condition framework
│   ├── material/                 Material database
│   ├── solver/                   Solver manager, time stepping, convergence
│   ├── coupling/                 Coupling manager + surface mapping
│   ├── io/                       Output channels (field stamps, CSV series)
│   ├── fluid/fdm3/               3D FDM solver (the coupled-CFD reference)
│   ├── engine/                   Slider-crank engines (Otto + steam)
│   └── fluid/reduced_order/bem/  BEM solver (the reference implementation)
├── tests/unit/                   doctest-based unit tests
├── docs/                         Architecture docs + this agent guide
└── data/airfoils/               Synthetic NACA polar CSVs
```

### The Public/Private header split

This is the most important structural rule in the repo:

- **`include/exd/physics/`** — Public API. Contains only type definitions,
  interfaces, and the solver entry point. Consumers of this library include
  only these headers. Everything here must be stable, well-documented, and
  free of implementation details.

- **`src/`** — Private implementation. The include path is `PRIVATE` in
  CMake so downstream projects cannot see these headers. This is where
  internal types, strategy interfaces, factory functions, and concrete
  models live.

The one exception: `bem_internal.hpp` is a private header that defines the
BEM solver's internal interfaces (`InductionModel`, `LossModel`,
`BladeGeometry`, etc.). Tests access it via a relative include when they
need to test internal components directly.

**Rule:** When adding a new subsystem, never put implementation types in
`include/`. Put them in `src/<subsystem>/` and keep the public surface
minimal.

---

## 3. Naming Conventions

These are enforced, not suggestions.

### Namespaces

Deeply nested, matching the directory path:

```
exd::physics::mesh
exd::physics::field
exd::physics::bc
exd::physics::material
exd::physics::solver
exd::physics::coupling
exd::physics::fluid::reduced_order::bem
```

Future domains: `exd::physics::thermal`, `exd::physics::structural`,
`exd::physics::electromagnetics`.

Reference docs per module: `fdm3_architecture.md`, `engine_architecture.md`,
`turbine_coupling_architecture.md`, `io_architecture.md`,
`output_channels.md`, `capability_matrix.md`.

**Rule:** Namespace nesting mirrors directory nesting. A file at
`src/foo/bar/baz.cpp` should use namespace `exd::physics::foo::bar`.

### File names

`snake_case` for everything:
- Headers: `bem_solver.hpp`, `plugin_interface.hpp`
- Sources: `bem_solver.cpp`, `induction_corrections.cpp`
- Tests: `induction_corrections_test.cpp`

### Types

`PascalCase` for all types: structs, classes, enums, type aliases.
Enum class values are also PascalCase: `InductionCorrection::GlauertIterative`.

### Functions

`snake_case` for all functions:
- Public: `solve_turbine()`, `max_shroud_radius_over_domain()`
- Factories: `make_prandtl_loss_model()`, `make_induction_model()`
- Internal: `build_blade_geometry()`, `compute_duct_state()`

### Members

`snake_case` for member variables. Private storage uses trailing underscore:
`materials_`, `storage_`. Public members (on value-type results) have no
trailing underscore: `axial_velocity`, `induction_axial`.

### Constants

`UPPER_SNAKE_CASE` for `constexpr` constants: `SOLVER_CREATE_SYMBOL`.
Also acceptable: `constexpr double` in PascalCase scope.

---

## 4. Architecture Patterns

### 4.1 Strategy + Factory (the primary pattern)

This is how every pluggable component works. The BEM solver is the
canonical example.

**Step 1: Define an abstract interface** (in a private header):
```cpp
// src/fluid/reduced_order/bem/bem_internal.hpp
class InductionModel {
public:
    virtual ~InductionModel() = default;
    virtual void solve(const BladeElementInput& elem, ...) const = 0;
};
```

**Step 2: Implement concrete models** in anonymous namespaces within `.cpp`
files:
```cpp
// src/fluid/reduced_order/bem/induction_corrections.cpp
namespace {
class GlauertIterativeInductionModel final : public InductionModel {
public:
    void solve(...) const override { ... }
};
} // anonymous namespace
```

**Step 3: Write a named factory function** in the BEM namespace:
```cpp
std::unique_ptr<InductionModel> make_glauert_iterative_induction_model();
```

**Step 4: Write a dispatcher factory** that switches on a config enum:
```cpp
std::unique_ptr<InductionModel> make_induction_model(InductionCorrection type) {
    switch (type) {
        case InductionCorrection::Standard:           return make_standard_induction_model();
        case InductionCorrection::GlauertIterative:   return make_glauert_iterative_induction_model();
        case InductionCorrection::Snel:               return make_snel_induction_model();
    }
}
```

**Step 5: The solver calls the dispatcher once**, then uses the polymorphic
interface in the hot loop:
```cpp
auto induction = make_induction_model(config.induction_correction);
auto loss = make_loss_model(config.loss_correction);
for (auto& elem : geometry.elements) {
    induction->solve(elem, ..., loss, ...);
}
```

**Why this pattern:**
- Adding a new model = new `.cpp` file + one line in the dispatcher.
- Models are testable in isolation via their interface.
- No hot-loop branching (vtable dispatch is fine for BEM; for CFD kernels,
  consider CRTP or compile-time dispatch — see §6).
- Anonymous namespace hides implementation details; only the factory is
  visible in the translation unit.

### 4.2 Value-type results (no exceptions)

All result types are plain structs with no exceptions. Error reporting is
done via fields on the result:

```cpp
struct TurbineResult {
    bool valid = false;
    std::string error;
    std::vector<std::string> warnings;
    // ... result data ...
};
```

**Rule:** Never throw exceptions in solver code. Use `valid` + `error` +
`warnings`. Callers check `valid` before accessing result data.

### 4.3 Config via struct (not maps)

Configuration is a typed struct with sensible defaults, not a string map:

```cpp
struct BEMSolverConfig {
    int element_count = 20;
    double k_duct = 0.0;
    InductionCorrection induction_correction = InductionCorrection::Standard;
    LossCorrection loss_correction = LossCorrection::Prandtl;
    // ...
};
```

**Rule:** Never use `ConfigNode` (string map) for internal solver config.
`ConfigNode` exists only for the `ISolverPlugin` interface where
cross-language flexibility is needed. Internal solvers use typed structs.

### 4.4 Pimpl for library boundaries

`SolverManager` uses the pimpl pattern because it's in the public API and
we want to hide implementation details from consumers:

```cpp
// Public header
class SolverManager {
    struct Impl;          // forward-declared
    std::unique_ptr<Impl> impl_;
public:
    static SolverManager& instance();
    // ...
};
```

**Rule:** Use pimpl only when the type is in the public API and you need
to hide implementation details. Internal types (like BEM models) don't
need pimpl.

### 4.5 Domain taxonomy

Organize code by physics domain:

```
src/mesh/          Mesh, I/O, partitioning, quality
src/field/         Scalar, vector, tensor fields, interpolation
src/bc/            Boundary condition framework
src/material/      Material database
src/solver/        Solver manager, time stepping, convergence
src/coupling/      Coupling manager, surface mapping
src/<domain>/      Domain-specific solvers
```

The `<domain>` directory contains solver implementations that implement
`ISolverPlugin` or, like the BEM solver, are reduced-order models that
use the physics data model directly.

---

## 5. Build System

CMake 3.21+, Ninja generator. C++23 required.

### Build commands

```bash
cmake -S . -B build -DEXT_PHYSICS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

### Local dependency overrides

```bash
cmake -S . -B build \
  -DEXT_CORE_DIR=/path/to/extropian-core \
  -DEXT_GEOMETRY_DIR=/path/to/extropian-geometry \
  -DEXT_PHYSICS_BUILD_TESTS=ON
```

### Adding a new source file

Add it to the `add_library(exd-physics STATIC ...)` block in the root
`CMakeLists.txt` under the appropriate section comment. Group files by
subsystem with section headers:

```cmake
    # Level-3 BEM turbine solver
    src/fluid/reduced_order/bem/bem_solver.cpp
    ...
```

### Adding a new test

In `tests/CMakeLists.txt`, use the `add_physics_test()` macro:

```cmake
add_physics_test(new_test_name unit/bem/new_test.cpp)
```

The macro links against `exd::physics` and `doctest::doctest_with_main`,
and sets `EXT_PHYSICS_DATA_DIR` for file-based tests.

---

## 6. Solver Evolution Strategy

Solvers progress through three tiers. Each tier is a superset of the
previous one, not a replacement.

### Tier 1: CPU single-threaded

Current state. All solvers are sequential. The BEM solver is a Tier 1
solver. Optimizations at this tier are algorithmic (convergence,
numerical stability) and structural (modularity, testability).

### Tier 2: CPU multi-threaded

Thread-parallel execution. The `Threads::Threads` dependency is already
linked. Strategies:
- **Element-level parallelism** (BEM): each blade element is independent
  after geometry setup. Use `std::for_each` or `std::async` over elements.
- **Block-level parallelism** (CFD): partition mesh into blocks, run
  block solvers in parallel with halo exchanges.

**Rule for Tier 2:** Thread parallelism must not change the public API.
The `solve_turbine()` function signature stays the same. Parallelism is
an internal implementation detail controlled by config:

```cpp
struct BEMSolverConfig {
    // ...
    int thread_count = 0;  // 0 = auto, 1 = single-threaded
};
```

### Tier 3: GPU compute

GPU-accelerated solvers. Strategies:
- **CUDA/HIP kernels** for element-level loops (BEM, FEM element
  assembly).
- **Batched execution** — keep data on GPU between timesteps. The
  coupling surface API already uses `std::span` which maps to device
  pointers.
- **Zero-copy where possible** — if the optimizer (extropian-optimization)
  calls the solver in a batch, the entire batch can live on GPU.

**Rule for Tier 3:** Design GPU-ready data layouts now. Specifically:
- Flat arrays-of-structures (AoS) or structures-of-arrays (SoA) — not
  nested `std::vector<struct>`.
- All result types should be trivially copyable to device memory.
- Avoid `std::string` in hot-path data types (use enums or fixed-size
  char arrays).

### When to add a new tier

Don't add GPU code until the CPU version is correct and tested. Don't add
threading until the sequential version is architecturally clean. The
modularity built in Tier 1 makes Tier 2 and 3 feasible without rewriting.

---

## 7. Modularity Principles

### 7.1 One file, one responsibility

Each `.cpp` file should implement one cohesive concept:
- `induction.cpp` → `StandardInductionModel`
- `induction_corrections.cpp` → `GlauertIterativeInductionModel`, `SnelInductionModel`
- `losses.cpp` → `PrandtlLossModel`, loss model dispatcher
- `loss_corrections.cpp` → `DuSeligLossModel`, `ChaviaropoulosLossModel`

**Rule:** If a file exceeds ~500 lines, split it. If it has two unrelated
concepts, split it.

### 7.2 Private headers for internal interfaces

Interfaces that are not part of the public API go in `src/`, not
`include/`. The file `bem_internal.hpp` is the model:

```cpp
// src/<subsystem>/<subsystem>_internal.hpp
// Contains: interfaces, internal types, factory declarations
// NOT included by any public header
```

**Rule:** Name it `<subsystem>_internal.hpp`. Never include it from
`include/exd/physics/`.

### 7.3 Factory functions over direct construction

Never construct model classes directly in solver code. Always go through
a factory function:

```cpp
// WRONG
GlauertIterativeInductionModel model;

// RIGHT
auto model = make_glauert_iterative_induction_model();
```

This ensures:
- New models can be added without touching solver code.
- Tests can inject mock models.
- The dispatcher pattern works cleanly.

### 7.4 No circular dependencies

The dependency graph must be a DAG:

```
public headers (include/)
    ↓
internal headers (src/*/)_internal.hpp
    ↓
implementation files (src/*//*.cpp)
```

Never have an internal header depend on a `.cpp` file. Never have a
public header depend on an internal header.

### 7.5 Coupling stubs

When adding a new solver, stub out coupling support from day one. The
`CouplingSurface` interface in `plugin_interface.hpp` defines what
coupling data a solver can provide. Implement it minimally:

```cpp
CouplingSurface get_coupling_surface(const std::string& name) override {
    // Return empty/default surface. Implement when coupling is needed.
    return {};
}
```

**Rule:** Every solver must implement `ISolverPlugin` and provide a
`get_coupling_surface()`, even if it returns empty data. This ensures the
coupling manager can discover and connect solvers without code changes.

---

## 8. Configuration Architecture

### 8.1 Config struct pattern

Every solver defines a config struct in the public API:

```cpp
// include/exd/physics/<domain>/<solver>_config.hpp
struct SolverConfig {
    // Time stepping
    double dt = 0.01;
    int max_steps = 1000;
    
    // Meshing
    int element_count = 20;
    
    // Output control
    bool output_flow_field = false;
    bool output_per_element = true;
    
    // Model selection (enums)
    SomeModel some_model = SomeModel::Default;
    
    // Future-proofing: extensible parameter bag
    // (avoid using this — prefer typed fields)
    // std::unordered_map<std::string, double> extra;
};
```

### 8.2 What to make configurable

Every solver should expose config for at minimum:

| Category | Examples |
|---|---|
| Time stepping | `dt`, `max_steps`, `convergence_tolerance` |
| Meshing | `element_count`, `mesh_resolution`, `refinement_level` |
| Output | `output_flow_field`, `output_per_element`, `output_filename` |
| Model selection | Enum fields for pluggable algorithms |
| Coupling | `coupling_interval`, `coupling_surface_names` |
| Threading | `thread_count` (0 = auto) |

### 8.3 Enum-based model selection

Always use `enum class` for model selection, never strings:

```cpp
// GOOD
enum class InductionCorrection { Standard, GlauertIterative, Snel };

// BAD
std::string induction_correction = "standard";
```

Enums give compile-time safety, enable `switch` statements (compiler
warns on missing cases), and are self-documenting.

---

## 9. Testing Standards

### 9.1 Framework

doctest v2.4.11. All tests link against `doctest::doctest_with_main`.

### 9.2 Test file structure

```cpp
#include <doctest/doctest.h>
#include <exd/physics/fluid/reduced_order/bem/bem_solver.hpp>

using namespace exd::physics::fluid::reduced_order::bem;

namespace {
// Shared fixtures in anonymous namespace
TurbineDefinition make_turbine() { ... }
BEMSolverConfig base_config() { ... }
} // anonymous namespace

TEST_CASE("descriptive name of what is tested") {
    // Arrange
    auto turbine = make_turbine();
    auto config = base_config();
    
    // Act
    auto result = solve_turbine(turbine, conditions, polars, config);
    
    // Assert
    REQUIRE(result.valid);
    CHECK(result.rotor.cp > 0.0);
    CHECK(result.rotor.cp < 0.6);
}
```

### 9.3 Naming conventions

- **File:** `<thing>_test.cpp`
- **TEST_CASE:** "descriptive name" — what behavior is being verified
- **Sections:** `SUBCASE("specific scenario")` for variations

### 9.4 Fixture duplication

Currently, test fixtures (`cylindrical_shroud()`, `rotor_betz()`, etc.)
are duplicated across test files. This is intentional for now — each test
file is self-contained. When the fixture set grows, extract a shared
`tests/unit/bem/test_fixtures.hpp` header. Until then, duplication is
acceptable to avoid coupling test files.

### 9.5 Floating point comparisons

Always use `doctest::Approx` with explicit tolerance:

```cpp
CHECK(result.rotor.cp == doctest::Approx(0.52).epsilon(0.05));
```

Never use exact equality for floating point results.

### 9.6 Testing internal components

When testing internal types (like specific induction models), include the
private header via relative path:

```cpp
#include "../../../../src/fluid/reduced_order/bem/bem_internal.hpp"
```

This is the only acceptable use of relative includes in tests.

---

## 10. Code Style

### 10.1 Formatting

- 4 spaces indentation (no tabs)
- Opening brace on same line as statement
- `#pragma once` for include guards
- ~100 character line limit (soft)

### 10.2 Includes

Order:
1. Corresponding header (if .cpp file)
2. Project headers (`<exd/physics/...>`)
3. Library headers (`<doctest/...>`)
4. Standard library headers

Use quotes for project headers, angle brackets for everything else.

### 10.3 Comments

- File-level comment explaining purpose
- `// ── Section Name ──` for major sections within a file
- Doc comments on public APIs (doxygen-style `///` or `/** */`)
- No commented-out code — delete it; git has history

### 10.4 Error handling

- No exceptions in solver code
- Return `valid`/`error`/`warnings` on result types
- `std::printf` for logging (no logging framework dependency)
- Warnings are non-fatal (e.g., "convergence not reached in 100
  iterations"); errors are fatal (e.g., "invalid config: rpm must be
  positive")

---

## 11. Future-Proofing Rules

### 11.1 GPU data readiness

All result structs must be trivially copyable. Avoid:
- `std::string` in hot-path result types (use enums, fixed-size arrays)
- `std::vector` in per-element results (use flat arrays with span views)
- Virtual methods in result types

### 11.2 Coupling stubs

Every new solver must implement `get_coupling_surface()` returning a
default/empty surface. This is zero-cost and ensures future coupling
works without refactoring.

### 11.3 Extensible config

Use typed enum fields for model selection. Avoid
`std::unordered_map<std::string, double>` config bags — they're
untestable and error-prone. If you need a new knob, add a typed field.

### 11.4 Optimizer integration

Solvers should be callable in batch mode. The `solve_turbine()` pattern
(function takes definition + conditions + config, returns result) is
optimizer-friendly. New solvers should follow the same pattern:

```cpp
Result solve_thing(const Definition& def, const Conditions& cond, const Config& config);
```

No global state, no side effects beyond the returned result.

### 11.5 Time stepping

The `TimeStepper` in `src/solver/time_stepping.cpp` provides CFL-based
adaptive dt. New solvers should accept a `TimeStepper` config in their
config struct and call `stepper.advance()` in their step loop, rather
than implementing their own dt logic.

---

## 12. Checklist for New Solver Implementation

Before writing any code, ensure:

- [ ] Public config struct in `include/exd/physics/<domain>/`
- [ ] Public result struct in `include/exd/physics/<domain>/`
- [ ] Public entry point function in `include/exd/physics/<domain>/`
- [ ] Private internal header in `src/<domain>/`
- [ ] Implementation files in `src/<domain>/`
- [ ] Factory functions for pluggable models
- [ ] Dispatcher factory on config enums
- [ ] Coupling surface stub
- [ ] Config for: time stepping, meshing, output, threading, model selection
- [ ] Unit tests in `tests/unit/<domain>/`
- [ ] Entry in root `CMakeLists.txt` source list
- [ ] Entry in `tests/CMakeLists.txt`
- [ ] Section in `docs/` architecture doc
- [ ] README update with usage example

---

## 13. Common Pitfalls

1. **Adding to public headers too early.** Move types to `include/` only
   when a solver outside this repo needs them. Everything else stays in
   `src/`.

2. **Using `std::string` in result types.** Use enums or fixed-size arrays
   in per-element results. Save strings for top-level error messages only.

3. **Skipping the factory pattern.** Direct construction of models in
   solver code makes the system rigid and untestable. Always use factories.

4. **Hard-coding mesh resolution.** Always make meshing configurable. The
   BEM solver uses `element_count`; CFD solvers should use
   `mesh_resolution` or similar.

5. **Ignoring coupling.** Even if your solver doesn't need coupling today,
   implement `get_coupling_surface()`. The coupling manager needs it to
   discover your solver's available fields.

6. **Throwing exceptions.** Solver code must never throw. Use the
   `valid`/`error`/`warnings` pattern on result types.

7. **Using relative includes in production code.** The one exception is
   tests accessing `bem_internal.hpp`. Production code should only use
   paths relative to `include/` or `src/`.

8. **Committing build artifacts.** The `build/` directory is committed for
   convenience but should not be treated as canonical. Clean builds should
   always work from source.

---

## 14. References

- `docs/bem_level3_architecture.md` — BEM solver implementation contract
- `docs/modular_solver_architecture.md` — Whole-system architecture + phased
  implementation plan (integrators, rigid bodies, engines, compressors,
  electrical/EM, coupling, FVM/FEM/mesh roadmap)
- `docs/agent_guide.md` — This document
- `include/exd/physics/solver/plugin_interface.hpp` — ISolverPlugin contract
- `src/fluid/reduced_order/bem/bem_internal.hpp` — BEM internal interfaces
- `src/fluid/reduced_order/bem/bem_solver.cpp` — BEM orchestration (reference impl)
