# extropian-physics

**Simulation physics library: data model, interfaces, and coupling framework.**

Contains the mesh data structures, field types, boundary condition framework, material database, solver plugin interface, and multiphysics coupling. Solver implementations live in separate `extropian-solver-*` repos.

Depends on `extropian-core` only.

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

## Building

```bash
cmake -S . -B build -G Ninja
cmake --build build
```

Requires: `extropian-core`.

## License

Business Source License 1.1 — see [LICENSE](LICENSE).
Converts to Apache 2.0 on 2029-05-26.
