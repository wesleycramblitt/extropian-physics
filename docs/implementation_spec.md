
# GPU-Portable Composable Multiphysics Simulation Engine
## Architecture & Implementation Specification — v0.3

---

# 1. Mission

Build a general-purpose, high-performance multiphysics simulation engine capable of solving a broad range of continuum and coupled physical problems while being designed for:

- CPU execution
- GPU execution
- future multi-GPU execution
- real-time and interactive simulation
- rapid low-fidelity prototyping
- high-fidelity engineering simulation
- novel/custom physics
- reusable numerical components
- composable multiphysics
- FDM, FVM, FEM and future discretizations
- PDEs and ODEs
- rigid bodies and particles
- strongly and weakly coupled systems
- structured and semi-structured meshes initially
- future extensibility to unstructured meshes

The engine must be a **composable numerical runtime**, not a collection of independent monolithic solvers.

The fundamental abstraction is:

```text
DATA
  +
OPERATORS
  +
PHYSICS
  +
DISCRETIZATION
  +
COUPLING
  +
EXECUTION
  +
HARDWARE
```

---

# 2. Product Scope

The engine owns:

- simulation state
- fields
- particles
- rigid bodies
- geometry
- mesh generation
- mesh topology
- physics models
- constitutive models
- boundary conditions
- initial conditions
- numerical operators
- discretization
- PDEs
- ODEs
- time integration
- linear solvers
- nonlinear solvers
- constraints
- coupling
- execution graphs
- CPU/GPU backends
- fidelity/performance control
- simulation presets
- basic diagnostics
- basic persistent output

The engine does NOT own:

- optimization
- inverse design
- machine-learning surrogates
- ROM systems
- advanced visualization
- advanced post-processing
- visualization applications
- optimization workflows

These may exist in separate repositories/products.

The engine must nevertheless expose efficient interfaces for those external systems.

---

# 3. Core Design Philosophy

## 3.1 Data-oriented

Numerical state must primarily consist of contiguous arrays.

Prefer:

```text
positions[N]
velocities[N]
temperature[N]
pressure[N]
density[N]
```

over:

```cpp
class Cell {
    ...
};
```

Avoid pointer-heavy object graphs in numerical hot paths.

---

## 3.2 Operators over objects

Numerical computation should be expressed as operators over fields/state.

Example:

```text
Gradient(T)
Divergence(u)
Laplacian(T)
Advection(u, T)
Diffusion(T, k)
PressureGradient(p)
```

rather than embedding numerical behavior inside mesh-cell objects.

---

## 3.3 Physics independent of hardware

Physics modules must not contain CUDA-specific or CPU-specific logic.

Physics declares mathematical requirements.

The execution system determines where and how those operations execute.

---

## 3.4 Physics and discretization are separate

Physics describes equations.

Discretization describes how equations become computable numerical operations.

Example:

```text
Navier-Stokes
      │
      ├── FDM
      └── FVM
```

and:

```text
Heat Equation
      │
      ├── FDM
      ├── FVM
      └── FEM
```

The engine must not duplicate physics implementations for each discretization.

---

# 4. Fundamental Architecture

```text
Simulation
    │
    ▼
Simulation State
    │
    ├── Fields
    ├── Particles
    ├── Rigid Bodies
    ├── Auxiliary State
    └── Solver State
    │
    ▼
Physics Modules
    │
    ▼
Numerical Operators
    │
    ▼
Discretization
    │
    ▼
Coupling Graph
    │
    ▼
Execution Graph
    │
    ├── CPU Backend
    └── GPU Backend
```

Surrounding infrastructure:

```text
Geometry
Mesh
Materials
Boundary Conditions
Initial Conditions
Constraints
Diagnostics
Fidelity Profiles
Presets
Output
```

---

# 5. Core Abstractions

## 5.1 State

State represents all evolving simulation data.

Examples:

```text
Fluid State
Thermal State
Structural State
EM State
Particle State
Rigid Body State
Reaction State
```

The complete simulation state may contain all of these simultaneously.

State must support:

- CPU residency
- GPU residency
- synchronization metadata
- field views
- checkpointing
- external access
- versioning/change tracking

---

# 6. Fields

A Field represents numerical data associated with a physical quantity.

Examples:

```text
Temperature
Pressure
Velocity
Density
Displacement
ElectricPotential
MagneticField
Concentration
Stress
Strain
```

A field must contain metadata such as:

```text
name
scalar/vector/tensor rank
components
units
location
association
precision
domain
mesh
boundary information
```

Possible locations:

```text
cell-centered
face-centered
node-centered
edge-centered
particle
rigid-body
quadrature
global
```

The numerical storage remains data-oriented.

---

# 7. Entity Sets

Support generalized collections of physical entities:

```text
Cells
Faces
Nodes
Edges
Particles
Rigid Bodies
Elements
DOFs
```

Each should be representable through contiguous/indexed data structures suitable for parallel execution.

---

# 8. Geometry

Geometry must be independent of numerical discretization.

Support:

- primitives
- CAD-like geometry interfaces
- implicit geometry
- analytical boundaries
- moving geometry
- parameterized geometry

Geometry must not assume FDM, FVM, or FEM.

---

# 9. Mesh

Initial mesh support:

```text
Cartesian
Structured
Block-structured
Semi-structured
```

Future:

```text
Unstructured
Adaptive
Octree
AMR
Hybrid
```

The API must not fundamentally assume a fixed Cartesian topology.

Separate:

```text
Geometry
Topology
Mesh Metrics
Discretization Metadata
```

Mesh should provide:

```text
cell volume
face area
face normal
centroid
neighbor relationships
boundary identifiers
node coordinates
connectivity
```

---

# 10. Mesh Generation

Provide automatic mesh generation for common cases.

Pipeline:

```text
Geometry
    ↓
Domain Analysis
    ↓
Mesh Generation
    ↓
Boundary Classification
    ↓
Mesh Validation
    ↓
Simulation Mesh
```

Initial capabilities:

- uniform refinement
- local refinement
- boundary refinement
- block decomposition
- resolution controls

Future:

- adaptive mesh refinement
- error-driven refinement
- automatic boundary-layer refinement
- geometry-aware refinement

---

# 11. Physics Modules

A Physics Module is a reusable physical component.

Examples:

```text
Fluid
Thermal
Structural
Electromagnetic
Reaction
Species Transport
Particle
Rigid Body
Radiation
Porous Media
```

A module should contain:

```text
State requirements
Fields
Equations
Constitutive models
Sources
Boundary requirements
Numerical requirements
Compatible discretizations
Coupling interfaces
Diagnostics
```

A module should NOT prescribe one immutable global solver configuration.

---

# 12. Module/Discretization Compatibility

Each module declares which discretizations it supports.

Example:

```text
Fluid:
    FDM: supported
    FVM: supported
    FEM: optional/future

Thermal:
    FDM: supported
    FVM: supported
    FEM: supported

Linear Elasticity:
    FEM: preferred
    FDM: supported
    FVM: optional

Electrostatics:
    FDM: supported
    FEM: supported
    FVM: supported
```

The engine must distinguish:

```text
Supported
Preferred
Default
Required
Experimental
Unsupported
```

A module may therefore declare:

```text
supported_discretizations = {FDM, FVM}
preferred_discretization = FVM
```

---

# 13. Discretization Selection

Discretization should be configurable **per module**, not globally forced across the simulation.

Example:

```text
Fluid       → FVM
Thermal     → FEM
RigidBody   → ODE
Controller  → ODE
```

This is essential for multiphysics.

However, not every combination is necessarily valid.

The engine must validate:

```text
Physics
    ×
Discretization
    ×
Mesh
    ×
Boundary Conditions
    ×
Coupling
    ×
Numerical Method
```

before execution.

---

# 14. Presets

Presets are preconfigured assemblies of reusable modules.

They are NOT separate solver implementations.

Examples:

```text
Incompressible CFD
Compressible CFD
Linear Static FEA
Modal Analysis
Electrostatics
Magnetostatics
Heat Transfer
Natural Convection
```

A preset defines:

```text
Modules
Default discretization
Default numerical methods
Default solver
Default BC conventions
Default timestep strategy
Default coupling strategy
Default fidelity profile
```

Users can override defaults where compatibility rules permit.

---

# 15. LEGO Architecture

The fundamental compositional model is:

```text
MODULES
   +
COUPLINGS
   +
CONSTRAINTS
   +
PRESETS
```

Modules are LEGO pieces.

Couplings connect them.

Presets are pre-built configurations.

Example:

```text
Fluid
Thermal
```

can become:

```text
Fluid
    │
    │ heat transfer
    ▼
Thermal
```

without requiring a dedicated:

```text
ConjugateHeatTransferSolver
```

Similarly:

```text
Fluid
+
RigidBody
```

creates fluid-structure interaction through coupling.

---

# 16. Coupling Engine

Coupling is a first-class subsystem.

The coupling engine represents relationships between modules.

Example:

```text
Fluid.velocity
    ↓
Thermal.convection
```

```text
Thermal.temperature
    ↓
Fluid.material_properties
```

```text
Fluid.force
    ↓
RigidBody.force
```

```text
RigidBody.velocity
    ↓
Fluid.boundary_velocity
```

```text
Electromagnetic.current
    ↓
Thermal.heat_source
```

---

# 17. Coupling Contracts

Every coupling must define a contract.

A coupling contract specifies:

```text
source
destination
quantity
units
field type
spatial association
temporal behavior
mapping
interpolation
conservation requirement
sign convention
execution frequency
coupling strength
```

Example:

```text
Source:
    Fluid.wall_shear

Destination:
    RigidBody.torque

Transformation:
    surface_integral

Units:
    Pa → N·m

Conservation:
    required

Execution:
    every timestep
```

---

# 18. Spatial Coupling Rules

Fields may live on different meshes.

The coupling engine must therefore support:

```text
same-mesh mapping
node → cell
cell → node
face → cell
surface → body
particle → field
field → particle
mesh → mesh
```

Mappings may use:

```text
nearest neighbor
linear interpolation
higher-order interpolation
conservative projection
L2 projection
weighted averaging
surface integration
volume integration
```

The coupling contract determines which is valid.

For conservation-critical transfers, arbitrary interpolation must not silently be allowed.

---

# 19. Temporal Coupling

Couplings must declare temporal behavior.

Examples:

```text
explicit
lagged
staggered
iterative
implicit
```

Example:

```text
Fluid(t)
    ↓
Force(t)
    ↓
RigidBody(t + dt)
```

versus:

```text
Fluid(t + dt)
    ↔
RigidBody(t + dt)
```

The latter requires iterative/implicit coupling.

---

# 20. Coupling Iteration

Strong coupling must support:

```text
Fluid
   ↓
Force
   ↓
Rigid Body
   ↓
Boundary Motion
   ↓
Fluid
```

repeated until:

```text
coupling_residual < tolerance
```

Parameters:

```text
max_iterations
tolerance
relaxation
under_relaxation
convergence_norm
```

---

# 21. Coupling Compatibility Rules

The engine must reject invalid combinations.

Examples:

```text
Temperature field
    → pressure field
```

should not be accepted without an explicitly defined transformation.

Similarly:

```text
velocity
    → temperature
```

requires a valid physical coupling such as advection/convection.

Rules must validate:

- quantity type
- dimensional units
- rank
- spatial association
- temporal compatibility
- physical meaning
- conservation requirements
- supported transformations

---

# 22. Coupling Constraint System

Create a general validation system:

```text
CompatibilityRule
```

Rules can inspect:

```text
Module
Discretization
Mesh
Field
Material
BoundaryCondition
Coupling
Solver
FidelityProfile
Backend
```

Validation occurs before execution.

Conceptually:

```text
Configuration
      ↓
Validation
      ↓
Compatibility Graph
      ↓
Execution Graph
```

Invalid configurations should produce explicit diagnostics.

Example:

```text
ERROR:
FEM thermal module uses nodal temperature.

FVM fluid module expects cell-centered temperature.

No valid coupling projection has been specified.

Suggested mappings:
    nodal → cell interpolation
    conservative projection
```

---

# 23. Physics Constraints

Modules may declare physical constraints.

Examples:

```text
mass conservation
momentum conservation
energy conservation
charge conservation
positivity
incompressibility
symmetry
boundary compatibility
```

The engine should distinguish:

```text
Hard Constraint
Soft Constraint
Diagnostic Constraint
```

---

# 24. Numerical Constraints

Numerical methods can declare requirements.

Examples:

```text
Explicit scheme:
    dt < CFL_limit

FEM:
    requires compatible element topology

FVM:
    requires face connectivity

Pressure projection:
    requires velocity-pressure coupling

Implicit solver:
    requires residual/Jacobian or Jacobian-vector interface
```

The engine should automatically validate these requirements.

---

# 25. Units and Dimensional Analysis

Physical quantities should carry units metadata.

The engine should support dimensional validation.

Example:

```text
force = mass × acceleration
```

must be dimensionally valid.

Invalid coupling:

```text
temperature → force
```

must require an explicit physical transformation.

Unit conversions should happen at interfaces, not inside performance-critical kernels where possible.

---

# 26. FDM

Support:

```text
gradient
divergence
curl
Laplacian
directional derivatives
stencils
ghost cells
boundary stencils
```

Optimize for:

- contiguous memory
- stencil locality
- kernel fusion
- GPU execution
- shared/local memory where useful

---

# 27. FVM

Represent:

```text
cells
faces
owner
neighbor
face area
face normal
cell volume
cell centroid
```

Residual:

```text
Residual =
    Face Fluxes
    +
    Sources
```

Support reusable flux schemes:

```text
Central
Upwind
Rusanov
HLL
HLLC
Roe
```

Prefer gather-based approaches where practical.

Avoid unnecessary atomics.

---

# 28. FEM

Represent:

```text
nodes
elements
DOFs
basis functions
quadrature
element connectivity
```

Support:

```text
element kernels
partial assembly
matrix-free operators
quadrature kernels
```

Global sparse matrix assembly must NOT be the fundamental abstraction.

Matrix-free execution is preferred.

---

# 29. ODE Framework

General form:

```text
dy/dt = f(t, y)
```

Initial integrators:

```text
Euler
RK2
RK3
RK4
```

Future:

```text
adaptive RK
BDF
Rosenbrock
implicit methods
symplectic methods
Newmark
generalized-alpha
```

ODE systems must use the same State/Operator infrastructure as PDE systems.

---

# 30. Rigid Bodies

Rigid bodies use SoA storage:

```text
position[N]
orientation[N]
linear_velocity[N]
angular_velocity[N]
mass[N]
inverse_mass[N]
inertia[N]
inverse_inertia[N]
force[N]
torque[N]
```

Integration should be highly parallel.

Collision/constraint solving may use:

```text
broad phase
narrow phase
spatial hashing
BVH
constraint graphs
graph coloring
iterative solvers
```

---

# 31. Particles

Support:

```text
position
velocity
mass
charge
temperature
species
```

Neighbor search:

```text
uniform grid
spatial hash
BVH
```

Particle ↔ field coupling must be supported.

Examples:

```text
particle → force field
particle → heat source
field → particle force
fluid → particle drag
```

---

# 32. Boundary Conditions

Boundary conditions are first-class objects.

Examples:

```text
Dirichlet
Neumann
Robin
Periodic
Wall
Slip
No-slip
Inlet
Outlet
Symmetry
Radiation
Flux
Convective
```

Boundary conditions must declare compatibility with:

- physics
- field
- discretization
- mesh location

---

# 33. Materials

Materials must be independent of discretization.

Examples:

```text
density
viscosity
conductivity
specific_heat
Young's modulus
Poisson ratio
permittivity
permeability
```

Support:

```text
constant
temperature-dependent
pressure-dependent
field-dependent
nonlinear
tabulated
user-defined
```

---

# 34. Numerical Operators

Operators are reusable mathematical transformations.

Examples:

```text
Gradient
Divergence
Curl
Laplacian
Advection
Diffusion
Interpolation
Projection
Integration
Reduction
Norm
Residual
Jacobian
JacobianVectorProduct
```

Operators must expose requirements and outputs.

---

# 35. Matrix-Free Numerical Architecture

Prefer:

```text
Input Fields
    ↓
Operator
    ↓
Output Field
```

instead of:

```text
Fields
    ↓
Global Matrix Assembly
    ↓
Sparse Matrix
    ↓
SpMV
```

where matrix-free formulation is practical.

Required interfaces should support:

```text
apply(x)
apply_transpose(x)
diagonal()
jacobian_vector_product(x, v)
```

where appropriate.

---

# 36. Linear Solvers

Initial:

```text
Jacobi
CG
BiCGSTAB
GMRES
```

Future:

```text
Multigrid
AMG
Block solvers
Domain decomposition
```

Solvers must operate on abstract operator interfaces.

---

# 37. Nonlinear Solvers

Initial:

```text
Fixed Point
Newton
Newton-Krylov
```

Support:

```text
residual evaluation
Jacobian
Jacobian-vector products
line search
damping
convergence criteria
```

---

# 38. Execution Graph

The execution graph represents computational dependencies.

Example:

```text
ComputeVelocity
       ↓
ComputeFlux
       ↓
ComputeResidual
       ↓
SolvePressure
       ↓
UpdateVelocity
       ↓
UpdateTemperature
```

The execution system handles:

- scheduling
- dependencies
- synchronization
- CPU/GPU placement
- asynchronous execution
- memory lifetime
- kernel fusion
- reductions
- transfers
- stream/event management

---

# 39. GPU Architecture

The GPU must be treated as a first-class execution target.

Prefer:

```text
Initialize
    ↓
Copy State → GPU
    ↓
Many Timesteps
    ↓
GPU-resident State
    ↓
Visualization / Analysis / Output
    ↓
Optional CPU synchronization
```

Do NOT do:

```text
GPU
 ↓
CPU
 ↓
GPU
```

every timestep.

---

# 40. Backend Abstraction

Backend interface should provide primitives such as:

```text
parallel_for
reduction
scan
sort
gather
scatter
stencil
synchronize
allocate
copy
async_copy
```

Initial:

```text
CPU
CUDA
```

Future:

```text
HIP
SYCL
Vulkan Compute
Metal
```

Physics must not depend on these backends.

---

# 41. GPU Context

Expose a GPU resource abstraction:

```text
GpuContext
GpuBuffer
GpuTexture
GpuMesh
GpuFieldView
GpuStream
GpuEvent
```

External rendering/post-processing systems must be able to consume GPU-resident simulation data.

Example:

```cpp
FieldView temperature =
    simulation.field("temperature");

renderer.render(temperature);
```

No CPU round-trip should be required.

---

# 42. Visualization/Post-Processing Boundary

Advanced visualization and post-processing may remain separate repositories.

However, they must share GPU resources where possible.

Architecture:

```text
Simulation Engine
       │
       ▼
GPU-resident State
       │
       ├── Renderer
       ├── Post Processing
       ├── Analysis
       └── Persistent Output
```

The simulation engine must provide a stable zero-copy/low-copy data access interface.

Derived quantities should preferably be computable directly on the GPU.

Example:

```text
Velocity
   ↓
Vorticity Operator
   ↓
Vorticity Field
   ↓
Renderer
```

No disk or CPU transfer required.

---

# 43. Persistent Output

Separate transient state from persistent output.

```text
Simulation State
      │
      ├── GPU Visualization
      ├── External Analysis
      ├── External Optimization
      └── File Output
```

Basic output formats:

```text
VTK / VTU
HDF5
XDMF
CSV
```

Output frequency must be configurable.

Real-time mode should minimize or disable persistent output.

---

# 44. External Consumer Contract

External systems should be able to consume:

```text
State
Fields
Mesh
Geometry
Observables
Diagnostics
Metadata
GPU resource handles/views
```

An external optimizer should ideally request:

```text
configuration
    ↓
simulation
    ↓
observable
```

rather than copying the complete simulation state.

---

# 45. Fidelity / Performance Profiles

Fidelity is a first-class engine concept.

The same simulation configuration should be capable of running at different fidelity levels.

Initial profiles:

```text
REALTIME
INTERACTIVE
FAST
BALANCED
ACCURATE
HIGH_FIDELITY
```

These are defaults, not hardcoded modes.

---

# 46. Fidelity Controls

A FidelityProfile may control:

```text
mesh resolution
temporal resolution
discretization order
floating-point precision
linear solver tolerance
nonlinear solver tolerance
maximum solver iterations
coupling iterations
coupling tolerance
physics model complexity
turbulence model
boundary-layer resolution
adaptive refinement
output frequency
```

The engine should permit fine-grained overrides.

Example:

```text
REALTIME:
    coarse mesh
    FP32
    few solver iterations
    loose tolerance
    simplified physics
    minimal output
```

versus:

```text
HIGH_FIDELITY:
    refined mesh
    FP64 where required
    strict convergence
    higher-order discretization
    more coupling iterations
    detailed physics
```

---

# 47. Runtime Performance Target

The engine should be able to answer:

```text
How fast can this simulation run?
```

and:

```text
How accurate is the current solution?
```

These are separate quantities.

Example:

```text
Simulation:
    18.7 ms / timestep
    53.5 timesteps/sec

Estimated accuracy:
    Temperature: ±2.8%
    Pressure: ±1.4%
```

---

# 48. Accuracy Estimation

Accuracy estimation is a first-class diagnostics subsystem.

Where mathematically appropriate, support:

```text
residual analysis
conservation error
mesh refinement comparison
temporal refinement comparison
Richardson extrapolation
local error indicators
discretization error estimation
solution convergence
constraint violation
physical invariant violation
```

The engine should distinguish:

```text
Numerical convergence
Discretization accuracy
Physical-model uncertainty
```

Do not represent these as the same metric.

---

# 49. Fidelity Escalation

The engine should eventually support:

```text
Run FAST
    ↓
Estimate uncertainty
    ↓
Identify insufficient regions
    ↓
Increase fidelity
    ↓
Run again
```

Potential future workflow:

```text
REALTIME
   ↓
INTERACTIVE
   ↓
BALANCED
   ↓
ACCURATE
```

This allows rapid prototyping followed by targeted refinement.

---

# 50. Diagnostics

Basic diagnostics belong in the engine.

Examples:

```text
residual
CFL number
mass conservation
energy conservation
momentum conservation
constraint violation
solver iterations
timestep
runtime
GPU utilization
memory consumption
```

Diagnostics should be accessible programmatically.

---

# 51. Recipes → Presets

Use the terminology:

```text
Preset
```

rather than Recipe for preconfigured engineering scenarios.

A Preset is an assembly of:

```text
Modules
Couplings
Discretizations
Solvers
Boundary conventions
Default numerical parameters
Fidelity defaults
```

---

# 52. Example Presets

## CFD

```text
Incompressible CFD
Compressible CFD
External Aerodynamics
Internal Flow
Heat Transfer CFD
Natural Convection
```

## FEA

```text
Linear Static
Linear Dynamic
Modal Analysis
Thermal Stress
Heat Conduction
Thermomechanical
```

## EM

```text
Electrostatics
Magnetostatics
Eddy Current
Joule Heating
```

## Thermal

```text
Heat Conduction
Transient Heat Transfer
Convection-Diffusion
```

Presets should compose modules rather than duplicate solver code.

---

# 53. Example: Generic User-Defined Simulation

A simulation engineer should be able to construct:

```text
Fluid Module
    └── Navier-Stokes
    └── FVM

Thermal Module
    └── Heat Equation
    └── FEM

Rigid Body Module
    └── Newton-Euler ODE

Coupling
    Fluid → Thermal
    Thermal → Fluid
    Fluid → RigidBody
    RigidBody → Fluid
```

The resulting system becomes:

```text
FVM Fluid
     ↕
FEM Thermal
     ↕
ODE Rigid Body
```

without requiring a special monolithic solver.

---

# 54. Configuration Pipeline

User configuration:

```text
Geometry
Materials
Modules
Discretizations
Boundary Conditions
Initial Conditions
Couplings
Numerical Settings
Fidelity
```

passes through:

```text
Parse
  ↓
Construct
  ↓
Validate
  ↓
Resolve Defaults
  ↓
Build Coupling Graph
  ↓
Build Execution Graph
  ↓
Allocate State
  ↓
Initialize Backend
  ↓
Execute
```

---

# 55. Validation System

Validation should occur before expensive simulation work.

Validation categories:

```text
Schema validation
Physical validation
Dimensional validation
Discretization validation
Mesh validation
Coupling validation
Boundary-condition validation
Solver validation
Fidelity validation
Backend validation
```

Diagnostics must explain:

```text
what is invalid
why it is invalid
which modules are involved
possible valid alternatives
```

---

# 56. Dependency Rules

Strict dependency direction:

```text
Core
 ↓
Numerics
 ↓
Discretization
 ↓
Physics
 ↓
Coupling
 ↓
Presets
```

Backend is orthogonal:

```text
Core
Numerics
Discretization
Physics
Coupling
      │
      └── Backend
```

Physics must never depend directly on CUDA.

FDM/FVM/FEM must not contain physics-specific duplication.

Presets must not implement numerical algorithms themselves.

---

# 57. Suggested Repository Structure

```text
engine/
│
├── core/
│   ├── state/
│   ├── fields/
│   ├── entities/
│   ├── memory/
│   ├── units/
│   ├── operators/
│   └── execution/
│
├── geometry/
│   ├── primitives/
│   ├── implicit/
│   └── interfaces/
│
├── mesh/
│   ├── topology/
│   ├── geometry/
│   ├── structured/
│   ├── block/
│   ├── generation/
│   └── validation/
│
├── discretization/
│   ├── fdm/
│   ├── fvm/
│   └── fem/
│
├── numerics/
│   ├── ode/
│   ├── time/
│   ├── linear/
│   ├── nonlinear/
│   ├── constraints/
│   ├── interpolation/
│   ├── projection/
│   └── multigrid/
│
├── physics/
│   ├── fluid/
│   ├── thermal/
│   ├── structural/
│   ├── electromagnetics/
│   ├── reaction/
│   ├── particles/
│   └── rigid_body/
│
├── coupling/
│   ├── graph/
│   ├── contracts/
│   ├── mapping/
│   ├── interpolation/
│   ├── conservation/
│   └── validation/
│
├── fidelity/
│   ├── profiles/
│   ├── estimation/
│   └── adaptation/
│
├── presets/
│   ├── cfd/
│   ├── fea/
│   ├── thermal/
│   ├── electromagnetics/
│   └── multiphysics/
│
├── backends/
│   ├── cpu/
│   └── cuda/
│
├── output/
│   ├── vtk/
│   ├── hdf5/
│   └── csv/
│
├── diagnostics/
│   ├── convergence/
│   ├── conservation/
│   ├── accuracy/
│   └── performance/
│
└── tests/
```

---

# 58. Fundamental Execution Primitives

The backend should eventually provide a small set of universal primitives:

```text
ParallelFor
Reduction
Scan
Sort
Gather
Scatter
Stencil
ElementKernel
Integration
Interpolation
SparseMatVec
ConstraintSolve
```

Higher-level numerical algorithms should compose these primitives.

---

# 59. Memory Architecture

Hot numerical paths must avoid:

```text
dynamic allocation
virtual dispatch
pointer chasing
unpredictable branching
unnecessary synchronization
host/device transfers
```

Fields should remain resident on the GPU across timesteps whenever practical.

Temporary memory should be managed through reusable arenas/pools.

---

# 60. Kernel Fusion

The execution graph should permit legal operations to be fused.

Example:

```text
Gradient(T)
    ↓
Multiply(k)
    ↓
Divergence(...)
```

may become a single GPU kernel where appropriate.

Fusion must be controlled by the execution layer, not manually duplicated in physics modules.

---

# 61. Asynchronous Execution

Support:

```text
GPU streams
events
dependency tracking
asynchronous transfers
overlapping communication/computation
```

Future multi-GPU support should be possible without redesigning physics modules.

---

# 62. Multi-GPU Preparation

Even before implementing multi-GPU, data structures should permit:

```text
domain decomposition
halo regions
owned data
ghost data
neighbor communication
```

Do not architect the system around a single monolithic global array that cannot later be partitioned.

---

# 63. Verification

Every numerical component requires verification.

Tests should include:

```text
analytical solutions
manufactured solutions
convergence tests
conservation tests
benchmark problems
cross-backend comparisons
CPU vs GPU
different fidelity levels
```

Example:

```text
FDM heat equation
    ↓
Known analytical solution
    ↓
Compute L2 error
    ↓
Refine mesh
    ↓
Verify expected convergence order
```

---

# 64. Preset Validation

Every preset must have validation cases.

For example:

```text
Incompressible CFD
    ├── lid-driven cavity
    ├── Poiseuille flow
    └── flow around cylinder

Linear FEA
    ├── cantilever beam
    └── patch test

Heat Transfer
    └── analytical conduction problem
```

---

# 65. Implementation Phases

## Phase 1 — Core Runtime

Implement:

```text
State
Field
EntitySet
Memory
Operators
Execution Graph
CPU Backend
```

---

## Phase 2 — Structured Mesh

Implement:

```text
Cartesian Mesh
Topology
Geometry
Boundary Identification
Mesh Generation
```

---

## Phase 3 — FDM

Implement:

```text
Stencil Operators
Boundary Conditions
Heat Equation
Diffusion
Advection-Diffusion
```

---

## Phase 4 — Basic Physics

Implement:

```text
Thermal
Fluid
Reaction-Diffusion
```

---

## Phase 5 — Numerical Solvers

Implement:

```text
ODE
Time Integrators
CG
GMRES
BiCGSTAB
Nonlinear Solvers
```

---

## Phase 6 — FVM

Implement:

```text
Face Topology
Fluxes
FVM Operators
Pressure-Velocity Coupling
Incompressible CFD
```

---

## Phase 7 — FEM

Implement:

```text
Elements
Basis Functions
Quadrature
DOFs
Element Kernels
Partial Assembly
Matrix-Free FEM
```

---

## Phase 8 — Bodies and Particles

Implement:

```text
Rigid Bodies
Particles
Neighbor Search
Constraints
```

---

## Phase 9 — Coupling Engine

Implement:

```text
Coupling Graph
Contracts
Field Mapping
Interpolation
Conservative Transfer
Temporal Coupling
Coupling Iteration
Validation
```

---

## Phase 10 — Presets

Implement:

```text
CFD
FEA
Thermal
EM
Multiphysics
```

using existing modules.

No duplicated solver implementations.

---

## Phase 11 — GPU

Implement CUDA backend.

Requirements:

```text
GPU-resident state
parallel kernels
reductions
stencils
FVM kernels
FEM kernels
matrix-free solvers
GPU execution graph
```

---

## Phase 12 — GPU-Aware External Interfaces

Implement:

```text
GpuContext
GpuBuffer
GpuFieldView
GpuMesh
GpuStream
GpuEvent
```

Allow external rendering/post-processing systems to consume GPU-resident state without mandatory CPU copies.

---

## Phase 13 — Fidelity and Accuracy

Implement:

```text
FidelityProfile
Performance Diagnostics
Residual Diagnostics
Conservation Diagnostics
Mesh Refinement Studies
Temporal Refinement
Error Estimation
```

---

# 66. Initial Physics Target

Prioritize:

```text
Heat Equation
Diffusion
Advection-Diffusion
Incompressible Navier-Stokes
Linear Elasticity
Reaction-Diffusion
Electrostatics
Rigid Body Dynamics
Particle Dynamics
```

These provide a broad test of:

```text
FDM
FVM
FEM
ODE
PDE
Coupling
GPU execution
```

---

# 67. Example Architecture

A generic simulation:

```text
Simulation
│
├── Fluid Module
│     ├── Navier-Stokes
│     ├── FVM
│     └── Material
│
├── Thermal Module
│     ├── Heat Equation
│     ├── FEM
│     └── Material
│
├── RigidBody Module
│     └── ODE
│
├── Coupling Graph
│     ├── Fluid ↔ Thermal
│     ├── Fluid ↔ RigidBody
│     └── Thermal → Material
│
├── Fidelity Profile
│     └── INTERACTIVE
│
└── Execution Graph
      └── CUDA
```

---

# 68. Preset Example

Selecting:

```text
Conjugate Heat Transfer
```

should internally produce something conceptually equivalent to:

```text
Fluid
    + Navier-Stokes
    + FVM
    + pressure solver

Thermal
    + Heat Equation
    + FVM

Coupling
    + fluid convection → thermal
    + thermal temperature → fluid properties
    + interface heat flux conservation

Numerics
    + timestep controller
    + linear solver

Fidelity
    + BALANCED
```

But the user can inspect and override the components.

---

# 69. Important Architectural Rule

Never create:

```text
FluidThermalSolver
FluidThermalRigidBodySolver
FluidThermalEMSolver
FluidThermalSpeciesSolver
```

as separate fundamental solver implementations.

Instead:

```text
Fluid
Thermal
RigidBody
EM
Species
```

are modules.

The coupling engine composes them.

Presets merely provide convenient validated configurations.

This prevents combinatorial explosion.

---

# 70. Generic vs Opinionated Usage

The engine must support both:

### Expert mode

```text
Create State
Create Physics
Select Operators
Select Discretization
Configure Coupling
Configure Solvers
Build Execution Graph
```

and:

### Preset mode

```text
Select:
    Incompressible CFD
```

then:

```text
Geometry
Materials
Boundary Conditions
Mesh
Fidelity
```

The same underlying infrastructure must execute both.

---

# 71. Ultimate Abstraction

The architecture should ultimately reduce to:

```text
                    SIMULATION
                         │
                         ▼
                       STATE
                         │
          ┌──────────────┼──────────────┐
          ▼              ▼              ▼
       PHYSICS       GEOMETRY        MATERIAL
          │              │              │
          ▼              ▼              │
      OPERATORS         MESH             │
          │              │              │
          └──────────────┼──────────────┘
                         ▼
                   DISCRETIZATION
                         │
                         ▼
                    COUPLING
                         │
                    VALIDATION
                         │
                  FIDELITY PROFILE
                         │
                  EXECUTION GRAPH
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
             CPU                   GPU
                                    │
                    ┌───────────────┼───────────────┐
                    ▼               ▼               ▼
                Renderer       Postprocess       Output
```

---

# 72. Non-Negotiable Design Rules

1. **Fields are arrays.**
2. **Operators operate on fields/state.**
3. **Physics is independent of hardware.**
4. **Discretization is independent of physics.**
5. **Discretization is selectable per module.**
6. **Presets provide defaults; they do not create duplicate solvers.**
7. **Modules declare compatibility; the validation system enforces it.**
8. **Couplings are explicit contracts, not arbitrary field assignments.**
9. **Units and dimensional consistency are validated.**
10. **Spatial and temporal coupling rules are explicit.**
11. **Conservation requirements are explicit.**
12. **Invalid module/discretization/coupling combinations fail before execution.**
13. **Matrix-free execution is preferred where practical.**
14. **Hot numerical paths contain no unnecessary dynamic allocation or virtual dispatch.**
15. **GPU state remains GPU-resident during simulation.**
16. **External rendering/post-processing can consume GPU state without mandatory CPU copies.**
17. **Real-time fidelity is a first-class execution mode, not a separate solver.**
18. **Accuracy/error estimation is separate from numerical convergence.**
19. **Multi-GPU/domain decomposition must remain architecturally possible.**
20. **New physics should be added by composing existing infrastructure whenever possible.**
21. **New multiphysics combinations should normally require new coupling configuration, not a new monolithic solver.**
22. **External optimization and advanced post-processing must be able to consume the engine without modifying its numerical core.**

---

# 73. Definition of Success

The architecture is successful if:

### Adding a new physical model

primarily requires:

```text
State
Equations
Operators
Constitutive Models
BCs
Compatibility Rules
```

rather than a new execution system.

### Adding a discretization

does not require rewriting physics.

### Adding a GPU backend

does not require rewriting physics or presets.

### Adding a new multiphysics scenario

primarily requires:

```text
Modules
+
Coupling Contracts
```

rather than a monolithic solver.

### Supporting a common engineering problem

primarily requires:

```text
Preset
+
Geometry
+
Materials
+
BCs
+
Mesh
+
Fidelity
```

### Supporting a novel research problem

allows the simulation engineer to manually construct:

```text
Modules
+
Discretizations
+
Operators
+
Couplings
+
Constraints
+
Solvers
```

### Running in real time

allows the same simulation to reduce fidelity without fundamentally changing the model architecture.

### Rendering in real time

allows:

```text
GPU Simulation State
        ↓
GPU Derived Fields
        ↓
GPU Renderer
```

without unnecessary:

```text
GPU → CPU → GPU
```

transfers.

---

# 74. Core Product Concept

The engine should ultimately feel like a **numerical LEGO system**.

Users should be able to think:

```text
I need fluid.
I need heat.
I need a rigid body.
Connect them.
```

while the expert can think:

```text
Navier-Stokes
+
FVM
+
custom flux
+
thermal PDE
+
FEM
+
custom constitutive law
+
strong coupling
+
custom constraint
+
GPU execution
```

Both workflows must resolve into the same underlying runtime.

The central architectural principle is therefore:

> **Build one composable numerical execution engine, then assemble physical systems from validated modules, discretizations, operators, constraints, and coupling contracts.**

Presets make common engineering problems easy.

The generic module system makes novel physics possible.

Fidelity profiles make the same engine useful from real-time interactive prototyping through high-fidelity simulation.

The GPU-resident state interface allows simulation, rendering, and external analysis to operate as one high-performance computational pipeline rather than independent applications.

