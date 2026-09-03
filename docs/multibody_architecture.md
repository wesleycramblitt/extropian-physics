# Multibody Architecture — Boundary, Model Contract, and Roadmap

**Status: PLAN (not implemented).** Authoritative statement of (a) the
external-geometry-repo boundary, (b) the format-neutral model-description
contract the engine owns, and (c) the roadmap for a general multibody solver
(arbitrary joint types and constraints). Companion to
`implementation_spec.md` (§30 rigid bodies, §16–22 coupling, §54 pipeline)
and `implementation_spec_conformance.md`.

---

## 1. The boundary (external geometry repo)

```
CAD / MuJoCo / any format        ← EXTERNAL REPO (owned elsewhere)
   │   geometry authoring
   │   mesh loaders (STL/OBJ/GLTF/MJCF/URDF/STEP…)
   │   converter + exporter
   ▼
ENGINE MODEL DESCRIPTIONS        ← THIS REPO owns the schemas (below)
   │   format-neutral, versioned, validation-ready
   ▼
extropian-physics engine         ← THIS REPO
   │   validate → build → allocate → execute (spec §54)
   ▼
solutions / observables
```

**Rules:**
1. The engine does **not** parse CAD formats. No MJCF/URDF/STL/OBJ parsers,
   no mesh loaders, no geometry kernels. All of that lives in the external
   repo ("unless absolutely necessary here" — it is not: the engine's own
   contract is the description, not the file formats).
2. The engine **owns the description schemas** — the data contracts the
   external repo emits. Schema versioning + validation live here (the
   `coupling::rules` compatibility machinery extends to model descriptions).
3. The external repo may target the engine through its **C++ API directly**
   (fill config structs) or through a **serialized model description** the
   engine reads with its own, versioned reader (a simple owned format; never
   a third-party one).
4. "Import CAD → boom the solver runs" therefore means: the external repo
   converts geometry/BCs/parameters into the engine's configs, then calls
   the class-level presets (`steady_conduction`, `manipulator_trajectory`,
   `incompressible`, …). The engine side already accepts data-driven
   configs everywhere (per-node fields, region lists, per-joint vectors).

**Existing contract points the external repo fills today:**
| Domain | Config to fill | Geometry source |
|---|---|---|
| Heat conduction | `physics::thermal::HeterogeneousConductionConfig` (per-node k/q, regions, face BCs) | CAD bodies → regions or per-node fields |
| Incompressible flow | `physics::fluid::fdm3::FDM3Config` (grid + face BCs) | CAD extents → grid, faces → BC kinds |
| Elastic structure | `ElasticityConfig` (grid, pins, tractions, thermal channel) | CAD → grid + restrained regions |
| Manipulator (serial, planar) | `physics::robotics::SerialManipulatorConfig` (links, joints, stops, gains) | CAD bodies → link lengths/masses; joints → limits |
| Multibody (general) | `MultibodyDescription` (planned, §3 below) | MuJoCo/URDF tree → bodies/joints/constraints |

---

## 2. Requirements for "any joint type and constraint"

The manipulator generalization must cover:

**Joint types (1 axis or more per joint):**
- revolute (hinge, 1R), prismatic (slide, 1T), screw (1 DOF, pitch-coupled
  rotation/translation), cylindrical (1R+1T), universal (2R), planar (2T+1R),
  ball (3R, quaternion), free (6 DOF), weld/fixed (0 DOF).

**Constraints beyond joints:**
- loop-closure joints (e.g., a closed four-bar: the tree joints + one
  closure constraint), distance constraints, equal/prescribed coordinate
  constraints, later contacts (collision).

**Robot-relevant features:**
- actuated joints with limits, damping, stiffness, friction (soft limits
  like the current stops, generalized);
- actuators: torque/position/velocity servos (beyond the current PD);
- loads: gravity, end-effector wrenches, motors on joints;
- initial state, reference trajectories.

---

## 3. The model description (the engine-owned contract)

**Authoritative schema: `docs/multibody_description_v0.1.md`** — v0.1 is
STABLE (2026-09): conventions (SI, Z-up, wxyz quaternions, parent-frame
anchoring), the exact struct layouts, the joint DOF table, validation rules
V1–V16, JSON mapping, MJCF/current-config correspondence, and the versioning
policy.  Summary of the shape:

```text
MultibodyDescription
├── version
├── gravity, integration hints (dt, method)
├── bodies[]                      // kinematic tree via parent links
│   ├── name, parent ("" = root)
│   ├── joint: { type, axis, anchor (in parent frame),
│   │            limits[], damping, stiffness, friction,
│   │            actuated: bool }
│   ├── mass, center_of_mass, inertia (principal + axes, or 3×3)
│   └── visual/geometry mesh name (metadata only — rendering lives outside)
├── constraints[]                 // non-tree constraints
│   ├── type: loop | distance | prescribed …
│   ├── bodyA, bodyB, frames (in-body anchors), limits
├── actuators[]
│   ├── type: torque | position | velocity (servo gains, limits)
│   └── target_joint
├── loads[]                       // applied wrenches, point forces
├── initial_state                // q, dq per joint (or per body)
└── controls[]                    // controller gains per actuator
```

Why this shape:
- **Tree + constraints, not just a tree**: every serial/parallel arm, robot
  hand, hexapod leg, vehicle suspension, or cable mechanism reduces to
  (tree joints) + (loop/contact constraints). MJCF's own model = bodies +
  joints + constraints, so the converter is near-mechanical.
- The current `SerialManipulatorConfig` becomes a degenerate case
  (one open chain of revolute joints).

Validation (before execution, spec §22/§55): joint-axis consistency, parent
frame existence, inertia positivity, graph is a tree (no cycles outside the
constraint list), constraint bodies exist, actuators reference joints.
Errors carry "what/why/which/alternatives" diagnostics.

---

## 4. Solver design decision (recorded here so implementation is mechanical)

- **Open kinematic trees: minimal (generalized) coordinates.** Extend the
  current N-link machinery (Jacobian/Christoffel or composite-rigid-body +
  recursive Newton–Euler) to a tree with per-joint DOF types. Exact, no
  constraint drift, no constraint machinery for open chains (the common
  robotics case). Keep the existing per-joint torque limits/stops/PD.
- **Closed loops + contacts + prescribed constraints: maximal coordinates
  with a soft-constraint iterative solver.** Bodies float in 6(7)-DOF
  maximal coordinates; joints AND loop closures become bilateral
  constraint equations; contacts are unilateral (complementarity);
  solve by projection (Gauss–Seidel/PGS on the constraint-Jacobian,
  graph coloring for parallel iteration — spec §30's "constraint graphs,
  graph coloring, iterative solvers"). This is the MuJoCo-proven design and
  the only way to get arbitrary constraint graphs without a symbolic
  reduction engine.
- **This hybrid is deliberate**: reduced coordinates are superior for open
  chains (robotics realtime), maximal+constraints are superior for loops and
  contacts — both share the model description, the state, the actuator
  layer and the validation.

Integration: `ManipulatorState` → `MultibodyState` (per-body pose +
per-joint coordinates); controls stay PD/actuator-classes; the existing
`presets::robotics::manipulator_trajectory` generalizes to the description
with trajectory + control configs.

---

## 5. Roadmap (not scheduled; implementation order when taken)

| Phase | Deliverable | Verification anchors |
|---|---|---|
| M1 | Tree generalization: arbitrary joint types in minimal coordinates (branching trees, per-joint DOFs, limits/damping/friction) | pendulum/generic double-pendulum energy, gyroscopic conservation, joint-type sanity (slide has zero rotational coupling), branching tree energy |
| M2 | Maximal-coordinate dynamics + bilateral constraint solver (loop closure, distance) | four-bar closure (loop remains closed, constraint drift bounded), planar loop energy, constraint-Jacobian symmetry checks |
| M3 | Unilateral contacts: spatial-hash broad phase + simplex narrow phase + PGS/LCP | collision restitution energy, stacking stability, friction cone behavior, contact determinism |
| M4 | Actuator/control layer beyond PD (position/velocity/torque servos, paths) + end-effector wrenches | servo step-response, torque-limit saturation, trajectory tracking on M1–M3 rigs |
| M5 | Model-description reader + validation pipeline; converter examples against MJCF-style models | schema round-trip, invalid-model diagnostics, capstone: an imported robot arm + payload + loop (e.g., four-bar gripper) on one description |

Symmetry with the rest of the engine: the description validates through the
`coupling::rules` machinery; the state feeds the `core::State`/diagnostics
layers; output/visualization consumes the same state (rendering is external).

---

## 6. What the external repo should prepare (so the boundary is real)

1. Author geometry in MuJoCo (MJCF) and/or any future format.
2. Mesh loaders + converters that emit `MultibodyDescription` data
   (version 0.1 schema; the field names above are the contract).
3. For non-multibody use cases, the same converter fills the existing
   thermal/fluid/structural configs (regions, BCs, parameters) from the
   geometry — "import CAD → boom the solver runs" per domain class.
4. Rendering/visualization consumes engine state (GPU-resident when the
   CUDA backend lands, spec Phase 11–12) — rendering never lives here.

The engine's obligation, tracked in the conformance doc: keep every public
config data-driven and format-neutral so the converter can always target
plain structs; never grow format parsers.
