# Multibody Description — Schema v0.1 (the engine-owned contract)

**Status: STABLE (agreed 2026-09).** The authoritative, format-neutral model
description the external geometry repo emits and the extropian-physics
multibody engine consumes.  Versioned; breaking changes require a version
bump.  Companion to `docs/multibody_architecture.md`; the C++ struct
layouts below are the contract (the engine may additionally ship a reader
for this OWN schema — never for third-party formats).

---

## 1. Conventions (non-negotiable)

| Item | Convention |
|---|---|
| Units | SI: meters, kilograms, seconds, radians, N, N·m, Pa, kg·m², W |
| World frame | right-handed Cartesian, **Z up**. Default gravity `{0, 0, −9.81}` m/s² |
| Angles | radians everywhere (joint ranges, initial q, axes free of units) |
| Quaternions | **active**, unit-norm, **order (w, x, y, z)**, w = cos(θ/2) — identical to `physics::rigid_body::quaternion` |
| Body frames | inertial data (COM, inertia) expressed in the BODY frame: origin at the joint anchor, axes initially aligned with the parent frame, rotated by the joint DOFs |
| Parent frame | every joint `anchor` and `axis` are expressed in the PARENT body frame; a `parent == ""` body is attached to the WORLD frame |
| Tree rule | bodies form a tree via `parent`; arrays need no order (the engine topologically sorts). Cycles are invalid |
| State flattening | q/dq vectors concatenate joint blocks in depth-first order, parents before children (algorithm in §6) |
| Text encoding | UTF-8; identifiers match `[A-Za-z_][A-Za-z0-9_]*`, unique within one scope (body names, joint names implied by bodies) |
| No exceptions | validation is a function of the description via `ModelStatus` |

---

## 2. Top-level record

```cpp
struct MultibodyDescription {
    uint32_t version = 1;                  // schema v0.1
    EnvironmentSpec environment;           // gravity, dt, integrator hints
    std::vector<BodySpec> bodies;          // the kinematic tree
    std::vector<ConstraintSpec> constraints;  // non-tree constraints (may be empty)
    std::vector<ActuatorSpec> actuators;   // may be empty
    std::vector<LoadSpec> loads;           // may be empty
    InitialStateSpec initial_state;        // per-joint q0/dq0
    std::string name;                      // model name (diagnostics)
};

struct EnvironmentSpec {
    std::array<double,3> gravity = {0, 0, -9.81};  // m/s², world frame
    double dt_hint = 1e-3;                 // s; the engine may adapt/validate
    // integrator hint: "rk4" (default) | "symplectic" | "implicit" — hints
    // only; the engine picks a validated method and may override.
    std::string integrator = "rk4";
};
```

## 3. Bodies and joints

```cpp
enum class JointType : uint8_t {
    Revolute,   // 1R about `axis`
    Prismatic,  // 1T along `axis`
    Screw,      // 1R about `axis`; translation = ang·pitch along `axis`
    Cylindrical,// 2 DOF: rotation + translation along `axis`
    Universal,  // 2R: about `axis`, then about `axis2` (perpendicular, rotates with the first hinge)
    Planar,     // 3 DOF: 2T in the plane ⊥ `axis` + 1R about `axis`
    Ball,       // 3R, quaternion (unit constrained)
    Free,       // 6 DOF rigidly free (quaternion, unit constrained)
    Weld,       // 0 DOF: rigidly attaches the child to the parent at `anchor`
};

struct JointSpec {
    JointType type = JointType::Revolute;
    std::array<double,3> anchor = {0,0,0};  // joint point in the PARENT frame (m)
    std::array<double,3> axis = {0,0,1};    // unit direction in the PARENT frame
    std::array<double,3> axis2 = {0,1,0};   // Universal only; ⊥ axis (validated)
    double pitch = 0.0;                     // Screw only: translation per radian (m/rad)

    // soft limits (1-DOF joints only in v0.1): range = [min, max] rad or m;
    // min >= max disables the stop.  stiffness/damping = the soft-stop
    // penalty (rotational: N·m/rad, N·m·s/rad; translational: N/m, N·s/m).
    std::array<double,2> range = {0, 0};
    double stiffness = 500.0;
    double damping = 20.0;
};

struct BodySpec {
    std::string name;
    std::string parent;                     // "" = world
    JointSpec joint;                        // default: Revolute at the anchor
    double mass = 1.0;                      // kg, > 0
    std::array<double,3> center_of_mass = {0,0,0};  // BODY frame (m)
    std::array<double,9> inertia = {0,0,0, 0,0,0, 0,0,0};
        // FULL symmetric 3×3 inertia tensor ABOUT THE COM, BODY frame, kg·m²,
        // row-major [ixx ixy ixz  ixy iyy iyz  ixz iyz izz].
        // All zeros = point mass at the COM (valid; the engine uses the
        // Jacobian-only dynamics).  Validated: symmetric, positive definite
        // when nonzero.
    std::string visual_mesh = "";           // METADATA ONLY: name the external
                                            // repo resolves to its own mesh
                                            // asset; the engine never loads it.
    std::vector<double> dof? = 0;           // (reserved; not in v0.1)
};
```

Root rule: `parent == ""` with `Weld` = bolted to the world; with `Free` =
floating base.  A rooted `Revolute/…` joint on `parent == ""` means the
joint anchors to the world frame (a hinge fixed in space).

## 4. Constraints (non-tree)

```cpp
enum class ConstraintType : uint8_t { Loop, Distance };

struct ConstraintSpec {
    std::string name;
    ConstraintType type = ConstraintType::Loop;
    std::string body_a;                     // must exist, != body_b
    std::string body_b;
    std::array<double,3> anchor_a = {0,0,0};  // local point in body_a
    std::array<double,3> anchor_b = {0,0,0};  // local point in body_b
    bool align_orientations = false;        // Loop only: also equate the body
                                            // orientations (their frames)
    std::array<double,2> distance_range = {0, 0};  // Distance: [min, max] m;
                                            // min == max = exact distance
    double stiffness = 1e6;   // soft-constraint penalty scale (N/m for
    double damping = 100.0;   // distance/loop position; N·m/rad for align)
};
```

`Loop` semantics: the world position of `anchor_a` equals the world position
of `anchor_b` (optionally also the orientations).  `Distance`: the world
distance between the anchors stays within `[min, max]`.  All constraints are
SOFT in v0.1 (penalty/PGS-style; the engine exposes the achieved residual as
a diagnostic).  These are the only non-tree constraint types in v0.1;
contact/collision constraints are a later version.

## 5. Actuators, loads, initial state

```cpp
enum class ActuatorType : uint8_t { Torque, Position, Velocity };

struct ActuatorSpec {
    std::string name;
    ActuatorType type = ActuatorType::Torque;
    std::string joint;                      // joint to actuate (the body name)
    uint8_t dof_index = 0;                  // 0 for 1-DOF joints; index into the joint's DOF list
    double gain = 1.0;                      // scaling on the computed command
    double limit = 0.0;                     // command clamp: torque N·m / force N;
                                            // 0 = unlimited
    // Position/Velocity servos (rotational: kp N·m/rad, kd N·m·s/rad;
    // translational: kp N/m, kd N·s/m, limit N):
    double kp = 100.0;
    double kd = 10.0;
    double reference = 0.0;                 // constant reference in v0.1
    // (trajectory references arrive with the control layer, M4)
};

struct LoadSpec {
    std::string body;
    bool in_world_frame = true;             // false = body frame
    std::array<double,3> force = {0,0,0};   // N, applied at `point`
    std::array<double,3> torque = {0,0,0};  // N·m about the body COM
    std::array<double,3> point = {0,0,0};   // application point in the chosen frame
};

struct InitialStateSpec {
    // per-joint q0/dq0, keyed by the body name; missing joints = zeros.
    // The block layout PER JOINT TYPE is fixed (§6) — the engine stitches
    // the flattened state in tree order.
    std::vector<JointInitialSpec> joints;

    struct JointInitialSpec {
        std::string body;
        std::vector<double> q0;    // e.g. Revolute: [θ]; Ball: [w,x,y,z];
                                   // Free: [x,y,z, w,x,y,z]
        std::vector<double> dq0;   // e.g. Revolute: [ω]; Free: [vx,vy,vz,ωx,ωy,ωz]
    };
};
```

## 6. Joint DOF table (the state-layout law)

| `JointType` | DOFs | q block (entries) | dq block (entries) | Notes |
|---|---|---|---|---|
| Revolute | 1 | `[θ]` (1) | `[ω]` (1) | axis positive = right-hand rotation |
| Prismatic | 1 | `[s]` (1) | `[ṡ]` (1) | positive along `axis` |
| Screw | 1 | `[θ]` (1) | `[ω]` (1) | translation along the axis = θ·pitch |
| Cylindrical | 2 | `[θ, s]` (2) | `[ω, v]` (2) | rotation then translation along `axis` |
| Universal | 2 | `[θ1, θ2]` (2) | `[ω1, ω2]` (2) | θ1 about `axis`; θ2 about `axis2` |
| Planar | 3 | `[sx, sy, θ]` (3) | `[vx, vy, ω]` (3) | 2T in the plane ⊥ `axis` (engine picks the in-plane basis from `axis`), then rotation about `axis` |
| Ball | 3 (4 stored) | `[w,x,y,z]` (4, unit) | `[ωx,ωy,ωz]` (3) | the engine projects to unit norm |
| Free | 6 (7 stored) | `[x,y,z, w,x,y,z]` (7, unit) | `[vx,vy,vz, ωx,ωy,ωz]` (6) | pose in the parent (= world for roots) frame |
| Weld | 0 | `[]` (0) | `[]` (0) | rigid attachment |

Flattening: depth-first, parents before children; append each joint's q
block (and dq block) in the table order.  Identifiers from this layout give
the per-joint slices for actuators, diagnostics, and output.

## 7. Validation rules (numbered; implemented as `coupling::rules`)

Failure wording style: *what is invalid, why, which entity, alternatives.*

| # | Rule | Failure example |
|---|---|---|
| V1 | version supported (== 1) | "multibody: description version 2 unsupported (schema v0.1)" |
| V2 | body names unique, non-empty | "multibody: duplicate body 'link2'" |
| V3 | `parent` exists (or "") | "multibody: body 'hand' parent 'wrist2' not found" |
| V4 | the `parent` graph is a tree (no cycles) | "multibody: parent cycle detected through 'link3'" |
| V5 | `mass > 0` | "multibody: body 'base' mass must be > 0" |
| V6 | `axis`/`axis2` are non-zero unit vectors (within 1e-6) | "multibody: body 'elbow' joint axis not a unit vector" |
| V7 | Universal `axis2` ⊥ `axis` (within 1e-3 rad) | "multibody: body 'wrist' universal axes not perpendicular" |
| V8 | `inertia` symmetric (|aᵢⱼ − aⱼᵢ| ≤ 1e-6·max) and positive definite when nonzero | "multibody: body 'gripper' inertia not positive definite" |
| V9 | `range` with `min > max` disables the stop; `min == max` is allowed | (no error; documented) |
| V10 | constraints: bodies exist, `body_a != body_b`, anchors finite | "multibody: constraint 'loop1' references unknown body 'thing'" |
| V11 | Distance `min ≤ max` | "multibody: constraint 'strut' distance min > max" |
| V12 | actuators: `joint` exists, `dof_index` within the joint's DOF count | "multibody: actuator 'motor2' dof_index 1 out of range for a Revolute joint" |
| V13 | loads: body exists | "multibody: load 'thrust' references unknown body" |
| V14 | initial state: body exists, q0/dq0 lengths match the joint table | "multibody: initial state for 'ball_joint' needs 4 q entries (w,x,y,z), got 3" |
| V15 | controller/actuator units match the DOF kind (rotational gains on rotational DOFs) | "multibody: actuator 'lift' uses rotational gains on Prismatic joint 'lift_joint'" |
| V16 | at least one body; the tree is connected through parents | "multibody: description has no bodies" |

## 8. JSON mapping (the external repo's serialization)

The C++ structs are authoritative; the JSON form uses the same names
(enum values serialized as their PascalCase names: `"Revolute"`,
`"Loop"`, `"Torque"`, …; quaternions as `[w,x,y,z]`; 3×3 inertia as a
9-element array).  Example:

```json
{
  "version": 1,
  "name": "two_link_arm",
  "environment": { "gravity": [0, 0, -9.81], "dt_hint": 0.001 },
  "bodies": [
    { "name": "base",  "parent": "",  "joint": { "type": "Weld" },
      "mass": 1.0, "center_of_mass": [0,0,0], "inertia": [0,0,0,0,0,0,0,0,0] },
    { "name": "link1", "parent": "base",
      "joint": { "type": "Revolute", "anchor": [0,0,0.2], "axis": [0,0,1],
                 "range": [-2.0, 2.0], "stiffness": 500, "damping": 20 },
      "mass": 1.0, "center_of_mass": [0,0,0.25], "inertia": [0,0,0,0,0,0,0,0,0] },
    { "name": "link2", "parent": "link1",
      "joint": { "type": "Revolute", "anchor": [0,0,0.5], "axis": [0,0,1] },
      "mass": 0.8, "center_of_mass": [0,0,0.2], "inertia": [0,0,0,0,0,0,0,0,0] }
  ],
  "constraints": [],
  "actuators": [
    { "name": "shoulder", "type": "Torque", "joint": "link1", "limit": 5.0 },
    { "name": "elbow",     "type": "Torque", "joint": "link2", "limit": 2.0 }
  ],
  "loads": [],
  "initial_state": { "joints": [
    { "body": "link1", "q0": [0.0], "dq0": [0.0] },
    { "body": "link2", "q0": [0.0], "dq0": [0.0] }
  ]}
}
```

## 9. Correspondence (so existing code and MJCF map near-mechanically)

**Current `SerialManipulatorConfig` → description:** `links[i]` → a body
named `link_i` with `mass = links[i].mass`, a Revolute joint
(`axis = (0,0,1)`; planar arms lie in the XY world plane, Z-up), a Weld
root body; `q_min/q_max` → `range`; `torque_max` → actuator `limit`;
`PidGains` → Position actuators.  The engine's existing tests remain valid
as the n-link regression suite.

**MJCF → description (converter guide):** `body` → `BodySpec` (`pos` →
`anchor` in the parent frame; `inertial` → `mass`, `pos` →
`center_of_mass`, `fullinertia`/`diaginertia` → `inertia`), `hinge` →
Revolute, `slide` → Prismatic, `ball` → Ball, `free` → Free, `axis` →
`axis`, `range` → `range`, `damping/stiffness` → `damping/stiffness`;
`motor` → Torque actuator, `position` → Position actuator, `velocity` →
Velocity actuator.  MJCF compound joints (universal/planar via paired
hinges) map to the native `Universal`/`Planar` types.  World frame and
quaternion order match exactly, so the converter is field-for-field.

## 10. Versioning & evolution policy

- **v0.1 is a stable contract**: the external converter targets it; the
  engine validates it (rules V1–V16).
- Additive changes (new joint types, new constraint types, new actuator
  types) bump the MINOR field (`version` stays 1, the doc gains a dated
  addendum) and must not break v0.1 consumers.
- Breaking changes (frame conventions, quaternion order, block layouts)
  bump `version` to 2 with a migration note; the engine validates and
  rejects unmatched versions with a pointer to the migration.
- The engine never parses third-party formats; the schema IS the
  interchange.  The external repo's MJCF/URDF/CAD pipeline ends at
  `MultibodyDescription` data (C++ structs or the JSON form).
