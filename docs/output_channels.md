# Output Channels — Real-Time Field & Machine-State Output

Authoritative spec for what extropian-physics emits and what a consumer
(animation/visualization repo, dashboard, optimizer) must parse.

Status: **v1 implemented** (2026-08-29). Wave 1 of the phased program
(`docs/modular_solver_architecture.md`).

---

## 1. Contract overview

| Channel | Purpose | Format |
|---|---|---|
| Field stamps | 3D scalar/vector fields (velocity, pressure, temperature…) per time step | binary `exd-fld v1` (`.fld` files) + `timeline.txt` manifest |
| Machine-state series | Scalars over time: ω, angle, piston x, torque, power, p_cyl… | CSV (one file per series) |

Drivers (`simulate_turbine`, `simulate_engine`, `run_coupled_turbine`) own
writer stamps. Solvers (`solve_*`, `step_*`) are pure and never write —
optimizer batchability is preserved. Writers default to a `null` sink;
output is on only when the caller supplies a writer/policy.

The typical real-time setup: `OutputPolicy{wall_clock_interval_s: 0.25}`
(≈4 fps) throttles stamps to a stable animation cadence regardless of the
compute step size; `every_n_steps` gives deterministic offline replay.

---

## 2. `exd-fld` v1 binary container (`.fld`)

One file per stamp. **All integers/floats little-endian.** Field payloads
are IEEE-754 float32 only.

### Layout

```
[0]     8   magic "EXDFLD01"
[8]     4   u32 version = 1
[12]    4   u32 field_count
[16]    8   f64 time (simulation time, s)
[24]    8   reserved (zero)
[32]        field table, 164 bytes per field, write order:
            [0]    64  char name[64]  (NUL-padded)
            [64]    1  u8 kind        (0 = scalar, 1 = vector)
            [65]    1  u8 element_type (0 = float32)
            [66]    2  reserved (zero)
            [68]    4  u32 count      (elements: scalar = nx·ny·nz,
                                       vector = 3·nx·ny·nz)
            [72]    4  u32 payload_offset (byte offset from file start)
            [76]    4  u32 payload_bytes  (= count·4)
            [80]   24  f64 origin[3]  (location of FIRST SAMPLE)
            [104]  24  f64 spacing[3] (per-axis sample spacing, > 0)
            [128]  12  u32 dims[3]    (sample point counts per axis)
            [140]  24  reserved (zero)
after table: payloads, raw float32, in table order
```

### Conventions (consumers MUST follow)

- **`dims` are sample point counts**, `origin` is the location of the first
  sample. FDM fields are **cell-centered**: `origin = (dx/2, dy/2, dz/2)`
  relative to the domain corner, `dims = (nx, ny, nz)` grid cells.
  Do NOT assume node-centered structured data.
- **Vector fields are interleaved** xyz per sample:
  `v[3*i+0]=vx, v[3*i+1]=vy, v[3*i+2]=vz`, index `i = ix + nx·(iy + ny·iz)`.
- Non-finite payload values are sanitized to `0.0f` by the writer
  (configurable off). A zero value is therefore "no data", not physics.
- Field names are ≤ 63 chars, unique within a stamp (duplicates reject the
  stamp).

### `timeline.txt` manifest

Appended one line per stamp in the output directory:

```
<time> <step> <filename>
```

e.g. `12.5 7 step_00000007.fld`. `time` is `%g`-formatted (9 sig. digits).
The animation repo reads this in order and loads each file; `time` is the
authoritative animation clock. `overwrite=true` (default) truncates the
manifest at writer construction.

---

## 3. CSV machine-state series

One CSV file per channel family, e.g. `engine_state.csv`:

```
time,omega,angle_rad,piston_x,piston_v,p_cyl,torque,power
0,0,0,0.05,0,101325,0,0
0.0004,18.3,0.00732,0.0499,-11.2,102100,8.4,154
```

- `time` first column, `%.9g` formatting, one row per stamp, header row
  always present. Column names are simple identifiers (no commas/quotes).
- `flush_each_row` option makes rows visible immediately for live
  consumers; otherwise the OS buffer batches (better throughput, slightly
  stale tail — flush on close).

## 4. Real-time guidance

- **Cadence**: 0.25–1.0 s wall-clock between stamps reads smoothly in
  renderers at 64³ grids (~1–4 MB/stamp float32). Never write every step
  at high resolution — I/O dominates the compute loop.
- **Throttle is clock-injected**: drivers pass an abstracted `now()`
  (steady clock), so tests are deterministic and the animation repo's
  time base can be simulation time.
- **Memory**: writers buffer one stamp; storage can be moved to a worker
  thread later (snapshot semantics already in place).

## 5. Interaction with optimization

`OutputPolicy`/`OutputScheduler` state is driver-owned; `solve_*` entry
points take no writer parameters. Optimizer batches call `solve_*` with
output disabled (null sink, no policy) — zero I/O overhead in the loop.

## 6. Future writers

`IFieldWriter` (`include/exd/physics/io/field_writer.hpp`) is the seam:
VTK/VTI, HDF5, or a shared-memory streaming transport implement it without
touching solvers. The animation repo only needs the exd-fld parser above.
