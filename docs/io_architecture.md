# Output Channels (io) — Architecture

Status: **implemented, validated** (2026-08-29)
Namespace: `exd::physics::io`
Sources: `src/io/*.cpp`, public headers `include/exd/physics/io/`
Format contract: `docs/output_channels.md` (authoritative for consumers)

## 0. Scope

Real-time-capable output for solvers and drivers: 3D field stamps (binary,
fast to write, trivial to parse) and machine-state time series (CSV), with
a cadence policy that supports deterministic offline replay AND
wall-clock-throttled "real-time if specified" output.

## 1. Public API

```
field_writer.hpp  FieldGeometry {origin, spacing, dims} — dims are SAMPLE
                  COUNTS, origin is the FIRST SAMPLE (cell centers for FDM)
                  IFieldWriter (leaf seam):
                    begin_stamp(t, step) / write_scalar_field(name, geo, f32)
                    / write_vector_field(name, geo, f32 interleaved xyz)
                    / end_stamp / name()
                  FldWriterConfig {directory, overwrite, sanitize_nonfinite}
                  make_fld_writer(config, status) → exd-fld v1 container
                  make_null_writer() → no-op sink (optimizer default)
series_writer.hpp CsvSeriesWriter {path, columns, flush_each_row, status}
                  — parent directories auto-created; one line per row
output_policy.hpp OutputPolicy {every_n_steps, wall_clock_interval_s} — PURE
                  OutputScheduler {policy + mutable last-emit state,
                  set_now(), should_emit(step)} — clock injected by the
                  caller (deterministic tests)
```

Design rule (repo doctrine): writers are LEAF infrastructure — solvers
(`solve_*`, `step_*`) never write; only drivers (`simulate_engine`,
`run_coupled_turbine`, `run_fdm3_simulation`) own stamps. `solve_*` stays
pure and batchable; the optimizer runs with null sinks (zero I/O).

## 2. exd-fld v1 binary container

One file per stamp, little-endian, fixed 32-byte header + per-field table
(164 B/field: name, kind, element type, count, payload offset/bytes,
per-field origin/spacing/dims) + raw float32 payloads. Non-finite values
sanitized to 0.0f (configurable). Timeline: `timeline.txt`, one line per
stamp (`time step filename`) — the animation clock. The consumer
implementation is ~150 lines; full byte-level spec in
`docs/output_channels.md`.

## 3. Cadence policy

`OutputScheduler` owns the last-emit bookkeeping; the caller advances the
clock (steady clock for wall-clock real-time, or sim time for
deterministic tests). Wall-clock throttle (e.g., 0.25 s ≈ 4 fps) keeps the
I/O load bounded regardless of compute step size; `every_n_steps` gives
byte-identical offline replay.

## 4. Validation (tests/unit/io/)

| Case | Result |
|---|---|
| Round-trip scalar+vector via a spec-faithful mini-reader | exact |
| NaN/Inf sanitization | 0.0f |
| Duplicate field names reject the stamp | ✓ |
| Geometry/size mismatch rejected; empty stamp allowed | ✓ |
| Bad directory → clean factory failure | ✓ |
| CSV header + `%.9g` rows; wrong column count rejected | ✓ |
| CSV parent-directory creation | ✓ |
| Step trigger deterministic; wall-clock with injected clock | ✓ |
| OR semantics; disabled policy never emits | ✓ |
```
