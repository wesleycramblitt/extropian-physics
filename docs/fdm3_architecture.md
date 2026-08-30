# 3D FDM Fluid Solver (fluid::fdm3) — Architecture

Status: **implemented, validated** (2026-08-29)
Namespace: `exd::physics::fluid::fdm3`
Sources: `src/fluid/fdm3/*.cpp`, public headers `include/exd/physics/fluid/fdm3/`

## 0. Scope

3D incompressible Navier–Stokes on a regular structured grid. Collocated
cell-centered variables (u, v, w, p at cell centers) with one ghost layer —
mirrors the legacy 2D FDM (`fluid::fdm`, which is kept untouched) rather
than a true MAC layout (documented deviation; the old docs' "staggered"
claim was wrong). Engineering-grade fidelity: pressure-projection (SIMPLE),
SOR Poisson, explicit time integration. The solver is the fluid partner of
the coupled turbine-in-grid (`docs/turbine_coupling_architecture.md`).

## 1. Public API

```
fdm3_config.hpp   BoundaryFace {XMin,XMax,YMin,YMax,ZMin,ZMax}
                  FDM3BoundaryCondition {face, type, u/v/w/p values, paired_face}
                  FDM3Config {nx,ny,nz, lx,ly,lz, rho, mu, dt, max_steps,
                              time_integration, advection_scheme, cfl_target,
                              adaptive_dt, pressure/SIMPLE knobs, BCs,
                              initial field, field_stamp_interval; validate()}
                  (reuses fluid::fdm enums: TimeIntegration, AdvectionScheme,
                   FDMBoundaryType)
fdm3_result.hpp   FDM3FieldData {u,v,w,p, x,y,z, index(i,j,k)}
                  FDM3StepResult {time, step, residuals, max_velocity, cfl, divergence}
                  FDM3Result {valid, error, warnings, field, history, steps_taken,
                              final_time, converged}
fdm3_solver.hpp   FDM3Solver — persistent-state stepper:
                    initialize(config, status) / step(dt, status) / time() /
                    step_count() / field() / last_step() / set_body_force(fx,fy,fz)
                  solve_fdm3(config)            — one-shot pure run
                  run_fdm3_simulation(config, writer?, scheduler?, status?)
                                                — driver: stamps fields at a cadence
```

Design rule (repo doctrine): `FDM3Solver::step` advances the SAME grid —
the 2D `step_fdm` defect (re-initializes every call) is not repeated. The
pure `solve_fdm3`/`run_fdm3_simulation` shapes keep optimizer batches free
of I/O (writers are driver-owned).

## 2. Method

- **Layout**: interior `nx×ny×nz` cells + one ghost layer per axis,
  `idx(i,j,k) = i + sx·(j + sy·k)`, `sx=nx+2`, `sy=ny+2`.
- **Predictor**: momentum equations integrated explicitly
  (ForwardEuler / Heun / RK4 / CrankNicolson-as-explicit-trapezoid —
  honest label in code; true implicit CN is a linear-solve follow-up).
  Advection Central/Upwind/Hybrid per axis. **Body force** (acceleration
  m/s² per interior cell, set via `set_body_force`) enters the predictor
  ONLY — the projection keeps the force divergence-consistent.
- **Projection (SIMPLE)**: `∇²p' = div(u*)/dt` (standard sign — the 2D
  solver's `+rhs` sign never actually decomposed divergence; corrected and
  documented here), 7-point SOR point relaxation with configurable ω,
  ghost Dirichlet from BCs, periodic ghosts refreshed during relaxation,
  then `u = u_old + α_u·(u* − u_old − dt·∇p')`, `p = p_old + α_p·p'`.
- **BCs**: all 6 faces via a generic per-face loop — Inlet (specified
  velocity), Outlet (zero-gradient), Wall (reflected ghosts),
  Symmetry (normal-reflected), Periodic (paired-opposite faces, verified
  naive pairing), FixedPressure (Dirichlet p ghost).

## 3. Validation (tests/unit/fdm3/, all green)

| Case | Result |
|---|---|
| Uniform-flow preservation (12³) | exact (u ≡ const, div ≡ 0) |
| 3D plane Poiseuille | centerline within 0.13% of analytic |
| Taylor–Green vortex decay | E(t)/E₀ within 4.7% of e^(−2νk²t) |
| Fixed actuator disk (momentum theory) | far-wake deficit within 22% of (1−2a); sign pinned |
| SOR Poisson vs analytic ∇²p = 6 | L∞ < 1e-4 at 12³ |

`-O2` is applied to these sources only (repo default build is unoptimized;
otherwise the 16–32³ kernels dominate test time).

## 4. Limitations (documented)

- Collocated scheme: no Rhie–Chow; checkerboard risk tamed by the
  projection + full velocity under-relaxation in coupled runs.
- Structured grids only; unstructured/FVM/FEM is Phase J.
- No turbulence model, single-phase incompressible; no free surface.
- Periodic pairing limited to natural opposite faces (v1).
```
