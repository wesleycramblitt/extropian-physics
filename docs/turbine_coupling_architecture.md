# Coupled Turbine-in-Grid (turbine) — Architecture

Status: **implemented, validated** (2026-08-29)
Namespaces: `exd::physics::turbine` (driver, builder),
`exd::physics::fluid::forces` (BladeElement evaluator)
Sources: `src/turbine/coupled_turbine.cpp`, `src/fluid/forces/blade_element.cpp`

## 0. Scope

Actuator-disk-style rotor inside the 3D FDM (`docs/fdm3_architecture.md`):
the fluid field resolves the induction, the rotor responds with torque, and
the exchange loop feeds the rotor's per-blade loads back as a smeared
volume force. This is `docs/real_solver_for_turbine.md`'s
CFD → forces → rotor → CFD loop, made real. Also supplies the parametric
`make_turbine_definition` builder so real runs need no hand-built geometry.

## 1. Public API

```
fluid::forces (additive):
  ForceEvaluatorType::BladeElement
  BladeElementConfig {}
  IForceEvaluator::compute(blade, flow, omega, axis, per_element, status)

turbine:
  make_turbine_definition(TurbineBuilderConfig, status) → exd-geometry
      TurbineDefinition (hub/tip radius, chord, linear twist, rpm, blades,
      duct length/radius, airfoil key; validated)

  default_grid_config(v_inf, n_per_axis, radius_mult, length_mult)
      → FDM3Config with the −Z inflow convention
        (Inlet at +z w=−v_inf, Outlet at −z, lateral Symmetry,
         uniform init w = −v_inf, CFL dt)

  CoupledTurbineConfig {grid, turbine, element_count, rotor_inertia,
      generator (CurveMomentConfig), rotor_origin, dt,
      fluid_steps_per_exchange, force_relaxation β, ramp_time_s τ_r,
      smear_cells ε, max_steps, history, field_writer/output_scheduler/
      field_stamp_interval, csv_path}
  CoupledTurbineStep {t, omega, angle_rad, torque, axial_force, power, exchange}
  CoupledTurbineResult {valid/error/warnings, fluid (fdm3 summary),
      history, final_omega/cp/tsr, exchanges, aero_work, load_work,
      rotor_ke_change, fluid_ke_change}
  run_coupled_turbine(config, status)
```

## 2. The coupling loop (per exchange window)

```
1. Sample   the CFD field at 8 azimuths × radius r_e in the rotor plane
            (world = grid coordinates, axis through rotor_origin, +Z).
2. Evaluate per-station forces with the BladeElement model:
            u_a = −(v_mean·e_z)   (azimuth-MEAN of the element's samples)
            V_t = ω·r_e,  W = √(u_a² + V_t²),  α = φ − twist,
            cl/cd from the polar database (Re-based selection),
            dL/dD per element, distributed per blade in the lab frame.
            NO induction solve — the field IS the induction
            (documented anti-double-count rule).
3. Integrate → per-station ring thrust T_e and torque Q_e.
4. Build the body force (fluid receives the NEGATED blade load):
            f_z   +=  +T_e·η(z̃)/(ρ·2π r_e dr_e)      (decelerates −Z inflow)
            f_θ   +=  −Q_e·η(z̃)/(ρ·2π r_e² dr_e)     (reaction swirl)
            ξ = e_t·(−sinφ, cosφ, 0);  ε = smear_cells·min(dx,dy,dz);
            η(z̃) = exp(−z̃²/2ε²)/(ε√(2π)), z̃ = z − z_rotor
5. Under-relax + ramp: f ← min(1, t/τ_r)·(β·f_new + (1−β)·f_old)
   (paper behind β: classic actuator-disk damping; τ_r ≥ 10·window
   prevents forcing-discontinuity excitation — validated).
6. set_body_force(solver, …)  (predictor-only, projection-consistent).
7. Advance the rotor over the window: net = ΣQ_e − T_load(ω);
   staggered explicit coupling, Heun rotor dynamics.
```

Stability guards enforced in validation: `dt < ε/(3·V_inf)`, β ∈ (0,1],
ramp ≥ 10·window, Inlet BC present with −Z convention (warning on
mismatch), rotor disc inside the box; force sample out-of-bounds →
stall-guard zero forces + warnings.

## 3. Why BladeElement (vs MomentumBalance in coupled mode)

The standalone `MomentumBalance` evaluator runs its own Buhl/Glauert
induction solve against one uniform inflow sample. Feeding it local
(already-decelerated) velocities would double-count the induction and
under-produce (driving the a-clamp). The coupled variant therefore uses
the local azimuth-mean velocity DIRECTLY — actuator-line practice.

## 4. Validation (tests/unit/turbine/coupled_turbine_test.cpp)

| Case | Result |
|---|---|
| BladeElement unit forces (finite, correct shape) | ✓ |
| Spin-up soak (20³, R_tip 0.4, V 1 m/s): settles, no NaN | ω* = 6.17 rad/s, spread 8% |
| Wake deceleration sign (positive disk slows −Z flow) | ~20% deficit |
| cp vs standalone MomentumBalance at same ω | ratio 0.63 ∈ [0.4, 2.0] |
| Energy balance: aero_work ≥ ΔKE_rotor + load_work + ΔKE_fluid | closed |
| Determinism (same config → identical trajectory) | 1e-12 |
| Validation guards (no Inlet, dt too large, ramp too short) | errors |
| Field writer stamps + timeline at cadence | ✓ |
| Rotor CSV: one row per fluid step | ✓ |
| Hydro: seawater ρ/μ with ρ-scaled load+inertia settles, Betz-bounded | ✓ |

## 5. Limitations (documented)

- Structured grid only; grid resolution/cost trade-off noted in
  `docs/real_run_guide.md` (20³–40³ seconds per run, 48³+ minutes).
- No blade structural response (per-element 3D forces ARE exported for a
  future structural solver); no tower/yaw/nacelle; no free-surface.
- Frozen-force windows: fields update every `fluid_steps_per_exchange`
  substeps (piecewise-constant field view — documented cosmetic).
- Formal `CouplingManager` exchange rewrite (Phase H) is deferred; the
  driver IS the working demo of the exchange pattern.
```
