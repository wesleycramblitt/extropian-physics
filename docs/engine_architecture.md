# Engine Simulation (engine) — Architecture

Status: **implemented, validated** (2026-08-29)
Namespace: `exd::physics::engine`
Sources: `src/engine/*.cpp`, public headers `include/exd/physics/engine/`,
plus the steam EOS in `exd::physics::thermo` (`include/exd/physics/thermo/steam.hpp`).

## 0. Scope

Single-cylinder slider-crank engine — combustion (Otto, 4-stroke, Wiebe
heat release) and steam (single-acting, Rankine-lite). 0D/1D fidelity
(industry-standard for cycle studies); the cylinder is a thermodynamic
state machine, not a mesh. Couples to loads (friction + generator curves),
a PI governor, and per-step time-series output (CSV).

## 1. Public API

```
engine_config.hpp   EngineCycleType {Otto, Steam}
                    EngineGeometryConfig {crank_radius, rod_length, bore,
                        clearance_volume, piston_mass, flywheel_inertia}
                    EngineThermoConfig {cycle, r_gas, p_intake, T_intake,
                        p_exhaust, gamma_compression, gamma_expansion,
                        q_in_cycle, Wiebe params (ignition, duration, a, m),
                        steam: p_boiler, p_condenser, cutoff, gamma=1.13,
                        steam_quality_cutoff}
                    EngineLoadConfig {friction_constant, friction_viscous,
                        generator curve}
                    EngineGovernorConfig {enabled, setpoint_omega, PI, clamps}
                    EngineConfig {geometry, thermo, load, governor,
                        integration (numerics::IntegratorConfig), dt, max_steps,
                        initial θ/ω, history, csv_path}
engine_result.hpp   EngineState {theta_rad, omega, cycles}
                    EngineStepResult {t, dt_used, state, piston_x/v, p_cyl,
                        T_cyl, V_cyl, gas_force, indicated_moment,
                        load_moment, power, throttle}
                    EngineSimResult {valid, error, warnings, final_step,
                        total_time, total_indicated_work, mean power/omega,
                        cycles, mean_throttle, efficiency_estimate, history}
engine_simulator.hpp step_engine(state, t, config, governor?, status)
                    simulate_engine(config, status)
```

## 2. Mechanism (crank_mechanism.cpp)

Analytic slider-crank kinematics (no integration error):

```
x(θ)   = r·cosθ + √(l² − r²sin²θ)          (from crank axis)
dx/dθ, d²x/dθ² analytic;  J_eq(θ) = J_f + m_p·(dx/dθ)²
ODE:   θ̇ = ω
       J(θ)·ω̇ = M_gas(θ) − M_load(ω) − ½·(dJ/dθ)·ω²
```

The ½·(dJ/dθ)·ω² inertia-torque term is REQUIRED for energy consistency
(without it, KE = ½Jω² drifts non-physically; pinned by the pure-inertia
test that nulls the gas via a huge clearance volume).

## 3. Gas force / cycle models (gas_force.cpp) — pure functions of θ

No persisted cycle state → restartable, deterministic, batchable.

- **Otto** (θ mod 4π): polytropic compression (p·V^γc), power stroke with
  Wiebe heat release ANCHORED AT TDC VOLUME and decaying polytropically —
  `(γe−1)·q·x_b/V_tdc·(V_tdc/V)^γe`. (The raw `(γ−1)·q·x_b/V` form decays
  ∝1/V and over-produces work past the Otto bound — caught and fixed;
  measured η = 0.384 vs the 0.489 bound.) sin² valve ramps (no hard
  discontinuities for RK4/adaptive). Equal per-phase γ default; a
  validation warning fires when they differ (that difference performs net
  work per cycle with zero heat — a heat-transfer stand-in only).
- **Steam** (θ mod 2π, single-acting): saturated admission at p_boiler to
  cutoff, wet-steam polytrope `p = p_b·(V_cut/V)^n`, n = 1.13, exhaust to
  p_condenser. Dryness `x = x_cut·(V_cut/V)^(n−1)`; while x < 1 the
  cylinder temperature sits ON THE SATURATION LINE (`physics::thermo::saturation_temperature(p)`),
  otherwise ideal-gas with trapped mass `ρ_g(p_b)·V_cut·x_cut`.

## 4. Energy accounting (Rankine-lite bookkeeping)

- Otto: `heat_sum = Σ q_in·u·Δθ/4π` (exact per-step rollup of the
  throttled heat); `efficiency_estimate = W_ind / heat_sum`.
- Steam: per-cycle boiler heat `m·(h_g(T_sat,p_b) − h_f(T_sat,p_c))`
  accumulated over Δθ/2π; same efficiency metric.

## 5. Steam EOS (physics::thermo::steam)

Clausius–Clapeyron saturation model, single latent-heat anchor
(p₀ = 101.325 kPa at 100 °C, h_fg = 2.257 MJ/kg, R = 461.5):

```
p_sat(T)  = p₀·exp(−(h_fg/R)·(1/T − 1/T₀))
T_sat(p)  = inverse (closed form)
h_f, h_fg, h_g, h_wet(T, x), ρ_g(p, T)  — all closed-form
```

~5–10% of IAPWS-IF97 across 0–200 °C (documented engineering model;
swap `SteamConstants` for tables without touching the cycle).

## 6. Control & loads

- Governor: PI (physics::control::IController) updated EXACTLY ONCE per step; the
  resulting throttle (heat fraction / steam admission scale) is held
  constant inside the integrator so RK4 stages never mutate controller
  state. Loads are pure T(ω) functions (friction const+viscous, generator
  curve) — no per-step allocation.

## 7. Validation (tests/unit/engine/, tests/unit/thermo/steam_test.cpp)

| Check | Result |
|---|---|
| TDC/BDC kinematics, J_eq, finite-difference agreement | exact |
| Compression polytrope p·V^γc | exact |
| Otto TDC pressure p_intake·r_c^γc; Wiebe zero at ignition | exact / verified |
| Phase continuity (ramps) | < 2% |
| Pure-inertia free spin: KE conserved (½·dJ/dθ·ω² term) | drift < 5% over cycles |
| Fired Otto: η ∈ (0.2, 0.55), positive work, no NaN | ✓ |
| Governor regulation to setpoint | < 0.1% with mid-range throttle |
| Steam: admission T = T_sat(p_boiler), positive work, η ∈ (0.01, 0.45) | ✓ |
| Steam EOS: anchor round-trips, monotonicity, enthalpy consistency | ✓ |
| CSV streaming rows == steps + header | ✓ |

## 8. Limitations (documented)

- 0D cylinder (no in-cylinder mesh; moving-piston 3D CFD is Phase J).
- No blow-by, heat transfer to walls, knock, or valve-flow dynamics.
- Steam: single-constant latent heat; no superheat schedule; exhaust
  wetness not reported.
- Crank won't self-start from rest at exactly TDC (zero lever arm) —
  provide initial_omega (starter momentum), documented.
```
