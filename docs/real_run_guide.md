# Real Runs — From Zero to Output Fields

How to actually run a simulation and get output fields, end to end.
Companion to `docs/output_channels.md` (the parsing contract).

---

## 1. Two runnable demos (fastest path)

```bash
cmake -S . -B build -DEXT_PHYSICS_BUILD_TESTS=ON
cmake --build build -j
./build/demo_engine                # → output/engine_state.csv
./build/demo_coupled_turbine out   # → out/turbine_field/*.fld + timeline.txt
                                   #   and out/turbine_rotor.csv
```

- `demo_engine` — governed 4-stroke Otto: crank angle, omega, piston x/v,
  p/T, torque, power per step.
- `demo_coupled_turbine` — parametric rotor inside the 3D FDM:
  velocity/pressure field stamps + rotor machine-state CSV.

## 2. Rolling your own run (what the demos do)

### Turbine + 3D CFD with fields

```cpp
#include <exd/physics/turbine/turbine_builder.hpp>   // no hand-built fixtures
#include <exd/physics/turbine/coupled_turbine.hpp>
#include <exd/physics/io/field_writer.hpp>

exd::physics::ModelStatus status;

// 1. Turbine from engineering parameters.
exd::physics::turbine::TurbineBuilderConfig tb;   // hub/tip radius, chord,
                                                  // twist hub→tip, rpm, blades
auto turbine = exd::physics::turbine::make_turbine_definition(tb, status);

// 2. Grid + run config. default_grid_config(v_inf, n) builds the box with
//    the −Z inflow convention (Inlet at +z face, Outlet at −z, lateral symmetry).
exd::physics::turbine::CoupledTurbineConfig c;
c.turbine = turbine;
c.grid = exd::physics::turbine::default_grid_config(3.0, 20);
c.rotor_origin = {1.5, 1.5, 1.5};     // axis point in grid coordinates
c.rotor_inertia = 0.1;
c.fluid_steps_per_exchange = 10;
c.max_steps = 2000;
c.csv_path = "out/turbine_rotor.csv"; // rotor states, one flushed row/step

// 3. Field stamps at a cadence (real-time: wall-clock throttle).
exd::physics::io::FldWriterConfig fw; fw.directory = "out/turbine_field";
auto writer = exd::physics::io::make_fld_writer(fw, status);
exd::physics::io::OutputScheduler sched({/*every_n_steps=*/50, /*wall_clock_s=*/0.0});
// → or {0, 0.25} for ~4 fps wall-clock real-time output

auto r = exd::physics::turbine::run_coupled_turbine(c, status);
// r.history: t, omega, angle, torque, thrust, power, exchange;
// r.final_omega/cp/tsr; r.aero_work/rotor_ke_change/load_work
```

### Engine with motion output

```cpp
#include <exd/physics/engine/engine_simulator.hpp>
exd::physics::engine::EngineConfig cfg;
cfg.thermo.q_in_cycle = 1500.0;      // J per cycle (Wiebe)
cfg.initial_omega = 50.0;            // starter momentum
cfg.governor.enabled = true;         // PI holds setpoint_omega
cfg.governor.setpoint_omega = 200.0;
cfg.csv_path = "out/engine_state.csv";
auto r = exd::physics::engine::simulate_engine(cfg, status);
// CSV: time,theta_rad,omega_rad_s,piston_x_m,piston_v_m_s,p_cyl_pa,T_cyl_K,
//      indicated_moment_Nm,load_moment_Nm,power_W,throttle,cycles
```

### Standalone 3D CFD case with fields

```cpp
#include <exd/physics/fluid/fdm3/fdm3_solver.hpp>
auto cfg = ...;                        // FDM3Config with your BCs
auto r = exd::physics::fluid::fdm3::run_fdm3_simulation(cfg, writer, &sched, &status);
// velocity + pressure stamps at the cadence; without a writer this is
// equivalent to solve_fdm3(config) (pure, optimizer-batchable)
```

## 3. What you get back

| Output | Location | Contents |
|---|---|---|
| Field stamps | `out/turbine_field/step_*.fld` | header + float32 velocity (interleaved xyz) and pressure; per-field origin/spacing/dims; cell-centered convention |
| Timeline | `out/turbine_field/timeline.txt` | `time step filename` per stamp — the animation clock |
| Rotor states | `out/turbine_rotor.csv` | t, ω, θ, torque, thrust, power, exchange |
| Engine kinematics | `out/engine_state.csv` | θ, ω, piston x/v, p/T, torque, power, throttle, cycles |

Parse spec: `docs/output_channels.md`. The animation repo needs only a
~150-line exd-fld reader (fixed header, raw float32 payloads, little-endian).

## 4. Notes for real use

- **Grid size**: 20³–40³ boxes run in seconds on one core; 48³+ gets
  minutes-per-run. `sor_omega = 1.5–1.7` and warm-started pressure help.
- **Coupling stability**: keep `dt < ε/(3·V_inf)` (ε = smear_cells·dz —
  enforced by validation), `ramp_time_s ≥ 10·window`, `force_relaxation`
  ≈ 0.3–0.5. These are the guards from `docs/modular_solver_architecture.md`.
- **Fidelity caveat**: the collocated scheme is engineering-grade, not
  higher-order RANS — use it for design-space exploration, coupling
  studies, and animation; push fidelity toward FVM/FEM (Phase J) for
  publication-grade numbers.
- **Optimizer loop**: `solve_*`/`simulate_*` take no writer parameters;
  batches run with null sinks (zero I/O overhead).
