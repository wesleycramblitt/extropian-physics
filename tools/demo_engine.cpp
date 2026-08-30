// demo_engine.cpp
// End-to-end real run: 4-stroke Otto engine with CSV motion output —
// crank angle, omega, piston x/v, cylinder p/T, torque, power per step.
//
//   build/bin/demo_engine [out.csv]   (default output/engine_state.csv)

#include <exd/physics/engine/engine_simulator.hpp>

#include <cstdio>
#include <string>

int main(int argc, char** argv)
{
    using namespace exd::physics::engine;
    const std::string out_csv = argc > 1 ? argv[1] : "output/engine_state.csv";

    EngineConfig cfg;
    cfg.geometry.crank_radius = 0.05;
    cfg.geometry.rod_length = 0.20;
    cfg.geometry.bore = 0.086;
    cfg.geometry.clearance_volume = 1.0e-4;
    cfg.geometry.piston_mass = 0.5;
    cfg.geometry.flywheel_inertia = 0.1;
    cfg.thermo.q_in_cycle = 1500.0;   // J heat release per cycle (Wiebe)
    cfg.initial_omega = 50.0;         // starter momentum
    cfg.dt = 2.0e-4;
    cfg.max_steps = 80000;            // 16 s
    cfg.csv_path = out_csv;
    // Speed governor: hold ~200 rad/s by throttling heat release.
    // (Efficiency uses the TIME-MEAN throttle, so the readout stays honest
    // under governing.)
    cfg.governor.enabled = true;
    cfg.governor.setpoint_omega = 200.0;
    cfg.governor.pi.kp = 3.0e-3;
    cfg.governor.pi.ki = 0.02;
    cfg.governor.pi.clamp_min = 0.0;
    cfg.governor.pi.clamp_max = 1.0;

    exd::physics::ModelStatus status;
    auto r = simulate_engine(cfg, status);
    if (!r.valid)
    {
        std::printf("engine run failed: %s\n", r.error.c_str());
        return 1;
    }
    std::printf("engine OK: %.3f s sim, %g cycles, mean indicated power %.1f W\n",
                r.total_time, r.cycles_completed, r.mean_indicated_power);
    std::printf("  final omega = %.2f rad/s (%.0f rpm), Otto efficiency est = %.3f\n",
                r.final_step.state.omega, r.final_step.state.omega * 60.0 / 6.283185307179586,
                r.efficiency_estimate);
    std::printf("output: %s (time,theta_rad,omega,piston_x,piston_v,p_cyl,T_cyl,...)\n",
                out_csv.c_str());
    return 0;
}
