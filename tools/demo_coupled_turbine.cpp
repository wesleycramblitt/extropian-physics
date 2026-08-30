// demo_coupled_turbine.cpp
// End-to-end real run: parametric turbine inside the 3D FDM with field
// stamps (exd-fld v1) + rotor machine-state CSV — the exact outputs the
// animation/visualization repo consumes (docs/output_channels.md).
//
//   build/bin/demo_coupled_turbine [out_dir]
//   → out_dir/turbine_field/ (step_*.fld + timeline.txt)
//   → out_dir/turbine_rotor.csv

#include <exd/physics/io/field_writer.hpp>
#include <exd/physics/io/output_policy.hpp>
#include <exd/physics/turbine/coupled_turbine.hpp>
#include <exd/physics/turbine/turbine_builder.hpp>

#include <cstdio>
#include <filesystem>
#include <string>

int main(int argc, char** argv)
{
    const std::string out_dir = argc > 1 ? argv[1] : "output";

    // 1. Build the turbine from engineering parameters.
    exd::physics::turbine::TurbineBuilderConfig tb;
    tb.hub_radius = 0.15;
    tb.tip_radius = 0.5;
    tb.chord = 0.12;
    tb.twist_hub_deg = 0.0;
    tb.twist_tip_deg = -3.0;
    tb.rpm = 60.0;
    tb.blade_count = 3;
    tb.section_count = 4;
    tb.leading_edge_z = 0.3;
    tb.duct_length = 2.0;

    exd::physics::ModelStatus status;
    const auto turbine = exd::physics::turbine::make_turbine_definition(tb, status);
    if (!status.ok)
    {
        std::printf("turbine builder: %s\n", status.error.c_str());
        return 1;
    }

    // 2. Coupled run configuration.
    using namespace exd::physics::turbine;
    CoupledTurbineConfig c;
    c.turbine = turbine;
    const double v_inf = 3.0;
    c.grid = default_grid_config(v_inf, 20, 3.0, 4.0); // 3R × 3R × 4R box
    c.rotor_origin = {1.5, 1.5, 1.5};
    c.element_count = 12;
    c.rotor_inertia = 0.1;
    c.fluid_steps_per_exchange = 10;
    c.force_relaxation = 0.4;
    c.ramp_time_s = 2.0;
    c.dt = 0.0;                     // grid dt (CFL-based default)
    c.max_steps = 2000;
    c.record_history = true;
    c.history_interval = 2;
    c.csv_path = out_dir + "/turbine_rotor.csv";

    // 3. Output: binary field stamps throttled to a cadence + CSV states.
    const std::string field_dir = out_dir + "/turbine_field";
    exd::physics::io::FldWriterConfig fw;
    fw.directory = field_dir;
    auto writer = exd::physics::io::make_fld_writer(fw, status);
    if (!writer)
    {
        std::printf("field writer: %s\n", status.error.c_str());
        return 1;
    }
    exd::physics::io::OutputScheduler sched({/*every_n_steps=*/50, /*wall_clock_s=*/0.0});

    // 4. Wire output into the run, then go.
    c.field_writer = writer.get();        // velocity + pressure stamps
    c.output_scheduler = &sched;          // every 50 fluid steps
    // (wall-clock real-time: OutputPolicy{0, 0.25} instead)
    auto r = run_coupled_turbine(c, status);
    if (!r.valid)
    {
        std::printf("coupled run failed: %s\n", r.error.c_str());
        return 1;
    }

    // 6. Summary.
    std::printf("coupled run OK: %llu fluid steps, %llu exchanges\n",
                static_cast<unsigned long long>(r.fluid.steps_taken),
                static_cast<unsigned long long>(r.exchanges));
    std::printf("  final omega = %.4f rad/s (%.1f rpm), tsr = %.3f, cp = %.4f\n",
                r.final_omega, r.final_omega * 60.0 / 6.283185307179586, r.final_tsr, r.final_cp);
    std::printf("  aero_work = %.4f J, rotor dKE = %.4f J, load_work = %.4f J\n",
                r.aero_work, r.rotor_ke_change, r.load_work);
    std::printf("outputs:\n  fields : %s/ (timeline.txt + step_*.fld)\n  rotor  : %s\n",
                field_dir.c_str(), c.csv_path.c_str());
    return 0;
}
