# CMakeLists.txt — doctest suites per engine layer (spec §63-64 verification).
# Test names are stable; source paths follow the spec §57 layout.

include(FetchContent)

if(NOT TARGET doctest::doctest_with_main)
    FetchContent_Declare(doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG v2.4.11
    )
    FetchContent_MakeAvailable(doctest)
endif()

set(TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/data")

function(add_physics_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE exd::physics doctest::doctest_with_main)
    target_compile_definitions(${name} PRIVATE EXT_PHYSICS_DATA_DIR="${TEST_DATA_DIR}")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_physics_test(polar_test unit/physics/fluid/bem/polar_test.cpp)
add_physics_test(blade_geometry_test unit/physics/fluid/bem/blade_geometry_test.cpp)
add_physics_test(bem_betz_test unit/physics/fluid/bem/bem_betz_test.cpp)
add_physics_test(induction_test unit/physics/fluid/bem/induction_test.cpp)
add_physics_test(duct_test unit/physics/fluid/bem/duct_test.cpp)
add_physics_test(flow_field_test unit/physics/fluid/bem/flow_field_test.cpp)
add_physics_test(hull_test unit/physics/fluid/bem/hull_test.cpp)
add_physics_test(integration_test unit/physics/fluid/bem/integration_test.cpp)
add_physics_test(induction_corrections_test unit/physics/fluid/bem/induction_corrections_test.cpp)
add_physics_test(loss_corrections_test unit/physics/fluid/bem/loss_corrections_test.cpp)
add_physics_test(fdm_grid_test unit/physics/fluid/fdm/fdm_grid_test.cpp)
add_physics_test(fdm_integration_test unit/physics/fluid/fdm/fdm_integration_test.cpp)
add_physics_test(fdm_pressure_test unit/physics/fluid/fdm/fdm_pressure_test.cpp)
add_physics_test(fdm_solver_test unit/physics/fluid/fdm/fdm_solver_test.cpp)
add_physics_test(fdm3_grid_test unit/physics/fluid/fdm3/fdm3_grid_test.cpp)
add_physics_test(fdm3_pressure_test unit/physics/fluid/fdm3/fdm3_pressure_test.cpp)
add_physics_test(fdm3_solver_test unit/physics/fluid/fdm3/fdm3_solver_test.cpp)
add_physics_test(fdm3_body_force_test unit/physics/fluid/fdm3/fdm3_body_force_test.cpp)
add_physics_test(integrator_test unit/numerics/integrator_test.cpp)
add_physics_test(time_stepping_test unit/numerics/time_stepping_test.cpp)
add_physics_test(moment_integration_test unit/physics/rigid_body/moment_integration_test.cpp)
add_physics_test(quaternion_test unit/physics/rigid_body/quaternion_test.cpp)
add_physics_test(rigid_body_test unit/physics/rigid_body/rigid_body_test.cpp)
add_physics_test(moment_models_test unit/physics/rigid_body/moment_models_test.cpp)
add_physics_test(rigid_dynamics_test unit/physics/rigid_body/rigid_dynamics_test.cpp)
add_physics_test(rotating_assembly_test unit/physics/rigid_body/rotating_assembly_test.cpp)
add_physics_test(blade_surface_test unit/physics/fluid/forces/blade_surface_test.cpp)
add_physics_test(pressure_integration_test unit/physics/fluid/forces/pressure_integration_test.cpp)
add_physics_test(momentum_balance_test unit/physics/fluid/forces/momentum_balance_test.cpp)
add_physics_test(field_sampler_test unit/coupling/field_sampler_test.cpp)
add_physics_test(turbine_simulator_test unit/presets/turbine/turbine_simulator_test.cpp)
add_physics_test(coupled_turbine_test unit/presets/turbine/coupled_turbine_test.cpp)
add_physics_test(eos_test unit/physics/thermo/eos_test.cpp)
add_physics_test(pi_controller_test unit/physics/control/pi_controller_test.cpp)
add_physics_test(field_channels_test unit/coupling/field_channels_test.cpp)
add_physics_test(static_fields_test unit/physics/electromagnetics/static_fields_test.cpp)
add_physics_test(fdtd_test unit/physics/electromagnetics/fdtd_test.cpp)
add_physics_test(circuit_test unit/physics/electromagnetics/circuit_test.cpp)
add_physics_test(binary_fld_writer_test unit/output/binary_fld_writer_test.cpp)
add_physics_test(csv_series_test unit/output/csv_series_test.cpp)
add_physics_test(output_policy_test unit/output/output_policy_test.cpp)
add_physics_test(crank_mechanism_test unit/presets/engine/crank_mechanism_test.cpp)
add_physics_test(gas_force_test unit/presets/engine/gas_force_test.cpp)
add_physics_test(engine_simulator_test unit/presets/engine/engine_simulator_test.cpp)
add_physics_test(turbine_builder_test unit/presets/turbine/turbine_builder_test.cpp)
add_physics_test(fdm3_output_driver_test unit/physics/fluid/fdm3/fdm3_output_driver_test.cpp)
add_physics_test(steam_test unit/physics/thermo/steam_test.cpp)
add_physics_test(polytropic_test unit/physics/thermo/polytropic_test.cpp)
add_physics_test(stage_test unit/physics/fluid/turbomachinery/stage_test.cpp)
add_physics_test(stage_stack_test unit/physics/fluid/turbomachinery/stage_stack_test.cpp)
add_physics_test(operating_map_test unit/physics/fluid/turbomachinery/operating_map_test.cpp)
add_physics_test(plenum_test unit/physics/fluid/lumped/plenum_test.cpp)
add_physics_test(compression_system_test unit/physics/fluid/turbomachinery/compression_system_test.cpp)
add_physics_test(surface_mapping_test unit/coupling/surface_mapping_test.cpp)
add_physics_test(coupling_manager_test unit/coupling/coupling_manager_test.cpp)
add_physics_test(drag_body_fsi_test unit/coupling/drag_body_fsi_test.cpp)
add_physics_test(thermal_solver_test unit/physics/thermal/thermal_solver_test.cpp)
add_physics_test(thermal_coupling_test unit/physics/thermal/thermal_coupling_test.cpp)
add_physics_test(wave_solver_test unit/physics/acoustics/wave_solver_test.cpp)
add_physics_test(elasticity_test unit/physics/structural/elasticity_test.cpp)
add_physics_test(particle_track_test unit/physics/particles/particle_track_test.cpp)
add_physics_test(reactor_test unit/physics/reaction/reactor_test.cpp)

# ── Spec-conformance core runtime (Phase 1) ──
add_physics_test(units_test unit/core/units_test.cpp)
add_physics_test(field_test unit/core/field_test.cpp)
add_physics_test(execution_graph_test unit/core/execution_test.cpp)

# ── Mesh layer (Phase 2) ──
add_physics_test(mesh_test unit/mesh/mesh_test.cpp)

# ── FDM operators (Phase 3) ──
add_physics_test(fdm_operators_test unit/discretization/fdm_operators_test.cpp)

# ── Linear + nonlinear solvers (Phase 5) ──
add_physics_test(linear_solvers_test unit/numerics/linear_solvers_test.cpp)
add_physics_test(nonlinear_test unit/numerics/nonlinear_test.cpp)

# ── Coupling contracts + compatibility rules (Phase 9) ──
add_physics_test(contracts_test unit/coupling_contracts/contracts_test.cpp)

# ── Fidelity + diagnostics ──
add_physics_test(fidelity_test unit/fidelity/fidelity_test.cpp)

# ── Presets (Phase 10) ──
add_physics_test(cht_preset_test unit/presets/cht_preset_test.cpp)
add_physics_test(cfd_preset_test unit/presets/cfd_preset_test.cpp)
add_physics_test(multiphysics_presets_test unit/presets/multiphysics_presets_test.cpp)
add_physics_test(steady_conduction_preset_test unit/presets/steady_conduction_preset_test.cpp)
add_physics_test(robotics_trajectory_test unit/presets/robotics_trajectory_test.cpp)

# ── Species transport (W14: FDM domain breadth) ──
add_physics_test(species_test unit/physics/species/species_test.cpp)
add_physics_test(porous_test unit/physics/porous/porous_test.cpp)
add_physics_test(manipulator_test unit/physics/robotics/manipulator_test.cpp)
add_physics_test(heterogeneous_test unit/physics/thermal/heterogeneous_test.cpp)

# ── Transient structural / elastic waves (W14) ──
add_physics_test(elastic_wave_test unit/physics/structural/elastic_wave_test.cpp)
