# tests/CMakeLists.txt

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

add_physics_test(polar_test unit/bem/polar_test.cpp)
add_physics_test(blade_geometry_test unit/bem/blade_geometry_test.cpp)
add_physics_test(bem_betz_test unit/bem/bem_betz_test.cpp)
add_physics_test(induction_test unit/bem/induction_test.cpp)
add_physics_test(duct_test unit/bem/duct_test.cpp)
add_physics_test(flow_field_test unit/bem/flow_field_test.cpp)
add_physics_test(hull_test unit/bem/hull_test.cpp)
add_physics_test(integration_test unit/bem/integration_test.cpp)

# Correction model tests.
add_physics_test(induction_corrections_test unit/bem/induction_corrections_test.cpp)
add_physics_test(loss_corrections_test unit/bem/loss_corrections_test.cpp)

# FDM CFD solver tests.
add_physics_test(fdm_grid_test unit/fdm/fdm_grid_test.cpp)
add_physics_test(fdm_integration_test unit/fdm/fdm_integration_test.cpp)
add_physics_test(fdm_pressure_test unit/fdm/fdm_pressure_test.cpp)
add_physics_test(fdm_solver_test unit/fdm/fdm_solver_test.cpp)

# 3D FDM CFD solver (fdm3) tests.
add_physics_test(fdm3_grid_test unit/fdm3/fdm3_grid_test.cpp)
add_physics_test(fdm3_pressure_test unit/fdm3/fdm3_pressure_test.cpp)
add_physics_test(fdm3_solver_test unit/fdm3/fdm3_solver_test.cpp)
add_physics_test(fdm3_body_force_test unit/fdm3/fdm3_body_force_test.cpp)

# Shared integrators + time stepping (Phase A)
add_physics_test(integrator_test unit/solver/integrator_test.cpp)
add_physics_test(time_stepping_test unit/solver/time_stepping_test.cpp)

# Generic mechanics
add_physics_test(moment_integration_test unit/mechanics/moment_integration_test.cpp)
add_physics_test(quaternion_test unit/mechanics/quaternion_test.cpp)
add_physics_test(rigid_body_test unit/mechanics/rigid_body_test.cpp)
add_physics_test(moment_models_test unit/mechanics/moment_models_test.cpp)
add_physics_test(rigid_dynamics_test unit/mechanics/rigid_dynamics_test.cpp)
add_physics_test(rotating_assembly_test unit/mechanics/rotating_assembly_test.cpp)

# Generic fluid force evaluation
add_physics_test(blade_surface_test unit/fluid/forces/blade_surface_test.cpp)
add_physics_test(pressure_integration_test unit/fluid/forces/pressure_integration_test.cpp)
add_physics_test(momentum_balance_test unit/fluid/forces/momentum_balance_test.cpp)

# Coupling field sampling
add_physics_test(field_sampler_test unit/coupling/field_sampler_test.cpp)

# Turbine application assembly
add_physics_test(turbine_simulator_test unit/turbine/turbine_simulator_test.cpp)
add_physics_test(coupled_turbine_test unit/turbine/coupled_turbine_test.cpp)

# Thermo + control (Phase C)
add_physics_test(eos_test unit/thermo/eos_test.cpp)
add_physics_test(pi_controller_test unit/control/pi_controller_test.cpp)

# Coupling field channels (Phase D)
add_physics_test(field_channels_test unit/coupling/field_channels_test.cpp)

# Electrical (Phase D)
add_physics_test(static_fields_test unit/electrical/static_fields_test.cpp)
add_physics_test(fdtd_test unit/electrical/fdtd_test.cpp)
add_physics_test(circuit_test unit/electrical/circuit_test.cpp)

# I/O output channel tests.
add_physics_test(binary_fld_writer_test unit/io/binary_fld_writer_test.cpp)
add_physics_test(csv_series_test unit/io/csv_series_test.cpp)
add_physics_test(output_policy_test unit/io/output_policy_test.cpp)

# Engine tests.
add_physics_test(crank_mechanism_test unit/engine/crank_mechanism_test.cpp)
add_physics_test(gas_force_test unit/engine/gas_force_test.cpp)
add_physics_test(engine_simulator_test unit/engine/engine_simulator_test.cpp)

# Turbine builder + fdm3 stamping-driver tests (Wave 7).
add_physics_test(turbine_builder_test unit/turbine/turbine_builder_test.cpp)
add_physics_test(fdm3_output_driver_test unit/fdm3/fdm3_output_driver_test.cpp)

# Steam properties + water turbine capability (Wave 8).
add_physics_test(steam_test unit/thermo/steam_test.cpp)

# Stagnation polytrope primitives + EOS density/entropy (W9).
add_physics_test(polytropic_test unit/thermo/polytropic_test.cpp)

# Turbomachinery generalization (W9): mean-line stage/stack/map
add_physics_test(stage_test unit/fluid/turbomachinery/stage_test.cpp)
add_physics_test(stage_stack_test unit/fluid/turbomachinery/stage_stack_test.cpp)
add_physics_test(operating_map_test unit/fluid/turbomachinery/operating_map_test.cpp)

# Greitzer plenum + compression system driver (W9)
add_physics_test(plenum_test unit/fluid/lumped/plenum_test.cpp)
add_physics_test(compression_system_test unit/fluid/turbomachinery/compression_system_test.cpp)

# Real coupling exchange (W10, Phase H-lite)
add_physics_test(surface_mapping_test unit/coupling/surface_mapping_test.cpp)
add_physics_test(coupling_manager_test unit/coupling/coupling_manager_test.cpp)

# Phase I grid-first domains (W10)
add_physics_test(thermal_solver_test unit/thermal/thermal_solver_test.cpp)
add_physics_test(thermal_coupling_test unit/thermal/thermal_coupling_test.cpp)
add_physics_test(wave_solver_test unit/acoustics/wave_solver_test.cpp)
add_physics_test(elasticity_test unit/structural/elasticity_test.cpp)
add_physics_test(particle_track_test unit/particles/particle_track_test.cpp)
add_physics_test(reactor_test unit/chemistry/reactor_test.cpp)
