# The engine (implementation_spec.md §57): one composable numerical runtime.
# Layers: core → numerics → discretization → physics → coupling → presets;
# backends are orthogonal. Geometry is the external `extropian-geometry` dep.

add_library(exd-physics STATIC

    # core (spec Phase 1 core runtime)
    src/core/convergence.cpp

    # coupling (spec Phase 9): channels, samplers, mapping, real exchange
    src/coupling/coupling_manager.cpp
    src/coupling/drag_body_solver.cpp
    src/coupling/field_channels.cpp
    src/coupling/field_sampler.cpp
    src/coupling/solver_manager.cpp
    src/coupling/surface_mapping.cpp

    # numerics (spec Phase 5): ODE integrators + time stepping
    src/numerics/integrators.cpp
    src/numerics/time_stepping.cpp

    # output (spec §43): field stamps, CSV series
    src/output/csv_series.cpp
    src/output/field_writer.cpp

    # physics — acoustics
    src/physics/acoustics/wave_solver.cpp

    # physics — control
    src/physics/control/pi_controller.cpp

    # physics — electromagnetics (static fields, FDTD, lumped machines)
    src/physics/electromagnetics/circuit.cpp
    src/physics/electromagnetics/fdtd.cpp
    src/physics/electromagnetics/static_fields.cpp

    # physics — fluid: FDM 2D (legacy), FDM3, forces, BEM, lumped, turbomachinery
    src/physics/fluid/fdm3/fdm3_boundary.cpp
    src/physics/fluid/fdm3/fdm3_grid.cpp
    src/physics/fluid/fdm3/fdm3_integration.cpp
    src/physics/fluid/fdm3/fdm3_pressure.cpp
    src/physics/fluid/fdm3/fdm3_solver.cpp
    src/physics/fluid/fdm3/fdm3_spatial.cpp
    src/physics/fluid/fdm/fdm_boundary.cpp
    src/physics/fluid/fdm/fdm_grid.cpp
    src/physics/fluid/fdm/fdm_integration.cpp
    src/physics/fluid/fdm/fdm_pressure.cpp
    src/physics/fluid/fdm/fdm_solver.cpp
    src/physics/fluid/fdm/fdm_spatial.cpp
    src/physics/fluid/forces/blade_element.cpp
    src/physics/fluid/forces/blade_surface.cpp
    src/physics/fluid/forces/force_evaluator.cpp
    src/physics/fluid/forces/momentum_balance.cpp
    src/physics/fluid/forces/pressure_integration.cpp
    src/physics/fluid/forces/table_lookup.cpp
    src/physics/fluid/lumped/plenum.cpp
    src/physics/fluid/reduced_order/bem/bem_solver.cpp
    src/physics/fluid/reduced_order/bem/blade_geometry.cpp
    src/physics/fluid/reduced_order/bem/duct_flow.cpp
    src/physics/fluid/reduced_order/bem/flow_field.cpp
    src/physics/fluid/reduced_order/bem/hull_drag.cpp
    src/physics/fluid/reduced_order/bem/induction_corrections.cpp
    src/physics/fluid/reduced_order/bem/induction.cpp
    src/physics/fluid/reduced_order/bem/loss_corrections.cpp
    src/physics/fluid/reduced_order/bem/losses.cpp
    src/physics/fluid/reduced_order/bem/polar_db.cpp
    src/physics/fluid/turbomachinery/compression_system.cpp
    src/physics/fluid/turbomachinery/operating_map.cpp
    src/physics/fluid/turbomachinery/stage.cpp
    src/physics/fluid/turbomachinery/stage_stack.cpp

    # physics — particles, reaction, species transport, porous media, rigid bodies, structural, thermal, thermo
    src/physics/particles/particle_track.cpp
    src/physics/species/species_solver.cpp
    src/physics/porous/porous_solver.cpp
    src/physics/reaction/reactor.cpp
    src/physics/rigid_body/moment_integration.cpp
    src/physics/rigid_body/moment_models.cpp
    src/physics/rigid_body/quaternion.cpp
    src/physics/rigid_body/rigid_body.cpp
    src/physics/rigid_body/rigid_dynamics.cpp
    src/physics/rigid_body/rotating_assembly.cpp
    src/physics/structural/elasticity.cpp
    src/physics/thermal/thermal_solver.cpp
    src/physics/thermo/eos.cpp
    src/physics/thermo/polytropic.cpp
    src/physics/thermo/properties.cpp

    # presets (spec §14/§69): validated assemblies — NOT solver implementations
    src/presets/multiphysics/aeroacoustics.cpp
    src/presets/multiphysics/joule_heating.cpp
    src/presets/multiphysics/species_in_flow.cpp
    src/presets/multiphysics/thermal_stress.cpp
    src/presets/engine/crank_mechanism.cpp
    src/presets/engine/engine_simulator.cpp
    src/presets/engine/gas_force.cpp
    src/presets/turbine/coupled_turbine.cpp
    src/presets/turbine/turbine_builder.cpp
    src/presets/turbine/turbine_simulator.cpp
)

# Hot kernel loops compile optimized even in developer builds (default
# CMAKE_BUILD_TYPE is intentionally empty in this repo).
set_source_files_properties(
    src/physics/fluid/fdm3/fdm3_grid.cpp
    src/physics/fluid/fdm3/fdm3_spatial.cpp
    src/physics/fluid/fdm3/fdm3_integration.cpp
    src/physics/fluid/fdm3/fdm3_pressure.cpp
    src/physics/fluid/fdm3/fdm3_boundary.cpp
    src/physics/fluid/fdm3/fdm3_solver.cpp
    src/physics/fluid/turbomachinery/stage.cpp
    src/physics/fluid/turbomachinery/stage_stack.cpp
    src/physics/fluid/turbomachinery/operating_map.cpp
    src/physics/fluid/turbomachinery/compression_system.cpp
    src/physics/fluid/lumped/plenum.cpp
    src/physics/thermal/thermal_solver.cpp
    src/physics/acoustics/wave_solver.cpp
    src/physics/structural/elasticity.cpp
    src/physics/particles/particle_track.cpp
    src/physics/reaction/reactor.cpp
    PROPERTIES COMPILE_OPTIONS "-O2"
)

target_include_directories(exd-physics PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_include_directories(exd-physics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(exd-physics PUBLIC
    exd::core
    exd::geometry
    Threads::Threads
)

add_library(exd::physics ALIAS exd-physics)

# ── Install (deferred; re-enable for packaging) ──
#include(GNUInstallDirs)
#install(TARGETS exd-physics EXPORT exd-physics-targets
#    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#)
#install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(EXPORT exd-physics-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exd-physics)
