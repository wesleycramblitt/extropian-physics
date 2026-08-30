cmake_minimum_required(VERSION 3.21)

# Allow fetched dependencies with old cmake_minimum_required() to configure
# under CMake >= 4.0 without extra command-line flags.
if(POLICY CMP0167)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Minimum policy version for fetched dependencies" FORCE)
endif()

project(extropian-physics VERSION 0.1.0 LANGUAGES CXX
    DESCRIPTION "Extropian Physics: mesh, fields, boundary conditions, solver plugin interface, coupling"
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Dependencies: use FetchContent when not already available ──
include(FetchContent)

if(NOT TARGET exd::core)
    FetchContent_Declare(exd-core
        GIT_REPOSITORY https://github.com/wesleycramblitt/extropian-core.git
        GIT_TAG main
    )
    FetchContent_MakeAvailable(exd-core)
endif()

if(NOT TARGET exd::geometry)
    if(DEFINED EXD_GEOMETRY_DIR)
        FetchContent_Declare(exd-geometry SOURCE_DIR "${EXD_GEOMETRY_DIR}")
    else()
        FetchContent_Declare(exd-geometry
            GIT_REPOSITORY https://github.com/wesleycramblitt/extropian-geometry.git
            GIT_TAG main
        )
    endif()
    if(DEFINED EXD_CORE_DIR)
        FetchContent_Declare(exd-core SOURCE_DIR "${EXD_CORE_DIR}")
    endif()
    set(ENABLE_TEXT  OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTS  OFF CACHE BOOL "" FORCE)
    set(BUILD_DEMO   OFF CACHE BOOL "" FORCE)
    set(EXT_CORE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(exd-geometry)
endif()

find_package(Threads REQUIRED)  # For solver threading

# ── Library target ────────────────────────────────────
# NOTE: This library contains the INTERFACE and DATA MODEL only.
# Solver implementations are in separate extropian-solver-* repos.

add_library(exd-physics STATIC
    # Mesh
    src/mesh/unstructured_volume.cpp
    src/mesh/surface.cpp
    src/mesh/structured.cpp
    src/mesh/quality.cpp
    src/mesh/partitioning.cpp
    src/mesh/io/vtk.cpp
    src/mesh/io/gmsh.cpp
    src/mesh/io/cgns.cpp
    src/mesh/io/stl.cpp
    src/mesh/io/obj.cpp
    src/mesh/io/openfoam.cpp

    # Fields
    src/field/scalar.cpp
    src/field/vector.cpp
    src/field/tensor.cpp
    src/field/interpolation.cpp

    # Boundary conditions
    src/bc/framework.cpp
    src/bc/fluid.cpp
    src/bc/solid.cpp
    src/bc/serialization.cpp

    # Materials
    src/material/material.cpp
    src/material/database.cpp

    # Solver interface (no implementations — those are in extropian-solver-*)
    src/solver/solver_manager.cpp
    src/solver/time_stepping.cpp
    src/solver/convergence.cpp
    src/solver/integrators.cpp

    # Coupling
    src/coupling/surface_mapping.cpp
    src/coupling/coupling_manager.cpp

    # Level-3 BEM turbine solver
    src/fluid/reduced_order/bem/bem_solver.cpp
    src/fluid/reduced_order/bem/blade_geometry.cpp
    src/fluid/reduced_order/bem/induction.cpp
    src/fluid/reduced_order/bem/losses.cpp
    src/fluid/reduced_order/bem/polar_db.cpp
    src/fluid/reduced_order/bem/duct_flow.cpp
    src/fluid/reduced_order/bem/hull_drag.cpp
    src/fluid/reduced_order/bem/flow_field.cpp
    src/fluid/reduced_order/bem/induction_corrections.cpp
    src/fluid/reduced_order/bem/loss_corrections.cpp
    # Level-1 FDM CFD solver
    src/fluid/fdm/fdm_grid.cpp
    src/fluid/fdm/fdm_spatial.cpp
    src/fluid/fdm/fdm_integration.cpp
    src/fluid/fdm/fdm_pressure.cpp
    src/fluid/fdm/fdm_boundary.cpp
    src/fluid/fdm/fdm_solver.cpp

    # Level-1 3D FDM CFD solver (fdm3)
    src/fluid/fdm3/fdm3_grid.cpp
    src/fluid/fdm3/fdm3_spatial.cpp
    src/fluid/fdm3/fdm3_integration.cpp
    src/fluid/fdm3/fdm3_pressure.cpp
    src/fluid/fdm3/fdm3_boundary.cpp
    src/fluid/fdm3/fdm3_solver.cpp

    # Generic mechanics (rigid-body rotational dynamics, moments, assembly)
    src/mechanics/moment_integration.cpp
    src/mechanics/moment_models.cpp
    src/mechanics/rigid_dynamics.cpp
    src/mechanics/rotating_assembly.cpp
    src/mechanics/quaternion.cpp
    src/mechanics/rigid_body.cpp

    # Generic fluid force evaluation (FSI coupling point)
    src/fluid/forces/blade_surface.cpp
    src/fluid/forces/pressure_integration.cpp
    src/fluid/forces/momentum_balance.cpp
    src/fluid/forces/table_lookup.cpp
    src/fluid/forces/blade_element.cpp
    src/fluid/forces/force_evaluator.cpp

    # Coupling field sampling (CFD field -> surface flow)
    src/coupling/field_sampler.cpp
    src/coupling/field_channels.cpp

    # Turbine application assembly
    src/turbine/turbine_simulator.cpp
    src/turbine/coupled_turbine.cpp

    # Thermo (EOS + properties)
    src/thermo/eos.cpp
    src/thermo/properties.cpp

    # Control (PI/pitch/governor)
    src/control/pi_controller.cpp

    # Electrical (static fields, FDTD, lumped machines)
    src/electrical/static_fields.cpp
    src/electrical/fdtd.cpp
    src/electrical/circuit.cpp

    # Engine application (slider-crank + cycle + governor)
    src/engine/crank_mechanism.cpp
    src/engine/gas_force.cpp
    src/engine/engine_simulator.cpp

    # Turbine builder (parametric TurbineDefinition)"
    src/turbine/turbine_builder.cpp

    # I/O output channels (field stamps, CSV machine-state series)
    src/io/field_writer.cpp
    src/io/csv_series.cpp
)

# The fdm3 kernels are hot loops over (nx+2)(ny+2)(nz+2) cells; compile
# them optimized even in developer builds (default CMAKE_BUILD_TYPE is
# intentionally empty in this repo).
set_source_files_properties(
    src/fluid/fdm3/fdm3_grid.cpp
    src/fluid/fdm3/fdm3_spatial.cpp
    src/fluid/fdm3/fdm3_integration.cpp
    src/fluid/fdm3/fdm3_pressure.cpp
    src/fluid/fdm3/fdm3_boundary.cpp
    src/fluid/fdm3/fdm3_solver.cpp
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

# ── Install ───────────────────────────────────────────
#include(GNUInstallDirs)
#install(TARGETS exd-physics EXPORT exd-physics-targets
#    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#)
#install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(EXPORT exd-physics-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exd-physics)

# ── Tests ─────────────────────────────────────────────
option(EXT_PHYSICS_BUILD_TESTS "Build tests" OFF)

if(EXT_PHYSICS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
