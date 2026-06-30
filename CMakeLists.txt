cmake_minimum_required(VERSION 3.21)
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

    # Coupling
    src/coupling/surface_mapping.cpp
    src/coupling/coupling_manager.cpp
)

target_include_directories(exd-physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(exd-physics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(exd-physics PUBLIC
    exd::core
    Threads::Threads
)

add_library(exd::physics ALIAS exd-physics)

# ── Install ───────────────────────────────────────────
include(GNUInstallDirs)
install(TARGETS exd-physics EXPORT exd-physics-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT exd-physics-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exd-physics)

# ── Tests ─────────────────────────────────────────────
option(EXT_PHYSICS_BUILD_TESTS "Build tests" OFF)
if(EXT_PHYSICS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
