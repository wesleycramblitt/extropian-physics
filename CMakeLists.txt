cmake_minimum_required(VERSION 3.21)

# Allow fetched dependencies with old cmake_minimum_required() to configure
# under CMake >= 4.0 without extra command-line flags.
if(POLICY CMP0167)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Minimum policy version for fetched dependencies" FORCE)
endif()

project(extropian-physics VERSION 0.1.0 LANGUAGES CXX
    DESCRIPTION "Extropian Physics engine: core runtime, mesh, numerics, physics modules, coupling, presets (implementation_spec v0.3 conformance)"
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

# ── The engine (spec §57 layout: engine/{core,mesh,discretization,numerics,physics,coupling,fidelity,presets,backends,output,diagnostics}) ──
add_subdirectory(engine)

# ── Tests ─────────────────────────────────────────────
option(EXT_PHYSICS_BUILD_TESTS "Build tests" OFF)

if(EXT_PHYSICS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(engine/tests)
endif()
