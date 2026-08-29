#pragma once

#include <exd/physics/coupling/field_channels.hpp>
#include <exd/physics/model_status.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::electrical {

// ─────────────────────────────────────────────────────
// Static field solvers on a regular 3D grid (finite
// differences, SOR Poisson): electrostatics (φ → E) and
// magnetostatics (A_z → B) — the grid-friendly first win
// of the electrical domain, sharing the SOR machinery
// pattern of the FDM pressure solver.
// ─────────────────────────────────────────────────────

enum class StaticFieldMode : uint8_t
{
    Electrostatic, //  ∇²φ = 0,        E = −∇φ
    Magnetostatic, //  ∇²A_z = −μ₀·J_z, B = ∇×A  (A = A_z·ẑ)
};

struct StaticFieldConfig
{
    StaticFieldMode mode = StaticFieldMode::Electrostatic;

    // Grid (interior nodes; Dirichlet on all 6 faces + electrode boxes)
    std::array<int32_t, 3> dims = {32, 32, 32};        // >= 2 per axis
    std::array<double, 3> spacing = {0.01, 0.01, 0.01}; // node spacing (m)

    // Dirichlet values on the faces in order {x−, x+, y−, y+, z−, z+}
    // (φ in V for electrostatic; A_z in T·m for magnetostatic)
    std::array<double, 6> face_values = {0, 0, 0, 0, 0, 0};

    // Internal Dirichlet patches (electrodes / conductors)
    struct BoxPatch
    {
        std::array<double, 3> center = {0, 0, 0};       // world coords (m)
        std::array<double, 3> half_extents = {0, 0, 0}; // half box size (m)
        double value = 0.0;                             // φ or A_z
    };
    std::vector<BoxPatch> patches;

    // Magnetostatic current sources: z-directed current over a box (A/m²)
    struct CurrentBox
    {
        std::array<double, 3> center = {0, 0, 0};
        std::array<double, 3> half_extents = {0, 0, 0};
        double jz = 0.0;   // current density along +z (A/m²)
    };
    std::vector<CurrentBox> currents;

    // SOR solver settings
    double sor_omega = 1.8;         // in (0, 2)
    double tolerance = 1e-8;
    int max_iterations = 200000;
};

struct StaticFieldResult
{
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    int iterations = 0;
    double residual = 0.0;
    // φ (V) or A_z (T·m) on the interior grid
    coupling::StructuredScalarGrid potential;
    // E = −∇φ (V/m) or B = ∇×A (T), cell-centered, finite-difference
    coupling::StructuredVectorGrid field_vector;
};

/// Solve the static field. Never throws; failures reported via `ok`/`error`.
StaticFieldResult solve_static_field(const StaticFieldConfig& config);

} // namespace exd::physics::electrical
