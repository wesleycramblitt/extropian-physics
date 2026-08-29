#pragma once

#include "field_sampler.hpp"   // IFlowField3D (for the fluid adapters)

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace exd::physics::coupling {

// ─────────────────────────────────────────────────────
// Generic field channels: the data currency for coupling
// between ANY domains (fluid, EM, thermal, structural…).
//
//   IVectorField3D — E, B, displacement, velocity…
//   IScalarField3D — temperature, potential, pressure…
//
// IFlowField3D (fluid) carries vector + scalar together
// and adapters expose its parts as generic channels.
// ─────────────────────────────────────────────────────

class IVectorField3D
{
public:
    virtual ~IVectorField3D() = default;
    /// Sample the vector field at p. Returns false when out of bounds.
    virtual bool sample(const std::array<double, 3>& p,
                        std::array<double, 3>& value_out) const = 0;
};

class IScalarField3D
{
public:
    virtual ~IScalarField3D() = default;
    /// Sample the scalar field at p. Returns false when out of bounds.
    virtual bool sample(const std::array<double, 3>& p,
                        double& value_out) const = 0;
};

// ── Structured regular grids (node-centered, trilinear sampling) ──

struct StructuredScalarGrid
{
    std::array<double, 3> origin = {0, 0, 0};   // node (0,0,0) (m)
    std::array<double, 3> spacing = {0.1, 0.1, 0.1}; // node spacing (m)
    std::array<int32_t, 3> dims = {0, 0, 0};    // node counts (>= 2 per axis)
    std::vector<double> values;                 // nx·ny·nz flat, index i + nx·(j + ny·k)
};

struct StructuredVectorGrid
{
    std::array<double, 3> origin = {0, 0, 0};
    std::array<double, 3> spacing = {0.1, 0.1, 0.1};
    std::array<int32_t, 3> dims = {0, 0, 0};
    std::vector<double> values;                 // 3·nx·ny·nz flat
};

/// Trilinear scalar channel over a structured grid. Returns nullptr on
/// invalid config (dims < 2, spacing <= 0, size mismatch).
std::unique_ptr<IScalarField3D> make_scalar_grid_field(const StructuredScalarGrid& grid);

/// Trilinear vector channel over a structured grid. Returns nullptr on
/// invalid config.
std::unique_ptr<IVectorField3D> make_vector_grid_field(const StructuredVectorGrid& grid);

// ── Fluid field adapters (IFlowField3D → generic channels) ──
// The caller must keep the underlying field alive.

/// Velocity component of a fluid field as a vector channel.
std::unique_ptr<IVectorField3D> make_velocity_field_adapter(const IFlowField3D& field);

/// Pressure component of a fluid field as a scalar channel.
std::unique_ptr<IScalarField3D> make_pressure_field_adapter(const IFlowField3D& field);

} // namespace exd::physics::coupling
