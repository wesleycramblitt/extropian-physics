#pragma once

#include <exd/engine/physics/fluid/forces/force_evaluator.hpp>
#include <exd/engine/physics/fluid/forces/flow_types.hpp>
#include <exd/engine/physics/fluid/fdm/fdm_result.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_solver.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <span>

namespace exd::engine::coupling {

// ─────────────────────────────────────────────────────
// Flow-field sampling: CFD field → flow at surface points.
// The coupling leg "CFD → forces" of a coupled simulation.
// ─────────────────────────────────────────────────────

/// A 3D flow field evaluable at arbitrary points.
class IFlowField3D
{
public:
    virtual ~IFlowField3D() = default;

    /// Sample velocity and pressure at `p`. Returns false when out of bounds.
    virtual bool sample(const std::array<double, 3>& p,
                        std::array<double, 3>& velocity_out,
                        double& pressure_out) const = 0;

    virtual double density() const = 0;
    virtual double viscosity() const = 0;
};

/// Uniform freestream field (used by standalone simulations).
struct UniformFieldConfig
{
    std::array<double, 3> velocity = {0.0, 0.0, 0.0}; // m/s
    double rho = 1.225;                               // kg/m³
    double mu = 1.81e-5;                              // Pa·s
    double p_ref = 101325.0;                          // Pa
};
std::unique_ptr<IFlowField3D> make_uniform_field(const UniformFieldConfig& config);

/// Structured regular grid field with trilinear interpolation.
/// Velocity stored flat: 3·nx·ny·nz, index(i,j,k) = i + nx·(j + ny·k).
struct StructuredGridConfig
{
    std::array<double, 3> origin = {0.0, 0.0, 0.0};  // node (0,0,0) (m)
    std::array<double, 3> spacing = {0.1, 0.1, 0.1}; // node spacing per axis (m)
    std::array<int32_t, 3> dims = {0, 0, 0};         // node counts (≥ 2 per axis)
    double rho = 1.225;
    double mu = 1.81e-5;
    std::vector<double> velocity;  // 3·nx·ny·nz
    std::vector<double> pressure;  // nx·ny·nz
};
std::unique_ptr<IFlowField3D> make_structured_grid_field(const StructuredGridConfig& config);

/// Adapter that samples a 2D FDM field at z ≈ 0 (bilinear in-plane,
/// vz = 0). Lets the current 2D FDM CFD solver feed the 3D pipeline.
std::unique_ptr<IFlowField3D> make_fdm_field_adapter(const exd::engine::physics::fluid::fdm::FDMFieldData& field,
                                                     double rho, double mu, double p_ref);

/// Adapter that samples a 3D FDM (fdm3) field using trilinear
/// interpolation with the cell-center convention
/// (x[i] = (i + 0.5) * dx, etc.).  Out-of-bounds queries return false.
std::unique_ptr<IFlowField3D> make_fdm3_field_adapter(const exd::engine::physics::fluid::fdm3::FDM3Solver& solver);

/// Sample a flow field onto surface arrays → SurfaceFlow for evaluators.
/// Shear traction is left zero (no wall model here).
exd::engine::physics::fluid::forces::SurfaceFlow sample_flow(const IFlowField3D& field,
                                       std::span<const std::array<double, 3>> points,
                                       std::span<const std::array<double, 3>> normals,
                                       std::span<const double> areas,
                                       std::span<const int32_t> element_index,
                                       double p_ref);

/// Convenience: sample onto a prebuilt blade surface.
exd::engine::physics::fluid::forces::SurfaceFlow sample_flow(const IFlowField3D& field,
                                       const exd::engine::physics::fluid::forces::BladeSurface& surface,
                                       double p_ref);

} // namespace exd::engine::coupling
