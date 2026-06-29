#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <span>
#include <vector>
#include <unordered_map>
#include <memory>

namespace ext::physics {

// ─────────────────────────────────────────────────────
// Domain classification
// ─────────────────────────────────────────────────────

enum class PhysicsDomain : uint8_t {
    FluidFlow,
    SolidMechanics,
    HeatTransfer,
    Electromagnetics,
    Acoustics,
    MultibodyDynamics,
    ParticleDynamics,
    Custom,
};

// ─────────────────────────────────────────────────────
// Typed field accessor — solver exposes fields through this
// ─────────────────────────────────────────────────────

enum class FieldDataType : uint8_t {
    ScalarFloat,
    ScalarDouble,
    VectorFloat,
    VectorDouble,
    TensorFloat,
    TensorDouble,
};

struct FieldAccessor {
    std::string name;
    FieldDataType type;
    uint64_t num_components; ///< 1=scalar, 3=vector, 6/9=symmetric/full tensor
    uint64_t num_elements;   ///< Number of nodes or cells

    /// Copy field data into a pre-allocated buffer. Returns bytes written.
    virtual uint64_t copy_to(void* buffer, uint64_t buffer_size_bytes) const = 0;
    virtual ~FieldAccessor() = default;
};

// ─────────────────────────────────────────────────────
// Coupling surface — data exchange between solvers
// ─────────────────────────────────────────────────────

struct CouplingSurface {
    std::string name;
    uint64_t num_vertices;
    uint64_t num_faces;

    /// Vertex coordinates (3 * num_vertices doubles)
    virtual void get_vertex_coords(double* out) const = 0;
    /// Face connectivity (face_count, then vertex indices per face)
    virtual void get_face_connectivity(int64_t* out, uint64_t max_size) const = 0;
    /// Field values on this surface
    virtual std::unique_ptr<FieldAccessor> get_field(const std::string& name) const = 0;
    /// Write field values to this surface (from another solver)
    virtual void set_field(const std::string& name, const void* data, uint64_t size_bytes) = 0;

    virtual ~CouplingSurface() = default;
};

// ─────────────────────────────────────────────────────
// THE SOLVER PLUGIN INTERFACE
// ─────────────────────────────────────────────────────

struct SolverMesh {
    /// Node coordinates (Nx3, double)
    std::span<const double> node_coords;
    /// Element connectivity (packed: cell_type, then node indices per cell)
    std::span<const int64_t> element_nodes;
    /// Element types per cell (4=tet, 5=hex, 6=prism, 7=pyramid, etc.)
    std::span<const uint8_t> element_types;
    /// Named boundary surface groups
    std::unordered_map<std::string, std::span<const int64_t>> boundary_faces;
};

struct BoundaryCondition {
    std::string region_name;
    PhysicsDomain domain;
    std::string type;       ///< "fixed_value", "zero_gradient", "slip", "inlet", "outlet", etc.
    std::unordered_map<std::string, double> parameters;
};

struct MaterialAssignment {
    std::string region_name;
    std::string material_name;
};

/// Configuration blob — solver-specific, parsed by the plugin.
using ConfigNode = std::unordered_map<std::string, std::string>;

/// @brief Every solver plugin MUST implement this interface.
///
/// Lifecycle: initialize() → step()* → finalize()
/// Field and coupling queries are valid between initialize() and finalize().
class ISolverPlugin {
public:
    virtual ~ISolverPlugin() = default;

    /// ── Identity ──────────────────────────────────────
    [[nodiscard]] virtual std::string_view name()    const = 0;
    [[nodiscard]] virtual std::string_view version() const = 0;
    [[nodiscard]] virtual PhysicsDomain domain()     const = 0;

    /// ── Capability query ──────────────────────────────
    [[nodiscard]] virtual std::vector<std::string_view> supported_bc_types() const = 0;

    /// ── Lifecycle ─────────────────────────────────────
    /// Called once before any step(). Mesh, BCs, materials, and solver params
    /// are stable for the lifetime of the solve.
    virtual void initialize(
        const SolverMesh& mesh,
        std::span<const BoundaryCondition> bcs,
        std::span<const MaterialAssignment> materials,
        const ConfigNode& solver_params
    ) = 0;

    /// Advance one timestep. Returns false when steady-state or converged.
    virtual bool step(double dt) = 0;

    /// Release all solver resources.
    virtual void finalize() = 0;

    /// ── Data extraction ───────────────────────────────
    /// Get a field by name (e.g., "velocity", "pressure", "stress").
    /// Valid after initialize(), invalid after finalize().
    [[nodiscard]] virtual std::unique_ptr<FieldAccessor> get_field(
        const std::string& field_name
    ) = 0;

    /// ── Coupling ──────────────────────────────────────
    /// Get the coupling surface by name for data exchange.
    [[nodiscard]] virtual std::unique_ptr<CouplingSurface> get_coupling_surface(
        const std::string& surface_name
    ) = 0;
};

/// Factory function type: each solver plugin shared library exports this symbol.
using CreateSolverFn = ISolverPlugin* (*)();
using DestroySolverFn = void (*)(ISolverPlugin*);

/// Exported symbol names expected by the plugin loader.
constexpr const char* SOLVER_CREATE_SYMBOL  = "ext_solver_create";
constexpr const char* SOLVER_DESTROY_SYMBOL = "ext_solver_destroy";

} // namespace ext::physics
