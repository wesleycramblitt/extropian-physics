#pragma once

// ─────────────────────────────────────────────────────
// Boundary conditions (implementation_spec §32).
//
// Boundary conditions are first-class objects.  Each
// carries its kind, patch, and scalar/vector parameters;
// conditions declare compatibility with physics, field,
// discretization, and mesh location (validated by the
// coupling rules layer before execution).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/mesh/structured.hpp>

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace exd::engine::mesh {

enum class BoundaryConditionKind : uint8_t
{
    Dirichlet,    // value prescribed
    Neumann,      // normal flux prescribed
    Robin,        // a·u + b·du/dn = c
    Periodic,
    Wall,         // no-slip wall (fluid)
    Slip,         // tangential slip (fluid)
    Inlet,        // inflow with prescribed velocity/temperature
    Outlet,       // outflow (often zero-gradient)
    Symmetry,
    Radiation,
    Flux,
    Convective,
};

constexpr const char* to_string(BoundaryConditionKind k)
{
    switch (k)
    {
    case BoundaryConditionKind::Dirichlet: return "dirichlet";
    case BoundaryConditionKind::Neumann: return "neumann";
    case BoundaryConditionKind::Robin: return "robin";
    case BoundaryConditionKind::Periodic: return "periodic";
    case BoundaryConditionKind::Wall: return "wall";
    case BoundaryConditionKind::Slip: return "slip";
    case BoundaryConditionKind::Inlet: return "inlet";
    case BoundaryConditionKind::Outlet: return "outlet";
    case BoundaryConditionKind::Symmetry: return "symmetry";
    case BoundaryConditionKind::Radiation: return "radiation";
    case BoundaryConditionKind::Flux: return "flux";
    case BoundaryConditionKind::Convective: return "convective";
    }
    return "?";
}

struct BoundaryCondition
{
    std::string name;                       // user-facing name
    BoundaryId patch = BoundaryId::XNeg;    // mesh patch this condition applies to
    BoundaryConditionKind kind = BoundaryConditionKind::Dirichlet;
    std::unordered_map<std::string, double> scalar_params;  // e.g. value, flux, h, T_inf
    std::vector<double> vector_param;       // e.g. inlet velocity vector
    std::string applies_to_field;           // field name ("" = all)

    bool validate(ModelStatus& status) const;
};

inline bool BoundaryCondition::validate(ModelStatus& status) const
{
    if (name.empty())
    {
        status.ok = false;
        status.error = "boundary condition: name must not be empty";
        return false;
    }
    return true;
}

} // namespace exd::engine::mesh
