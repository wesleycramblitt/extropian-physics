#pragma once

// ─────────────────────────────────────────────────────
// Materials (implementation_spec §33).
//
// Materials are independent of discretization: they are
// constitutive data (density, viscosity, conductivity,
// specific heat, modulus, permittivity …) with optional
// field/temperature dependence (v1: constant properties;
// tabulated/field-dependent = roadmap).  Physics modules
// consume materials; the coupling layer validates the
// material ↔ field association.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/core/units.hpp>

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace exd::engine::core {

enum class MaterialType : uint8_t { Fluid, Solid, Electromagnetic };

/// Constant constitutive properties (SI).  Field/temperature-dependent
/// materials are roadmap (spec §33 "future").
struct MaterialProperties
{
    std::string name;
    MaterialType type = MaterialType::Fluid;

    // mechanical / fluid
    double density = 1.0;                 // kg/m³
    double dynamic_viscosity = 0.001;     // Pa·s
    double youngs_modulus = 0.0;          // Pa
    double poisson_ratio = 0.3;           // —

    // thermal
    double thermal_conductivity = 0.0;    // W/(m·K)
    double specific_heat = 0.0;           // J/(kg·K)
    double thermal_expansion = 0.0;       // 1/K

    // electromagnetic
    double electrical_conductivity = 0.0; // S/m
    double relative_permittivity = 1.0;   // —
    double relative_permeability = 1.0;   // —

    [[nodiscard]] double kinematic_viscosity() const { return density > 0.0 ? dynamic_viscosity / density : 0.0; }
    [[nodiscard]] double thermal_diffusivity() const
    {
        return (density > 0.0 && specific_heat > 0.0) ? thermal_conductivity / (density * specific_heat) : 0.0;
    }
    [[nodiscard]] bool validate(ModelStatus& status) const;
};

inline bool MaterialProperties::validate(ModelStatus& status) const
{
    if (name.empty())
    {
        status.ok = false;
        status.error = "material: name must not be empty";
        return false;
    }
    if (!(density > 0.0))
    {
        status.ok = false;
        status.error = "material '" + name + "': density must be > 0";
        return false;
    }
    return true;
}

/// Named material registry (§33: materials independent of discretization).
class MaterialDatabase
{
public:
    bool add(MaterialProperties material, ModelStatus& status);
    const MaterialProperties* find(std::string_view name) const;
    size_t size() const { return materials_.size(); }

private:
    std::unordered_map<std::string, MaterialProperties> materials_;
};

inline bool MaterialDatabase::add(MaterialProperties material, ModelStatus& status)
{
    if (!material.validate(status)) return false;
    if (materials_.count(material.name))
    {
        status.ok = false;
        status.error = "material database: duplicate '" + material.name + "'";
        return false;
    }
    materials_.emplace(material.name, std::move(material));
    return true;
}

inline const MaterialProperties* MaterialDatabase::find(std::string_view name) const
{
    auto it = materials_.find(std::string(name));
    return it == materials_.end() ? nullptr : &it->second;
}

} // namespace exd::engine::core
