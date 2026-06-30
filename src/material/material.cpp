// Material property database

#include <string>
#include <unordered_map>
#include <vector>

namespace exd::physics::material {

enum class MaterialType { Fluid, Solid, Electromagnetic };

struct MaterialProperties {
    std::string name;
    MaterialType type = MaterialType::Fluid;

    // Fluid properties
    double density = 1.0;            // kg/m^3
    double dynamic_viscosity = 0.001;// Pa·s
    double kinematic_viscosity = 0;  // m^2/s (computed from density+viscosity)

    // Solid properties
    double youngs_modulus = 0;       // Pa
    double poisson_ratio = 0.3;
    double yield_strength = 0;       // Pa
    double thermal_expansion = 0;    // 1/K
    double thermal_conductivity = 0; // W/(m·K)
    double specific_heat = 0;        // J/(kg·K)

    // Electromagnetic
    double electrical_conductivity = 0; // S/m
    double relative_permittivity = 1.0;
    double relative_permeability = 1.0;

    void compute_derived() {
        if (density > 0)
            kinematic_viscosity = dynamic_viscosity / density;
    }
};

/// Built-in material database
class MaterialDatabase {
public:
    MaterialDatabase() { load_defaults(); }

    const MaterialProperties* find(const std::string& name) const {
        auto it = materials_.find(name);
        return (it != materials_.end()) ? &it->second : nullptr;
    }

    void add(MaterialProperties mat) {
        mat.compute_derived();
        materials_[mat.name] = std::move(mat);
    }

    const auto& all() const { return materials_; }

private:
    std::unordered_map<std::string, MaterialProperties> materials_;

    void load_defaults() {
        add({"Air", MaterialType::Fluid, 1.225, 1.81e-5});
        add({"Water", MaterialType::Fluid, 1000.0, 0.001});
        add({"Aluminum 6061", MaterialType::Solid, 2700.0,
             .youngs_modulus=68.9e9, .poisson_ratio=0.33, .yield_strength=276e6,
             .thermal_expansion=23.1e-6, .thermal_conductivity=167, .specific_heat=896});
        add({"Steel A36", MaterialType::Solid, 7850.0,
             .youngs_modulus=200e9, .poisson_ratio=0.26, .yield_strength=250e6,
             .thermal_expansion=11.7e-6, .thermal_conductivity=50, .specific_heat=480});
        add({"Copper", MaterialType::Electromagnetic, 8960.0,
             .electrical_conductivity=5.96e7, .relative_permittivity=1.0,
             .relative_permeability=0.999994, .thermal_conductivity=401});
    }
};

} // namespace exd::physics::material
