#include <string>
#include <unordered_map>

namespace exd::physics::material {

enum class MaterialType { Fluid, Solid, Electromagnetic };

struct MaterialProperties {
    std::string name;
    MaterialType type = MaterialType::Fluid;
    double density = 1.0;
    double dynamic_viscosity = 0.001;
    double kinematic_viscosity = 0;
    double youngs_modulus = 0;
    double poisson_ratio = 0.3;
    double yield_strength = 0;
    double thermal_expansion = 0;
    double thermal_conductivity = 0;
    double specific_heat = 0;
    double electrical_conductivity = 0;
    double relative_permittivity = 1.0;
    double relative_permeability = 1.0;

    void compute_derived() {
        if (density > 0) kinematic_viscosity = dynamic_viscosity / density;
    }
};

class MaterialDatabase {
    std::unordered_map<std::string, MaterialProperties> materials_;

    void add(MaterialProperties m) { m.compute_derived(); materials_[m.name] = m; }

public:
    MaterialDatabase() {
        MaterialProperties air;
        air.name = "Air"; air.type = MaterialType::Fluid;
        air.density = 1.225; air.dynamic_viscosity = 1.81e-5;
        add(air);

        MaterialProperties water;
        water.name = "Water"; water.type = MaterialType::Fluid;
        water.density = 1000.0; water.dynamic_viscosity = 0.001;
        add(water);

        MaterialProperties al;
        al.name = "Aluminum 6061"; al.type = MaterialType::Solid;
        al.density = 2700; al.youngs_modulus = 68.9e9; al.poisson_ratio = 0.33;
        al.yield_strength = 276e6; al.thermal_expansion = 23.1e-6;
        al.thermal_conductivity = 167; al.specific_heat = 896;
        add(al);

        MaterialProperties steel;
        steel.name = "Steel A36"; steel.type = MaterialType::Solid;
        steel.density = 7850; steel.youngs_modulus = 200e9; steel.poisson_ratio = 0.26;
        steel.yield_strength = 250e6; steel.thermal_expansion = 11.7e-6;
        steel.thermal_conductivity = 50; steel.specific_heat = 480;
        add(steel);

        MaterialProperties cu;
        cu.name = "Copper"; cu.type = MaterialType::Electromagnetic;
        cu.density = 8960; cu.electrical_conductivity = 5.96e7;
        cu.relative_permittivity = 1.0; cu.relative_permeability = 0.999994;
        cu.thermal_conductivity = 401;
        add(cu);
    }

    const MaterialProperties* find(const std::string& name) const {
        auto it = materials_.find(name);
        return (it != materials_.end()) ? &it->second : nullptr;
    }
};

} // namespace exd::physics::material
