// Boundary condition framework

#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace exd::physics::bc {

enum class BCType {
    FixedValue,       // Dirichlet: u = specified
    ZeroGradient,     // Neumann: du/dn = 0
    Slip,             // No normal velocity
    NoSlip,           // u = 0
    InletVelocity,    // u = specified velocity vector
    OutletPressure,   // p = specified pressure
    Periodic,         // Values match on paired surfaces
    Symmetry,         // Normal gradient = 0, tangential = free
    FixedDisplacement,// u = specified displacement (solid mechanics)
    FixedTraction,    // t = specified traction vector
};

struct BoundaryCondition {
    std::string name;               // user-facing name
    std::string region_name;        // mesh boundary patch name
    BCType type = BCType::FixedValue;

    // Parameter values (depends on type)
    std::unordered_map<std::string, double> scalar_params;
    std::vector<double> vector_param; // for velocity, traction, etc.

    static BoundaryCondition fixed_value(const std::string& region, double val) {
        BoundaryCondition bc{region + "_fixed", region, BCType::FixedValue};
        bc.scalar_params["value"] = val;
        return bc;
    }

    static BoundaryCondition inlet_velocity(const std::string& region,
                                             double vx, double vy, double vz) {
        BoundaryCondition bc{"inlet", region, BCType::InletVelocity};
        bc.vector_param = {vx, vy, vz};
        return bc;
    }

    static BoundaryCondition outlet_pressure(const std::string& region, double p) {
        BoundaryCondition bc{"outlet", region, BCType::OutletPressure};
        bc.scalar_params["pressure"] = p;
        return bc;
    }
};

/// Serialize/deserialize BCs to a string for project files
inline std::string bc_type_to_string(BCType t) {
    switch (t) {
        case BCType::FixedValue:       return "fixed_value";
        case BCType::ZeroGradient:     return "zero_gradient";
        case BCType::Slip:             return "slip";
        case BCType::NoSlip:           return "no_slip";
        case BCType::InletVelocity:    return "inlet_velocity";
        case BCType::OutletPressure:   return "outlet_pressure";
        case BCType::Periodic:         return "periodic";
        case BCType::Symmetry:         return "symmetry";
        case BCType::FixedDisplacement:return "fixed_displacement";
        case BCType::FixedTraction:     return "fixed_traction";
    }
    return "unknown";
}

} // namespace exd::physics::bc
