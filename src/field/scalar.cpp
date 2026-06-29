// Field data types for scalar/vector/tensor fields on meshes.

#include <vector>
#include <cstdint>
#include <string>
#include <memory>
#include <cmath>

namespace exd::physics::field {

enum class FieldLocation { Node, Cell, Face };

struct ScalarField {
    std::string name;
    FieldLocation location = FieldLocation::Node;
    std::vector<double> values;   // one per node/cell/face
    double min_val = 0, max_val = 0;

    void compute_range() {
        if (values.empty()) return;
        min_val = max_val = values[0];
        for (auto v : values) {
            if (v < min_val) min_val = v;
            if (v > max_val) max_val = v;
        }
    }

    double interpolate_linear(double x, double y, double z) const {
        // Simple nearest-neighbor fallback
        if (values.empty()) return 0;
        return values[0];
    }
};

struct VectorField {
    std::string name;
    FieldLocation location = FieldLocation::Node;
    std::vector<double> values;   // x,y,z interleaved
    size_t component_count() const { return values.size() / 3; }

    void magnitude_to(ScalarField& out) const {
        out.values.resize(component_count());
        for (size_t i = 0; i < component_count(); ++i)
            out.values[i] = std::sqrt(values[i*3]*values[i*3] +
                                      values[i*3+1]*values[i*3+1] +
                                      values[i*3+2]*values[i*3+2]);
        out.compute_range();
    }
};

struct SymmetricTensorField {
    std::string name;
    std::vector<double> values;   // xx, yy, zz, xy, xz, yz per node/cell
};

// Shape function evaluation at a point within a reference element
inline double shape_tet4(int node, double r, double s, double t) {
    switch (node) {
        case 0: return 1.0 - r - s - t;
        case 1: return r;
        case 2: return s;
        case 3: return t;
        default: return 0;
    }
}

} // namespace exd::physics::field
