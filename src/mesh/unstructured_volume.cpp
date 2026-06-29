#include <vector>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace exd::physics::mesh {

using NodeId = int64_t;
using CellId = int64_t;
using FaceId = int64_t;

/// Cell types
enum class CellType : uint8_t {
    Tet4 = 4, Hex8 = 8, Prism6 = 6, Pyramid5 = 5,
    Tri3 = 3, Quad4 = 4
};

/// 3D unstructured volume mesh
struct UnstructuredVolume {
    // Geometry
    std::vector<double> node_coords;      // x,y,z interleaved, size = N*3
    std::vector<int64_t> cell_nodes;      // packed: cell_type, n0, n1, ...
    std::vector<int64_t> cell_offsets;    // start index per cell in cell_nodes
    std::vector<CellType> cell_types;

    // Named boundary surfaces
    std::unordered_map<std::string, std::vector<FaceId>> boundary_faces;

    // Metrics
    size_t node_count() const { return node_coords.size() / 3; }
    size_t cell_count() const { return cell_types.size(); }

    void add_cell(CellType type, std::initializer_list<int64_t> nodes) {
        cell_offsets.push_back(static_cast<int64_t>(cell_nodes.size()));
        cell_nodes.push_back(static_cast<int64_t>(type));
        for (auto n : nodes) cell_nodes.push_back(n);
        cell_types.push_back(type);
    }
};

/// Triangle surface mesh
struct SurfaceMesh {
    std::vector<double> node_coords;       // x,y,z interleaved
    std::vector<int64_t> triangles;        // 3 indices per tri
    size_t tri_count() const { return triangles.size() / 3; }
};

/// Structured (uniform) grid
struct StructuredGrid {
    int nx, ny, nz;
    double dx, dy, dz;
    double origin_x, origin_y, origin_z;

    size_t cell_count() const { return (size_t)nx * ny * nz; }
    size_t node_count() const { return (size_t)(nx+1)*(ny+1)*(nz+1); }

    void node_ijk_to_xyz(int i, int j, int k, double& x, double& y, double& z) const {
        x = origin_x + i * dx;
        y = origin_y + j * dy;
        z = origin_z + k * dz;
    }
};

} // namespace exd::physics::mesh
