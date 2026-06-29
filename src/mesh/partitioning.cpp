#include <vector>
#include <cstdint>
namespace exd::physics::mesh {
void partition_uniform(int num_parts, int num_cells, std::vector<int>& part_ids) {
    part_ids.resize(num_cells);
    for (int i = 0; i < num_cells; ++i) part_ids[i] = i % num_parts;
}
}
