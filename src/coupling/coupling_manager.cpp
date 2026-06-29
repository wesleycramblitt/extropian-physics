// Coupling manager — orchestrates data exchange between solvers.

#include <string>
#include <unordered_map>
#include <memory>
#include <cstdio>

namespace exd::physics {

/// Surface mapping: transfer field data between non-matching surface meshes.
struct SurfaceMapper {
    std::string source_surface;
    std::string target_surface;
    std::string field_name;

    // Transfer using nearest-neighbor interpolation (TODO: radial basis, MLS)
    void transfer(const double* source_data, double* target_data,
                  size_t source_count, size_t target_count) const {
        // Simplified: copy first value (placeholder for real interpolation)
        if (target_count > 0) target_data[0] = source_data[0];
        (void)source_count;
    }
};

/// Coupling manager: handles two-way coupled simulations.
class CouplingManager {
public:
    struct CoupledPair {
        std::string solver_a;
        std::string solver_b;
        std::string coupling_type; // "fsi", "cht", "general"
        std::vector<SurfaceMapper> mappers;
        double coupling_interval = 0.01; // seconds between data exchanges
        double last_exchange = 0.0;
    };

    void add_coupling(const std::string& solver_a, const std::string& solver_b,
                      const std::string& type) {
        pairs_.push_back({solver_a, solver_b, type});
        std::printf("[Coupling] Added %s ↔ %s (%s)\n",
                    solver_a.c_str(), solver_b.c_str(), type.c_str());
    }

    bool should_exchange(const CoupledPair& pair, double current_time) const {
        return (current_time - pair.last_exchange) >= pair.coupling_interval;
    }

    void exchange(const std::string& solver_name, double current_time) {
        for (auto& pair : pairs_) {
            if ((pair.solver_a == solver_name || pair.solver_b == solver_name)
                && should_exchange(pair, current_time)) {
                for (auto& mapper : pair.mappers)
                    mapper.transfer(nullptr, nullptr, 0, 0); // placeholder
                pair.last_exchange = current_time;
            }
        }
    }

    const auto& pairs() const { return pairs_; }

private:
    std::vector<CoupledPair> pairs_;
};

} // namespace exd::physics
