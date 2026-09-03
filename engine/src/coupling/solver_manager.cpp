// Solver manager — discovery, lifecycle, and orchestration of ISolverPlugin instances.

#include <exd/engine/coupling/plugin_interface.hpp>
#include <exd/engine/coupling/solver_manager.hpp>

#include <vector>
#include <memory>
#include <string>
#include <cstdio>
#include <unordered_map>

namespace exd::engine::coupling {

struct SolverManager::Impl {
    std::unordered_map<std::string, std::unique_ptr<ISolverPlugin>> instances;
    std::vector<CreateSolverFn> registered_factories;
};

SolverManager::SolverManager() : impl_(std::make_unique<Impl>()) {}
SolverManager::~SolverManager() = default;

SolverManager& SolverManager::instance() {
    static SolverManager mgr;
    return mgr;
}

void SolverManager::register_factory(const std::string& /*name*/, CreateSolverFn factory) {
    impl_->registered_factories.push_back(factory);
}

ISolverPlugin* SolverManager::create(const std::string& name) {
    auto it = impl_->instances.find(name);
    if (it != impl_->instances.end()) return it->second.get();

    // Try registered factories
    for (auto& factory : impl_->registered_factories) {
        auto* plugin = factory();
        if (plugin && plugin->name() == name) {
            impl_->instances[name] = std::unique_ptr<ISolverPlugin>(plugin);
            std::printf("[SolverManager] Created solver: %s\n", name.c_str());
            return plugin;
        }
        delete plugin; // wrong plugin, discard
    }
    return nullptr;
}

void SolverManager::destroy(const std::string& name) {
    impl_->instances.erase(name);
}

void SolverManager::destroy_all() {
    impl_->instances.clear();
}

std::vector<std::string> SolverManager::available_solvers() const {
    std::vector<std::string> names;
    for (auto& [name, _] : impl_->instances) names.push_back(name);
    return names;
}

} // namespace exd::engine::coupling
