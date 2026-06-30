#pragma once
#include <string>
#include <vector>
#include <memory>
namespace exd::physics {
class ISolverPlugin;
using CreateSolverFn = ISolverPlugin* (*)();
class SolverManager {
public:
    SolverManager(); ~SolverManager();
    static SolverManager& instance();
    void register_factory(const std::string& name, CreateSolverFn);
    ISolverPlugin* create(const std::string& name);
    void destroy(const std::string& name);
    void destroy_all();
    std::vector<std::string> available_solvers() const;
private:
    struct Impl; std::unique_ptr<Impl> impl_;
};
}
