#pragma once

// ─────────────────────────────────────────────────────
// Compatibility rules (implementation_spec §21–§22, §55).
//
// A general validation system: CompatibilityRule objects
// inspect Module / Discretization / Mesh / Field /
// Material / BoundaryCondition / Coupling / Solver /
// FidelityProfile / Backend and reject invalid
// configurations with explicit diagnostics BEFORE
// execution.  Rules are composable and extensible.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/coupling/contracts.hpp>

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace exd::engine::coupling {

/// What a rule may inspect (spec §22).  v1 keeps the context light: the
/// contract plus free-text module/mesh/backend declarations.
struct ValidationContext
{
    const CouplingContract* contract = nullptr;      // when validating a coupling
    std::string module_name;                         // owning module
    std::string discretization;                      // e.g. "FDM"
    std::string mesh_family;                         // e.g. "structured"
    std::string backend;                             // e.g. "cpu"
};

class CompatibilityRule
{
public:
    virtual ~CompatibilityRule() = default;
    virtual std::string_view name() const = 0;
    /// Returns true when the configuration is compatible; on failure sets
    /// status.ok = false with an explicit diagnostic (spec §22 example
    /// style: what is invalid, why, which modules, alternatives).
    virtual bool check(const ValidationContext& ctx, core::ModelStatus& status) const = 0;
};

/// Built-in rule: a coupling contract must be internally consistent and
/// dimensionally closed (spec §21: quantity type, units, rank).
class ContractRule final : public CompatibilityRule
{
public:
    std::string_view name() const override { return "coupling_contract"; }
    bool check(const ValidationContext& ctx, core::ModelStatus& status) const override
    {
        if (!ctx.contract) return true;
        return ctx.contract->validate(status);
    }
};

/// Built-in rule: module declares the discretization it actually runs
/// (spec §12–§13: per-module discretization support).
class DiscretizationRule final : public CompatibilityRule
{
public:
    explicit DiscretizationRule(std::vector<std::string> supported)
        : supported_(std::move(supported)) {}
    std::string_view name() const override { return "module_discretization"; }
    bool check(const ValidationContext& ctx, core::ModelStatus& status) const override
    {
        if (ctx.discretization.empty()) return true;
        for (auto& s : supported_)
            if (s == ctx.discretization)
                return true;
        status.ok = false;
        status.error = "module '" + ctx.module_name + "': discretization '" +
                       ctx.discretization + "' is not supported";
        for (auto& s : supported_)
            status.warnings.push_back("supported: " + s);
        return false;
    }

private:
    std::vector<std::string> supported_;
};

/// Built-in rule: mesh family vs discretization (FDM requires a structured
/// lattice; unstructured is roadmap — spec §9).
class MeshDiscretizationRule final : public CompatibilityRule
{
public:
    std::string_view name() const override { return "mesh_discretization"; }
    bool check(const ValidationContext& ctx, core::ModelStatus& status) const override
    {
        if (ctx.mesh_family.empty() || ctx.discretization.empty()) return true;
        if (ctx.discretization == "FDM" && ctx.mesh_family == "unstructured")
        {
            status.ok = false;
            status.error = "FDM requires a structured mesh (unstructured is roadmap)";
            return false;
        }
        return true;
    }
};

/// Rule registry: run all registered rules over a context (spec §22:
/// Configuration → Validation → Compatibility Graph → Execution Graph).
class RuleRegistry
{
public:
    void add(std::unique_ptr<CompatibilityRule> rule) { rules_.push_back(std::move(rule)); }
    bool validate(const ValidationContext& ctx, core::ModelStatus& status) const
    {
        for (auto& rule : rules_)
        {
            if (!rule->check(ctx, status))
                return false;
        }
        return true;
    }
    size_t rule_count() const { return rules_.size(); }

private:
    std::vector<std::unique_ptr<CompatibilityRule>> rules_;
};

/// Registry preloaded with the built-in rules.
inline RuleRegistry default_rule_registry()
{
    RuleRegistry reg;
    reg.add(std::make_unique<ContractRule>());
    reg.add(std::make_unique<MeshDiscretizationRule>());
    return reg;
}

} // namespace exd::engine::coupling
