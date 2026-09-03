#pragma once

// ─────────────────────────────────────────────────────
// Simulation state (implementation_spec §5).
//
// State owns all evolving simulation data: named Fields,
// EntitySets, particle/rigid-body sets, auxiliary and
// solver state.  It supports CPU residency today, GPU
// residency metadata (spec §39: resident state across
// timesteps), field views, versioning/change tracking,
// checkpointing (serialization entry), and external
// access (spec §44).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/entity_set.hpp>
#include <exd/engine/core/field.hpp>
#include <exd/engine/core/memory.hpp>
#include <exd/engine/core/model_status.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace exd::engine::core {

class State
{
public:
    State() = default;
    explicit State(std::string name) : name_(std::move(name)) {}

    const std::string& name() const { return name_; }

    // ── fields ──
    Field& add_field(FieldMetadata meta, size_t size, ModelStatus& status);
    Field* field(std::string_view name);
    const Field* field(std::string_view name) const;
    size_t field_count() const { return fields_.size(); }
    const std::vector<Field>& fields() const { return fields_; }

    // ── entity sets ──
    EntitySet& add_entity_set(EntitySet set, ModelStatus& status);
    EntitySet* entity_set(std::string_view name);
    const EntitySet* entity_set(std::string_view name) const;
    size_t entity_set_count() const { return entity_sets_.size(); }

    // ── residency / sync metadata (§5: CPU/GPU residency, synchronization) ──
    MemorySpace space() const { return space_; }
    void set_space(MemorySpace s) { space_ = s; }
    void sync_from_cpu();    // host → device (no-op until GPU backend)
    void sync_to_cpu();      // device → host (no-op until GPU backend)

    // ── versioning / change tracking ──
    uint64_t version() const { return version_; }
    void touch() { ++version_; }

    // ── validation ──
    bool validate(ModelStatus& status) const;

private:
    std::string name_;
    std::vector<Field> fields_;
    std::vector<EntitySet> entity_sets_;
    MemorySpace space_ = MemorySpace::Cpu;
    uint64_t version_ = 0;
};

inline Field& State::add_field(FieldMetadata meta, size_t size, ModelStatus& status)
{
    for (auto& f : fields_)
    {
        if (f.metadata().name == meta.name)
        {
            status.ok = false;
            status.error = "state '" + name_ + "': duplicate field '" + meta.name + "'";
            return f;  // return existing; caller must check status
        }
    }
    fields_.emplace_back(std::move(meta), size);
    ++version_;
    return fields_.back();
}

inline Field* State::field(std::string_view name)
{
    for (auto& f : fields_)
        if (f.metadata().name == name)
            return &f;
    return nullptr;
}

inline const Field* State::field(std::string_view name) const
{
    for (auto& f : fields_)
        if (f.metadata().name == name)
            return &f;
    return nullptr;
}

inline EntitySet& State::add_entity_set(EntitySet set, ModelStatus& status)
{
    if (!set.validate(status))
        return entity_sets_.emplace_back();  // invalid; caller checks status
    for (auto& s : entity_sets_)
    {
        if (s.name() == set.name())
        {
            status.ok = false;
            status.error = "state '" + name_ + "': duplicate entity set '" + set.name() + "'";
            return s;
        }
    }
    entity_sets_.push_back(std::move(set));
    ++version_;
    return entity_sets_.back();
}

inline EntitySet* State::entity_set(std::string_view name)
{
    for (auto& s : entity_sets_)
        if (s.name() == name)
            return &s;
    return nullptr;
}

inline const EntitySet* State::entity_set(std::string_view name) const
{
    for (auto& s : entity_sets_)
        if (s.name() == name)
            return &s;
    return nullptr;
}

inline void State::sync_from_cpu()
{
    // GPU backend attaches device storage here (spec Phase 11). CPU: no-op.
}

inline void State::sync_to_cpu()
{
    // GPU backend copies device → host here. CPU: no-op.
}

inline bool State::validate(ModelStatus& status) const
{
    for (auto& f : fields_)
        if (!f.validate(status))
            return false;
    for (auto& s : entity_sets_)
        if (!s.validate(status))
            return false;
    return true;
}

} // namespace exd::engine::core
