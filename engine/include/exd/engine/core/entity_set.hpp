#pragma once

// ─────────────────────────────────────────────────────
// Entity sets (implementation_spec §7).
//
// Generalized collections of physical entities — cells,
// faces, nodes, edges, particles, rigid bodies, elements,
// DOFs — each representable through contiguous/indexed
// data structures suitable for parallel execution.
// Every entity in a set shares one index space, which is
// the precondition for later domain decomposition
// (owned/ghost regions, spec §62).
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::core {

enum class EntityKind : uint8_t
{
    Cells,
    Faces,
    Nodes,
    Edges,
    Particles,
    RigidBodies,
    Elements,
    Dofs,
};

constexpr const char* to_string(EntityKind k)
{
    switch (k)
    {
    case EntityKind::Cells: return "cells";
    case EntityKind::Faces: return "faces";
    case EntityKind::Nodes: return "nodes";
    case EntityKind::Edges: return "edges";
    case EntityKind::Particles: return "particles";
    case EntityKind::RigidBodies: return "rigid_bodies";
    case EntityKind::Elements: return "elements";
    case EntityKind::Dofs: return "dofs";
    }
    return "?";
}

/// A homogeneous, contiguous, indexed collection of entities.
/// Entity `i` owns index range [offset + i, offset + i + 1) of any field
/// associated with the set (fields store one value per entity).
class EntitySet
{
public:
    EntitySet() = default;
    EntitySet(std::string name, EntityKind kind, size_t count);

    const std::string& name() const { return name_; }
    EntityKind kind() const { return kind_; }
    size_t count() const { return count_; }
    bool empty() const { return count_ == 0; }

    /// Offset into the shared index space (0 for the first set of a domain;
    /// nonzero when multiple sets concatenate — decomposition-ready, §62).
    size_t offset() const { return offset_; }
    void set_offset(size_t o) { offset_ = o; }

    /// Optional per-entity scalar attribute rows (contiguous, count() each).
    void add_attribute(std::string attribute_name, std::vector<double> values, ModelStatus& status);
    const std::vector<double>* attribute(std::string_view attribute_name) const;

    bool validate(ModelStatus& status) const;

private:
    std::string name_;
    EntityKind kind_ = EntityKind::Cells;
    size_t count_ = 0;
    size_t offset_ = 0;
    std::vector<std::pair<std::string, std::vector<double>>> attributes_;
};

inline EntitySet::EntitySet(std::string name, EntityKind kind, size_t count)
    : name_(std::move(name)), kind_(kind), count_(count)
{
}

inline void EntitySet::add_attribute(std::string attribute_name,
                                     std::vector<double> values, ModelStatus& status)
{
    if (values.size() != count_)
    {
        status.ok = false;
        status.error = "entity set '" + name_ + "': attribute '" + attribute_name +
                       "' has " + std::to_string(values.size()) + " values, expected " +
                       std::to_string(count_);
        return;
    }
    for (auto& attr : attributes_)
    {
        if (attr.first == attribute_name)
        {
            attr.second = std::move(values);
            return;
        }
    }
    attributes_.emplace_back(std::move(attribute_name), std::move(values));
}

inline const std::vector<double>* EntitySet::attribute(std::string_view attribute_name) const
{
    for (auto& attr : attributes_)
        if (attr.first == attribute_name)
            return &attr.second;
    return nullptr;
}

inline bool EntitySet::validate(ModelStatus& status) const
{
    if (name_.empty()) { status.ok = false; status.error = "entity set: name must not be empty"; return false; }
    for (auto& attr : attributes_)
    {
        if (attr.second.size() != count_)
        {
            status.ok = false;
            status.error = "entity set '" + name_ + "': attribute '" + attr.first +
                           "' size mismatch";
            return false;
        }
    }
    return true;
}

} // namespace exd::engine::core
