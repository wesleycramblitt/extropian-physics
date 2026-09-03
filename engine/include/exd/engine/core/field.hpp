#pragma once

// ─────────────────────────────────────────────────────
// Fields (implementation_spec §6, §3.1).
//
// A Field is a named, contiguous array of numerical data
// associated with a physical quantity, carrying the full
// metadata set the spec requires: name, rank, components,
// units, location, association, precision, domain, mesh,
// boundary information.  Storage is data-oriented (one
// flat buffer — never object-per-cell), CPU-resident
// today with residency metadata for the GPU path.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/memory.hpp>
#include <exd/engine/core/model_status.hpp>
#include <exd/engine/core/units.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace exd::engine::core {

enum class FieldRank : uint8_t { Scalar, Vector, Tensor };

enum class FieldLocation : uint8_t
{
    Cell,        // cell-centered
    Face,        // face-centered
    Node,        // node-centered
    Edge,        // edge-centered
    Particle,    // per-particle
    RigidBody,   // per-rigid-body
    Quadrature,  // quadrature points
    Global,      // single value (0-D)
};

enum class FieldPrecision : uint8_t { F32, F64 };

/// Full metadata set for a field (spec §6).
struct FieldMetadata
{
    std::string name;                // e.g. "temperature"
    FieldRank rank = FieldRank::Scalar;
    uint8_t components = 1;          // 1 scalar, 3 vector, 6/9 tensor
    Unit units;                      // dimensional metadata (§25)
    FieldLocation location = FieldLocation::Cell;
    FieldPrecision precision = FieldPrecision::F64;
    std::string association;         // semantic role, e.g. "fluid.velocity"
    std::string domain;              // owning module, e.g. "thermal"
    std::string mesh;                // mesh/domain name ("" = global)
    bool on_boundary = false;        // boundary information flag
    std::string boundary_name;       // boundary patch when on_boundary

    [[nodiscard]] bool is_valid(ModelStatus& status) const;
};

/// The numerical storage stays data-oriented: one contiguous buffer,
/// `size()` elements, `components` per entity.  F64 today; F32 precision
/// is a fidelity switch on the same layout (metadata records it).
class Field
{
public:
    Field() = default;
    Field(FieldMetadata meta, size_t size);

    const FieldMetadata& metadata() const { return meta_; }
    FieldMetadata& metadata() { return meta_; }

    size_t size() const { return data_.size(); }          // data points (size == n_entities * components)
    size_t entity_count() const { return meta_.components ? size() / meta_.components : 0; }
    bool empty() const { return data_.empty(); }

    std::span<double> data() { return data_; }
    std::span<const double> data() const { return data_; }

    double at(size_t i) const { return data_[i]; }
    void set(size_t i, double v)
    {
        data_[i] = v;
        ++version_;
    }

    void assign(double v) { std::fill(data_.begin(), data_.end(), v); ++version_; }

    /// Change tracking (§5: versioning/change tracking).
    uint64_t version() const { return version_; }

    MemorySpace space() const { return space_; }
    void set_space(MemorySpace s) { space_ = s; }         // residency metadata (§5/§39)

    /// Interpret the flat buffer as entities of `components` values each.
    std::span<const double> entity(size_t e) const;
    std::span<double> entity(size_t e);

    [[nodiscard]] bool validate(ModelStatus& status) const;

private:
    FieldMetadata meta_;
    std::vector<double> data_;
    MemorySpace space_ = MemorySpace::Cpu;
    uint64_t version_ = 0;
};

// ── inline implementations ──────────────────────────

inline bool FieldMetadata::is_valid(ModelStatus& status) const
{
    if (name.empty()) { status.ok = false; status.error = "field: name must not be empty"; return false; }
    if (components == 0 ||
        (rank == FieldRank::Scalar && components != 1) ||
        (rank == FieldRank::Vector && components != 3) ||
        (rank == FieldRank::Tensor && components != 6 && components != 9))
    {
        status.ok = false;
        status.error = "field: rank/components mismatch (" + name + ")";
        return false;
    }
    return true;
}

inline Field::Field(FieldMetadata meta, size_t size)
    : meta_(std::move(meta)), data_(size)
{
    if (size && size % meta_.components != 0)
    {
        data_.clear();
        // caller-facing error handled by validate(); keep construction cheap
    }
}

inline std::span<const double> Field::entity(size_t e) const
{
    const size_t c = meta_.components ? meta_.components : 1;
    const size_t off = e * c;
    return (off + c <= data_.size()) ? std::span<const double>(data_.data() + off, c)
                                     : std::span<const double>{};
}

inline std::span<double> Field::entity(size_t e)
{
    const size_t c = meta_.components ? meta_.components : 1;
    const size_t off = e * c;
    return (off + c <= data_.size()) ? std::span<double>(data_.data() + off, c)
                                     : std::span<double>{};
}

inline bool Field::validate(ModelStatus& status) const
{
    if (!meta_.is_valid(status)) return false;
    if (!data_.empty() && data_.size() % meta_.components != 0)
    {
        status.ok = false;
        status.error = "field: size " + std::to_string(data_.size()) +
                       " not a multiple of components " + std::to_string(meta_.components);
        return false;
    }
    return true;
}

} // namespace exd::engine::core
