#pragma once

// ─────────────────────────────────────────────────────
// Memory architecture (implementation_spec §59, §40).
//
// Fields remain resident in one MemorySpace; transfers are
// explicit.  Hot numerical paths avoid dynamic allocation:
// temporaries belong to caller-owned Buffers / arenas
// (exd::core provides pool allocators — the engine wraps
// them where needed).  GPU residency is declared now,
// implemented by the CUDA backend (spec Phase 11).
// ─────────────────────────────────────────────────────

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::core {

enum class MemorySpace : uint8_t
{
    Cpu,
    Gpu,   // reserved (CUDA/HIP backend, spec Phase 11)
};

constexpr const char* to_string(MemorySpace s)
{
    return s == MemorySpace::Cpu ? "cpu" : "gpu";
}

/// Owning contiguous byte buffer in a memory space.  For CPU residency the
/// storage is a plain vector; a GPU backend attaches device storage.
/// Copy/synchronize are explicit (spec §39: no per-timestep host↔device
/// round trips).
class Buffer
{
public:
    Buffer() = default;
    Buffer(MemorySpace space, size_t bytes)
        : space_(space), size_(bytes)
    {
        if (space_ == MemorySpace::Cpu) cpu_.resize(bytes);
    }

    MemorySpace space() const { return space_; }
    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }

    std::byte* data() { return cpu_.data(); }
    const std::byte* data() const { return cpu_.data(); }

    /// Device pointer handle — null until a GPU backend attaches storage.
    void* device_handle() const { return device_; }
    void set_device_handle(void* h) { device_ = h; }

    bool copy_from(const Buffer& src);          // host copy (device path: roadmap)
    bool copy_to(Buffer& dst) const { return dst.copy_from(*this); }

private:
    MemorySpace space_ = MemorySpace::Cpu;
    size_t size_ = 0;
    std::vector<std::byte> cpu_;
    void* device_ = nullptr;                    // GPU backend-owned
};

inline bool Buffer::copy_from(const Buffer& src)
{
    if (src.size_ != size_ || src.space_ != space_)
        return false;
    if (space_ == MemorySpace::Cpu)
        std::copy(src.cpu_.begin(), src.cpu_.end(), cpu_.begin());
    // GPU copies are backend operations (roadmap).
    return true;
}

} // namespace exd::engine::core
