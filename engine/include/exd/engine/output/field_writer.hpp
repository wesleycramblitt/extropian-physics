#pragma once

// ─────────────────────────────────────────────────────
// Field output channels.
//
// Format-agnostic stamping interface for 3D scalar and
// vector fields plus containers behind it. The default
// container is the compact binary "exd-fld" format
// (see docs/output_channels.md) consumed by external
// renderers/animation apps; other writers (VTK, HDF5…)
// implement IFieldWriter without touching solvers.
//
// Writers are LEAF infrastructure: solvers (solve_*,
// step_*) stay pure and never write; only application
// drivers (simulate_turbine, simulate_engine,
// run_coupled_turbine…) own writer stamps.
// ─────────────────────────────────────────────────────

#include <exd/engine/core/model_status.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>

namespace exd::engine::output {

// ─────────────────────────────────────────────────────
// Field geometry
// ─────────────────────────────────────────────────────

/// Description of a regular (structured) field sample layout.
/// Convention: `dims` are SAMPLE POINT COUNTS per axis and
/// `origin` is the location OF THE FIRST SAMPLE. For
/// cell-centered data (FDM), origin = cell center of the
/// first cell = (dx/2, dy/2, dz/2) offset from the domain
/// corner. Consumers must not assume node-centered layout.
struct FieldGeometry
{
    std::array<double, 3> origin = {0.0, 0.0, 0.0};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<uint32_t, 3> dims = {0, 0, 0};

    /// True when all dims > 0 and all spacing > 0.
    [[nodiscard]] bool valid() const;
};

// ─────────────────────────────────────────────────────
// Writer interface
// ─────────────────────────────────────────────────────

/// One stamp = one point in (simulation) time with one or
/// more named scalar/vector fields. Snapshots are
/// copy-on-write: the writer consumes the data during
/// write_* calls; a storage thread may be added later.
class IFieldWriter
{
public:
    virtual ~IFieldWriter() = default;

    virtual std::string_view name() const = 0;

    /// Begin a new stamp. Returns false when the stamp
    /// cannot be opened (I/O error). Between begin_stamp
    /// and end_stamp, write_scalar_field/write_vector_field
    /// may be called in any order.
    virtual bool begin_stamp(double t, uint64_t step) = 0;

    /// Scalar field: `data` has prod(dims) float32 values,
    /// row-major x fastest (i + nx*(j + ny*k)).
    virtual bool write_scalar_field(std::string_view name,
                                    const FieldGeometry& geo,
                                    std::span<const float> data) = 0;

    /// Vector field: `data` has 3*prod(dims) values,
    /// interleaved xyz, same layout as (x,y,z) of the
    /// scalar case.
    virtual bool write_vector_field(std::string_view name,
                                    const FieldGeometry& geo,
                                    std::span<const float> data) = 0;

    /// Finish the current stamp. Returns false when
    /// finalization (manifest append) fails.
    virtual bool end_stamp() = 0;
};

// ─────────────────────────────────────────────────────
// Config + factory
// ─────────────────────────────────────────────────────

struct FldWriterConfig
{
    std::string directory = "output";   // created if missing
    bool overwrite = true;              // stamp files may overwrite
    bool sanitize_nonfinite = true;     // NaN/Inf -> 0.0f before write
};

/// Factory for the binary "exd-fld" container writer.
/// Returns nullptr with `status.error` set when the
/// directory cannot be created/opened. A non-finite value
/// reported via `status.warnings` if sanitization engages.
std::unique_ptr<IFieldWriter> make_fld_writer(const FldWriterConfig& config,
                                              ModelStatus& status);

/// No-op writer: accepts everything, writes nothing.
/// Default sink for drivers that do not request output.
std::unique_ptr<IFieldWriter> make_null_writer();

} // namespace exd::engine::output
