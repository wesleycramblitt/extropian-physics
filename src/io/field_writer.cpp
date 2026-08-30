// field_writer.cpp
// Binary "exd-fld" container writer + null writer.
// Format spec: docs/output_channels.md.

#include <exd/physics/io/field_writer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <vector>

namespace exd::physics::io {

// ── FieldGeometry ─────────────────────────────────────

bool FieldGeometry::valid() const
{
    if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) return false;
    if (spacing[0] <= 0.0 || spacing[1] <= 0.0 || spacing[2] <= 0.0) return false;
    return std::isfinite(origin[0]) && std::isfinite(origin[1]) && std::isfinite(origin[2]);
}

// ── exd-fld v1 container ──────────────────────────────

namespace {

constexpr char MAGIC[8] = {'E', 'X', 'D', 'F', 'L', 'D', '0', '1'};
constexpr uint32_t VERSION = 1;
constexpr size_t HEADER_BYTES = 32;
constexpr size_t FIELD_TABLE_BYTES = 164;
constexpr size_t MAX_NAME_BYTES = 64;

enum class FieldKind : uint8_t { Scalar = 0, Vector = 1 };

struct PendingField
{
    std::string name;
    FieldKind kind = FieldKind::Scalar;
    FieldGeometry geo;
    std::vector<float> payload; // float32 = the only element type in v1
};

class FldWriter final : public IFieldWriter
{
public:
    FldWriter(FldWriterConfig config, ModelStatus& status)
        : config_(std::move(config))
    {
        std::error_code ec;
        if (!std::filesystem::create_directories(config_.directory, ec) && ec)
        {
            status.ok = false;
            status.error = "fld writer: cannot create directory '" + config_.directory
                           + "': " + ec.message();
            return;
        }
        timeline_path_ = (std::filesystem::path(config_.directory) / "timeline.txt").string();
        if (config_.overwrite)
        {
            FILE* f = std::fopen(timeline_path_.c_str(), "w");
            if (f) std::fclose(f);
        }
    }

    std::string_view name() const override { return "exd-fld-v1"; }

    bool begin_stamp(double t, uint64_t step) override
    {
        fields_.clear();
        t_ = t;
        step_ = step;
        return true;
    }

    bool write_scalar_field(std::string_view name, const FieldGeometry& geo,
                            std::span<const float> data) override
    {
        return append_field(name, FieldKind::Scalar, geo, data, false);
    }

    bool write_vector_field(std::string_view name, const FieldGeometry& geo,
                            std::span<const float> data) override
    {
        return append_field(name, FieldKind::Vector, geo, data, true);
    }

    bool end_stamp() override
    {
        if (timeline_path_.empty()) return false; // failed ctor
        if (!has_duplicate_names())
        {
            // Duplicate field names are a caller contract error; report and skip the stamp.
            has_error_ = true;
            last_error_ = "fld writer: duplicate field name in stamp";
            return false;
        }

        char fname[64];
        std::snprintf(fname, sizeof fname, "step_%08llu.fld",
                      static_cast<unsigned long long>(step_));
        const std::string path =
            (std::filesystem::path(config_.directory) / fname).string();

        // Header + per-field table (each field carries its own geometry).
        std::vector<uint8_t> head(HEADER_BYTES + fields_.size() * FIELD_TABLE_BYTES, 0);
        std::memcpy(head.data(), MAGIC, 8);
        uint32_t fcount = static_cast<uint32_t>(fields_.size());
        std::memcpy(head.data() + 8, &VERSION, 4);
        std::memcpy(head.data() + 12, &fcount, 4);
        std::memcpy(head.data() + 16, &t_, 8);
        // [24..32) reserved

        uint32_t offset = static_cast<uint32_t>(HEADER_BYTES + fcount * FIELD_TABLE_BYTES);
        for (size_t i = 0; i < fields_.size(); ++i)
        {
            const PendingField& f = fields_[i];
            const size_t base = HEADER_BYTES + i * FIELD_TABLE_BYTES;
            std::memcpy(head.data() + base, f.name.c_str(),
                        std::min(f.name.size(), MAX_NAME_BYTES));
            head[base + 64] = static_cast<uint8_t>(f.kind);
            head[base + 65] = 0; // element_type: float32
            uint32_t count = static_cast<uint32_t>(f.payload.size());
            uint32_t bytes = count * 4u;
            std::memcpy(head.data() + base + 68, &count, 4);
            std::memcpy(head.data() + base + 72, &offset, 4);
            std::memcpy(head.data() + base + 76, &bytes, 4);
            std::memcpy(head.data() + base + 80, f.geo.origin.data(), 24);
            std::memcpy(head.data() + base + 104, f.geo.spacing.data(), 24);
            std::memcpy(head.data() + base + 128, f.geo.dims.data(), 12);
            offset += bytes;
        }

        FILE* out = std::fopen(path.c_str(), "wb");
        if (!out)
        {
            has_error_ = true;
            last_error_ = "fld writer: cannot open '" + path + "'";
            return false;
        }
        bool ok = std::fwrite(head.data(), 1, head.size(), out) == head.size();
        for (const PendingField& f : fields_)
        {
            if (!ok) break;
            ok = std::fwrite(f.payload.data(), sizeof(float), f.payload.size(), out)
                 == f.payload.size();
        }
        if (std::fclose(out) != 0) ok = false;

        if (ok)
        {
            FILE* tl = std::fopen(timeline_path_.c_str(), "a");
            if (tl)
            {
                std::fprintf(tl, "%.9g %llu %s\n", t_,
                             static_cast<unsigned long long>(step_), fname);
                std::fclose(tl);
            }
        }
        else
        {
            has_error_ = true;
            last_error_ = "fld writer: write failed for '" + path + "'";
        }
        return ok;
    }

private:
    bool append_field(std::string_view name, FieldKind kind, const FieldGeometry& geo,
                      std::span<const float> data, bool is_vector)
    {
        if (!geo.valid())
        {
            has_error_ = true;
            last_error_ = "fld writer: invalid geometry";
            return false;
        }
        const uint64_t expect = static_cast<uint64_t>(geo.dims[0]) * geo.dims[1] * geo.dims[2]
                                * (is_vector ? 3u : 1u);
        if (data.size() != expect)
        {
            has_error_ = true;
            last_error_ = "fld writer: field '" + std::string(name) + "' size mismatch";
            return false;
        }

        PendingField f;
        f.name = name.substr(0, MAX_NAME_BYTES - 1);
        f.kind = kind;
        f.geo = geo;
        f.payload.assign(data.begin(), data.end());
        if (config_.sanitize_nonfinite)
        {
            for (float& v : f.payload)
            {
                if (!std::isfinite(v)) { v = 0.0f; sanitized_ = true; }
            }
        }
        fields_.push_back(std::move(f));
        return true;
    }

    bool has_duplicate_names() const
    {
        for (size_t i = 0; i < fields_.size(); ++i)
            for (size_t j = i + 1; j < fields_.size(); ++j)
                if (fields_[i].name == fields_[j].name) return false;
        return true;
    }

    FldWriterConfig config_;
    std::string timeline_path_;
    std::vector<PendingField> fields_;
    double t_ = 0.0;
    uint64_t step_ = 0;
    bool sanitized_ = false;
    bool has_error_ = false;
    std::string last_error_;
};

class NullWriter final : public IFieldWriter
{
public:
    std::string_view name() const override { return "null"; }
    bool begin_stamp(double, uint64_t) override { return true; }
    bool write_scalar_field(std::string_view, const FieldGeometry&,
                            std::span<const float>) override { return true; }
    bool write_vector_field(std::string_view, const FieldGeometry&,
                            std::span<const float>) override { return true; }
    bool end_stamp() override { return true; }
};

} // anonymous namespace

std::unique_ptr<IFieldWriter> make_fld_writer(const FldWriterConfig& config,
                                              ModelStatus& status)
{
    auto w = std::make_unique<FldWriter>(config, status);
    if (!status.ok) return nullptr;
    return w;
}

std::unique_ptr<IFieldWriter> make_null_writer()
{
    return std::make_unique<NullWriter>();
}

} // namespace exd::physics::io
