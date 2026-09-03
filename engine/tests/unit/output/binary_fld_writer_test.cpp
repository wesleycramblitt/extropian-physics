// binary_fld_writer_test.cpp
// Round-trip verification of the exd-fld v1 container via a
// minimal spec-faithful reader (the same contract the
// animation repo must implement).

#include <exd/engine/output/field_writer.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>

using namespace exd::engine::output;
using exd::engine::core::ModelStatus;

namespace
{

constexpr size_t HEADER_BYTES = 32;
constexpr size_t FIELD_TABLE_BYTES = 164;

/// Minimal spec reader: returns false on any contract violation.
struct ParsedField
{
    std::string name;
    uint8_t kind = 0;
    uint8_t element_type = 0;
    uint32_t count = 0;
    uint32_t payload_offset = 0;
    uint32_t payload_bytes = 0;
    std::array<double, 3> origin{};
    std::array<double, 3> spacing{};
    std::array<uint32_t, 3> dims{};
    std::vector<float> payload;
};

struct ParsedStamp
{
    double t = 0.0;
    uint32_t field_count = 0;
    std::vector<ParsedField> fields;
};

bool parse_fld(const std::string& path, ParsedStamp& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    uint8_t h[HEADER_BYTES];
    f.read(reinterpret_cast<char*>(h), HEADER_BYTES);
    if (std::memcmp(h, "EXDFLD01", 8) != 0) return false;   // magic
    uint32_t version, fcount;
    std::memcpy(&version, h + 8, 4);
    std::memcpy(&fcount, h + 12, 4);
    if (version != 1) return false;
    double t;
    std::memcpy(&t, h + 16, 8);
    out.t = t;
    out.field_count = fcount;

    out.fields.resize(fcount);
    for (uint32_t i = 0; i < fcount; ++i)
    {
        uint8_t ft[FIELD_TABLE_BYTES];
        f.read(reinterpret_cast<char*>(ft), FIELD_TABLE_BYTES);
        ParsedField& pf = out.fields[i];
        pf.name.assign(reinterpret_cast<char*>(ft), strnlen(reinterpret_cast<char*>(ft), 64));
        pf.kind = ft[64];
        pf.element_type = ft[65];
        std::memcpy(&pf.count, ft + 68, 4);
        std::memcpy(&pf.payload_offset, ft + 72, 4);
        std::memcpy(&pf.payload_bytes, ft + 76, 4);
        std::memcpy(pf.origin.data(), ft + 80, 24);
        std::memcpy(pf.spacing.data(), ft + 104, 24);
        std::memcpy(pf.dims.data(), ft + 128, 12);
        if (pf.element_type != 0) return false; // float32 only in v1
        if (pf.payload_bytes != pf.count * 4u) return false;
    }
    for (auto& pf : out.fields)
    {
        f.seekg(pf.payload_offset, std::ios::beg);
        pf.payload.resize(pf.count);
        f.read(reinterpret_cast<char*>(pf.payload.data()), pf.payload_bytes);
    }
    return true;
}

std::string temp_dir()
{
    auto base = std::filesystem::temp_directory_path();
    auto d = base / ("exd_fld_test_" + std::to_string(::getpid()));
    std::error_code ec;
    std::filesystem::remove_all(d, ec);
    std::filesystem::create_directories(d, ec);
    return d.string();
}

} // anonymous namespace

TEST_CASE("fld writer: round-trip of scalar + vector fields")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(status.ok);
    REQUIRE(writer);

    FieldGeometry geo;
    geo.origin = {0.5, 0.25, 0.125};
    geo.spacing = {0.1, 0.2, 0.3};
    geo.dims = {2, 3, 4}; // 24 scalar samples

    std::vector<float> scalar(24);
    std::vector<float> vector(24 * 3);
    for (size_t i = 0; i < scalar.size(); ++i)
    {
        scalar[i] = static_cast<float>(i) * 0.5f;
        vector[3 * i + 0] = static_cast<float>(i);
        vector[3 * i + 1] = -static_cast<float>(i);
        vector[3 * i + 2] = 100.0f + static_cast<float>(i);
    }

    REQUIRE(writer->begin_stamp(12.5, 7));
    REQUIRE(writer->write_scalar_field("pressure", geo, scalar));
    REQUIRE(writer->write_vector_field("velocity", geo, vector));
    REQUIRE(writer->end_stamp());

    ParsedStamp stamp;
    const std::string path = dir + "/step_00000007.fld";
    REQUIRE(parse_fld(path, stamp));
    CHECK(stamp.t == doctest::Approx(12.5));
    REQUIRE(stamp.field_count == 2);
    REQUIRE(stamp.fields.size() == 2);

    const auto& ps = stamp.fields[0];
    CHECK(ps.name == "pressure");
    CHECK(ps.kind == 0);
    CHECK(ps.element_type == 0);
    REQUIRE(ps.count == 24);
    CHECK(ps.dims == std::array<uint32_t, 3>{2, 3, 4});
    CHECK(ps.origin[0] == doctest::Approx(0.5));
    CHECK(ps.spacing[2] == doctest::Approx(0.3));
    for (uint32_t i = 0; i < ps.count; ++i)
        CHECK(ps.payload[i] == doctest::Approx(static_cast<double>(i) * 0.5));

    const auto& pv = stamp.fields[1];
    CHECK(pv.name == "velocity");
    CHECK(pv.kind == 1);
    REQUIRE(pv.count == 72);
    CHECK(pv.payload[0] == doctest::Approx(0.0));
    CHECK(pv.payload[1] == doctest::Approx(0.0));
    CHECK(pv.payload[2] == doctest::Approx(100.0));
    CHECK(pv.payload[71] == doctest::Approx(123.0));

    // Timeline manifest: one line "t step filename".
    std::ifstream tl(dir + "/timeline.txt");
    REQUIRE(tl.good());
    std::string line;
    std::getline(tl, line);
    CHECK(line.find("12.5") != std::string::npos);
    CHECK(line.find("00000007") != std::string::npos);
    CHECK_FALSE(std::getline(tl, line)); // exactly one line
}

TEST_CASE("fld writer: NaN/Inf payloads are sanitized to zero")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(writer);

    FieldGeometry geo;
    geo.dims = {2, 1, 1};
    geo.spacing = {1, 1, 1};
    std::vector<float> v = {1.0f, std::numeric_limits<float>::quiet_NaN()};

    REQUIRE(writer->begin_stamp(0.0, 0));
    REQUIRE(writer->write_scalar_field("p", geo, v));
    REQUIRE(writer->end_stamp());

    ParsedStamp stamp;
    REQUIRE(parse_fld(dir + "/step_00000000.fld", stamp));
    REQUIRE(stamp.fields.size() == 1);
    CHECK(stamp.fields[0].payload[0] == doctest::Approx(1.0));
    CHECK(stamp.fields[0].payload[1] == doctest::Approx(0.0));
}

TEST_CASE("fld writer: duplicate field names fail the stamp")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(writer);

    FieldGeometry geo;
    geo.dims = {1, 1, 1};
    geo.spacing = {1, 1, 1};
    std::vector<float> one = {1.0f};

    REQUIRE(writer->begin_stamp(0.0, 1));
    REQUIRE(writer->write_scalar_field("p", geo, one));
    REQUIRE(writer->write_scalar_field("p", geo, one));
    CHECK_FALSE(writer->end_stamp());
}

TEST_CASE("fld writer: invalid geometry or mismatched size is rejected")
{
    const std::string dir = temp_dir();
    ModelStatus status;
    auto writer = make_fld_writer({dir}, status);
    REQUIRE(writer);

    FieldGeometry geo;
    geo.dims = {2, 2, 2};
    geo.spacing = {1, 1, 1};
    std::vector<float> v(8);
    REQUIRE(writer->begin_stamp(0.0, 2));
    CHECK_FALSE(writer->write_scalar_field("p", geo, std::span<const float>(v.data(), 7)));
    CHECK(writer->end_stamp()); // stamp with zero fields still allowed
}

TEST_CASE("fld writer: bad directory fails the factory cleanly")
{
    ModelStatus status;
    auto writer = make_fld_writer({"/proc/definitely/not/creatable/exd_test", true, true}, status);
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
    CHECK(writer == nullptr);
}

TEST_CASE("null writer: accepts everything")
{
    auto writer = make_null_writer();
    CHECK(writer->name() == "null");
    CHECK(writer->begin_stamp(0.0, 0));
    FieldGeometry geo;
    geo.dims = {1, 1, 1};
    geo.spacing = {1, 1, 1};
    const float one_f = 1.0f;
    const std::array<float, 3> v3 = {1.0f, 2.0f, 3.0f};
    CHECK(writer->write_scalar_field("p", geo, std::span<const float>(&one_f, 1)));
    CHECK(writer->write_vector_field("v", geo, std::span<const float>(v3.data(), 3)));
    CHECK(writer->end_stamp());
}
