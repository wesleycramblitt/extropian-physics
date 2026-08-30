// csv_series.cpp
// Append-style CSV time-series writer for machine states.

#include <exd/physics/io/series_writer.hpp>

#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace exd::physics::io {

namespace {

/// Column names must be simple identifiers: no comma, quote,
/// CR or LF (keeps CSV trivial to parse — consumer contract).
bool valid_column_name(const std::string& c)
{
    if (c.empty()) return false;
    for (char ch : c)
        if (ch == ',' || ch == '"' || ch == '\n' || ch == '\r') return false;
    return true;
}

} // anonymous namespace

struct CsvSeriesWriter::Impl
{
    std::FILE* file = nullptr;
    std::size_t columns = 0;
    bool flush_each_row = false;
};

CsvSeriesWriter::CsvSeriesWriter(std::string path,
                                 std::vector<std::string> columns,
                                 bool flush_each_row,
                                 ModelStatus* status)
    : impl_(new Impl()), ok_(false)
{
    if (columns.empty() || !impl_)
    {
        if (status) { status->ok = false; status->error = "csv series: no columns"; }
        return;
    }
    for (const auto& c : columns)
    {
        if (!valid_column_name(c))
        {
            if (status)
            {
                status->ok = false;
                status->error = "csv series: invalid column name '" + c + "'";
            }
            return;
        }
    }
    impl_->file = std::fopen(path.c_str(), "w");
    if (!impl_->file)
    {
        if (status) { status->ok = false; status->error = "csv series: cannot open '" + path + "'"; }
        return;
    }
    impl_->columns = columns.size();
    impl_->flush_each_row = flush_each_row;
    ok_ = true;

    // Header row: time,<columns...>
    std::fprintf(impl_->file, "time");
    for (const auto& c : columns) std::fprintf(impl_->file, ",%s", c.c_str());
    std::fprintf(impl_->file, "\n");
    if (flush_each_row) std::fflush(impl_->file);
}

CsvSeriesWriter::~CsvSeriesWriter()
{
    close();
    delete impl_;
}

CsvSeriesWriter::CsvSeriesWriter(CsvSeriesWriter&& other) noexcept
    : impl_(other.impl_), ok_(other.ok_)
{
    other.impl_ = nullptr;
    other.ok_ = false;
}

CsvSeriesWriter& CsvSeriesWriter::operator=(CsvSeriesWriter&& other) noexcept
{
    if (this != &other)
    {
        close();
        delete impl_;
        impl_ = other.impl_;
        ok_ = other.ok_;
        other.impl_ = nullptr;
        other.ok_ = false;
    }
    return *this;
}

void CsvSeriesWriter::close()
{
    if (impl_ && impl_->file)
    {
        std::fflush(impl_->file);
        std::fclose(impl_->file);
        impl_->file = nullptr;
    }
    ok_ = false;
}

bool CsvSeriesWriter::write_row(double t, std::span<const double> values)
{
    if (!ok_ || !impl_ || !impl_->file) return false;
    if (values.size() != impl_->columns) return false;
    std::fprintf(impl_->file, "%.9g", t);
    for (double v : values) std::fprintf(impl_->file, ",%.9g", v);
    std::fprintf(impl_->file, "\n");
    if (impl_->flush_each_row) std::fflush(impl_->file);
    return true;
}

bool CsvSeriesWriter::write_row(double t, const std::vector<double>& values)
{
    return write_row(t, std::span<const double>(values.data(), values.size()));
}

} // namespace exd::physics::io
