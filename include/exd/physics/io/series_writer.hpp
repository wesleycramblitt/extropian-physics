#pragma once

// ─────────────────────────────────────────────────────
// Machine-state time series output (CSV).
//
// One row per stamp: time + named scalar columns. This is
// the channel for rotor/crank/piston kinematics (omega,
// angle, position, torque, power, temperatures…) consumed
// by dashboards, plots, and the animation repo's part
// motion (linear + angular position over time).
// ─────────────────────────────────────────────────────

#include <exd/physics/model_status.hpp>

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace exd::physics::io {

/// Append-style CSV time series with a header row.
/// One line per write_row call; `flush_each_row` makes a
/// row visible immediately (real-time consumers) at the
/// cost of per-row flush.
class CsvSeriesWriter
{
public:
    /// Opens `path` (truncates if it exists) and writes
    /// the header row `time,<columns...>`.
    /// On failure `status->ok` is set false (status may be
    /// null, in which case failures are silent).
    CsvSeriesWriter(std::string path,
                    std::vector<std::string> columns,
                    bool flush_each_row = false,
                    ModelStatus* status = nullptr);

    explicit operator bool() const { return ok_; }

    /// Write one row: `values` must have exactly
    /// columns.size() entries. Returns false when the file
    /// is not open or the row is rejected.
    bool write_row(double t, std::span<const double> values);

    /// Convenience: doubles only (use for fixed columns).
    bool write_row(double t, const std::vector<double>& values);

    /// Explicit close (also flushes). Idempotent; the
    /// destructor closes as well.
    void close();

    ~CsvSeriesWriter();

    CsvSeriesWriter(const CsvSeriesWriter&) = delete;
    CsvSeriesWriter& operator=(const CsvSeriesWriter&) = delete;
    CsvSeriesWriter(CsvSeriesWriter&&) noexcept;
    CsvSeriesWriter& operator=(CsvSeriesWriter&&) noexcept;

private:
    struct Impl;
    Impl* impl_ = nullptr;
    bool ok_ = false;
};

} // namespace exd::physics::io
