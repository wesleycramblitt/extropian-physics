// csv_series_test.cpp
// Machine-state CSV time series: header, rows, validation.

#include <exd/physics/io/series_writer.hpp>

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <string>
#include <vector>

using namespace exd::physics::io;
using exd::physics::ModelStatus;

namespace
{

std::string temp_path(const std::string& name)
{
    auto base = std::filesystem::temp_directory_path();
    auto d = base / ("exd_csv_test_" + std::to_string(::getpid()));
    std::error_code ec;
    std::filesystem::create_directories(d, ec);
    return (d / name).string();
}

std::string read_all(const std::string& path)
{
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

} // anonymous namespace

TEST_CASE("csv series: writes header and rows with 9-significant-digit formatting")
{
    const std::string path = temp_path("series1.csv");
    ModelStatus status;
    CsvSeriesWriter w(path, {"omega", "angle_rad", "torque"}, false, &status);
    REQUIRE(status.ok);
    REQUIRE(w);

    CHECK(w.write_row(0.0, {1.5, 0.25, 10.0}));
    CHECK(w.write_row(0.004, {1.75, 0.257, 9.9}));
    w.close();

    const std::string content = read_all(path);
    const std::string expected =
        "time,omega,angle_rad,torque\n"
        "0,1.5,0.25,10\n"
        "0.004,1.75,0.257,9.9\n";
    CHECK(content == expected);
}

TEST_CASE("csv series: wrong column count is rejected")
{
    const std::string path = temp_path("series2.csv");
    CsvSeriesWriter w(path, {"a", "b"}, false, nullptr);
    REQUIRE(w);
    CHECK_FALSE(w.write_row(0.0, {1.0}));     // wrong count
    CHECK_FALSE(w.write_row(0.0, {1.0, 2.0, 3.0})); // wrong count
    CHECK(w.write_row(1.0, {1.0, 2.0}));      // correct
    CHECK(w.write_row(1.0, {3.0, 4.0}));      // correct
    w.close();
    CHECK_FALSE(w.write_row(2.0, {1.0, 2.0})); // closed
}

TEST_CASE("csv series: invalid column names fail construction")
{
    ModelStatus status;
    CsvSeriesWriter bad1("unused.csv", {"a,b"}, false, &status);
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());

    ModelStatus status2;
    CsvSeriesWriter bad2("unused.csv", {"a\nb"}, false, &status2);
    CHECK_FALSE(status2.ok);
}

TEST_CASE("csv series: vector<double> convenience overload")
{
    const std::string path = temp_path("series3.csv");
    CsvSeriesWriter w(path, {"x"}, false, nullptr);
    REQUIRE(w);
    CHECK(w.write_row(0.0, std::vector<double>{42.5}));
    w.close();
    CHECK(read_all(path) == "time,x\n0,42.5\n");
}
