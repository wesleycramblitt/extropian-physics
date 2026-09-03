// Static field solver tests: parallel-plate capacitor (electrostatic →
// E field and linear potential) and an infinite-wire-like current line
// (magnetostatic → Ampère's-law B field), plus validation and result sanity.

#include <exd/engine/physics/electromagnetics/static_fields.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace exd::engine;
using namespace exd::engine::coupling;
using namespace exd::engine::physics::electromagnetics;

namespace
{

// ── Fixtures ──────────────────────────────────────────────────────

// Parallel-plate capacitor: plates at x = ±0.1, 100 V (left) and 0 V (right),
// spanning the full y/z cross-section.  The domain is wide enough in y/z that
// the grounded box walls do not drag the central potential down.
StaticFieldConfig capacitor_config()
{
    StaticFieldConfig cfg;
    cfg.mode = StaticFieldMode::Electrostatic;
    cfg.dims = {48, 49, 49};
    cfg.spacing = {0.005, 0.01, 0.01};
    cfg.face_values = {0, 0, 0, 0, 0, 0};
    StaticFieldConfig::BoxPatch plate_left;
    plate_left.center = {-0.1, 0.0, 0.0};
    plate_left.half_extents = {0.001, 1.0, 1.0};
    plate_left.value = 100.0;
    StaticFieldConfig::BoxPatch plate_right;
    plate_right.center = {0.1, 0.0, 0.0};
    plate_right.half_extents = {0.001, 1.0, 1.0};
    plate_right.value = 0.0;
    cfg.patches = {plate_left, plate_right};
    return cfg;
}

// Current line along z at the box center.  The box half-extent in z spans the
// whole domain so the source is uniform in z (an effective infinite wire);
// the x/y half-extents equal half a cell, so exactly one node per z-plane
// carries the source and the discrete current is I = jz·(2·hx)(2·hy) = 40 A.
StaticFieldConfig wire_config()
{
    StaticFieldConfig cfg;
    cfg.mode = StaticFieldMode::Magnetostatic;
    cfg.dims = {65, 65, 65};
    cfg.spacing = {0.002, 0.002, 0.002};
    cfg.face_values = {0, 0, 0, 0, 0, 0};
    StaticFieldConfig::CurrentBox wire;
    wire.center = {0.0, 0.0, 0.0};
    wire.half_extents = {0.001, 0.001, 1.0};
    wire.jz = 1e7;
    cfg.currents = {wire};
    return cfg;
}

constexpr double MU0 = 1.25663706212e-6;               // H/m
constexpr double PI = 3.14159265358979323846;

} // anonymous namespace

// ── Electrostatic capacitor ───────────────────────────────────────

TEST_CASE("parallel-plate capacitor gives uniform E and linear potential")
{
    auto config = capacitor_config();
    auto result = solve_static_field(config);
    REQUIRE(result.ok);
    REQUIRE(result.warnings.empty());

    auto E = make_vector_grid_field(result.field_vector);
    REQUIRE(E != nullptr);

    // E = −∇φ points from the high-potential plate (x = −0.1) to the grounded
    // plate (x = +0.1), so E_x ≈ +ΔV/d = +100/0.2 = +500 V/m in the gap.
    std::array<double, 3> e{0.0, 0.0, 0.0};
    REQUIRE(E->sample({0.0, 0.0, 0.0}, e));
    CHECK(e[0] == doctest::Approx(500.0).epsilon(0.10));
    CHECK(std::fabs(e[1]) < 1e-3);
    CHECK(std::fabs(e[2]) < 1e-3);

    // Potential is linear in x between the plates: φ(x = −0.05) ≈ 75 V.
    auto potential = make_scalar_grid_field(result.potential);
    REQUIRE(potential != nullptr);
    double phi_value = 0.0;
    REQUIRE(potential->sample({-0.05, 0.0, 0.0}, phi_value));
    CHECK(phi_value == doctest::Approx(75.0).epsilon(0.05));
}

// ── Magnetostatic wire ────────────────────────────────────────────

TEST_CASE("infinite wire field matches Ampere's law B = mu0 I / (2 pi r)")
{
    auto config = wire_config();
    auto result = solve_static_field(config);
    REQUIRE(result.ok);
    REQUIRE(result.warnings.empty());

    const double I = 1e7 * (2.0 * 0.001) * (2.0 * 0.001); // jz·(2hx)(2hy) = 40 A
    const double r = 0.02;                                // m
    const double B_expected = MU0 * I / (2.0 * PI * r);   // ≈ 4e-4 T

    auto B = make_vector_grid_field(result.field_vector);
    REQUIRE(B != nullptr);

    // Far-field magnitude on the x and y axes (15% band absorbs the finite
    // box and the discrete source).
    for (const std::array<double, 3>& p :
         std::array<std::array<double, 3>, 2>{{{r, 0.0, 0.0}, {0.0, r, 0.0}}})
    {
        std::array<double, 3> b{0.0, 0.0, 0.0};
        REQUIRE(B->sample(p, b));
        const double magnitude =
            std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        CHECK(magnitude == doctest::Approx(B_expected).epsilon(0.15));
        // A = A_z·ẑ has no x/y component, so B_z ≡ 0.
        CHECK(std::fabs(b[2]) < 1e-8);
    }
}

// ── Validation ────────────────────────────────────────────────────

TEST_CASE("static field solver rejects invalid configurations")
{
    SUBCASE("dims below 2 per axis")
    {
        StaticFieldConfig cfg;
        cfg.dims = {1, 5, 5};
        auto result = solve_static_field(cfg);
        CHECK_FALSE(result.ok);
        CHECK_FALSE(result.error.empty());
    }

    SUBCASE("invalid SOR relaxation")
    {
        StaticFieldConfig cfg;
        cfg.sor_omega = 2.5;
        auto result = solve_static_field(cfg);
        CHECK_FALSE(result.ok);
        CHECK_FALSE(result.error.empty());
    }

    SUBCASE("non-positive spacing")
    {
        StaticFieldConfig cfg;
        cfg.spacing = {0.0, 0.01, 0.01};
        auto result = solve_static_field(cfg);
        CHECK_FALSE(result.ok);
        CHECK_FALSE(result.error.empty());
    }
}

// ── Result sanity ─────────────────────────────────────────────────

TEST_CASE("static field result carries consistent grids and iteration data")
{
    auto config = capacitor_config();
    auto result = solve_static_field(config);
    REQUIRE(result.ok);

    CHECK(result.iterations > 0);
    CHECK(result.residual > 0.0);
    CHECK(result.potential.dims == config.dims);
    CHECK(result.potential.spacing == config.spacing);
    CHECK(result.field_vector.dims == config.dims);
    CHECK(result.field_vector.spacing == config.spacing);

    const std::size_t node_count = static_cast<std::size_t>(config.dims[0]) *
                                   config.dims[1] * config.dims[2];
    CHECK(result.potential.values.size() == node_count);
    CHECK(result.field_vector.values.size() == 3 * node_count);
}