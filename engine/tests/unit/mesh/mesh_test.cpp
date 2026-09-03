// mesh_test.cpp — spec §9/§10: structured topology, metrics, generation.
#include <exd/engine/mesh/boundary.hpp>
#include <exd/engine/mesh/generation.hpp>
#include <exd/engine/mesh/structured.hpp>
#include <exd/engine/mesh/validation.hpp>
#include <doctest/doctest.h>

using namespace exd::engine;
using namespace exd::engine::mesh;

TEST_CASE("StructuredGrid: metrics and neighbors")
{
    StructuredGrid g;
    g.origin = {0, 0, 0};
    g.spacing = {0.1, 0.2, 0.4};
    g.dims = {11, 6, 4};
    ModelStatus st;
    REQUIRE(g.validate(st));

    CHECK(g.node_count() == 11 * 6 * 4ull);
    CHECK(g.cell_count() == 10 * 5 * 3ull);
    CHECK(g.cell_volume() == doctest::Approx(0.008).epsilon(1e-12));
    CHECK(g.face_area(Axis::X) == doctest::Approx(0.08).epsilon(1e-12));
    CHECK(g.face_area(Axis::Y) == doctest::Approx(0.04).epsilon(1e-12));
    CHECK(g.face_area(Axis::Z) == doctest::Approx(0.02).epsilon(1e-12));

    const auto n = g.face_normal(Axis::X, true);
    CHECK(n[0] == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(n[1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(n[2] == doctest::Approx(0.0).epsilon(1e-12));
    const auto nm = g.face_normal(Axis::Z, false);
    CHECK(nm[2] == doctest::Approx(-1.0).epsilon(1e-12));

    // neighbor at the boundary returns -1
    CHECK(g.neighbor(0, 0, 0, Axis::X, false) == -1);
    CHECK(g.neighbor(10, 0, 0, Axis::X, true) == -1);
    CHECK(g.neighbor(5, 2, 1, Axis::Y, true) == 3);
    CHECK(g.neighbor(5, 2, 1, Axis::X, false) == 4);
}

TEST_CASE("Mesh generation: uniform refinement preserves linear fields")
{
    const StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {5, 5, 5});
    ModelStatus st;
    REQUIRE(g.validate(st));

    // linear ramp field
    const size_t n = g.node_count();
    std::vector<double> v(n);
    for (int k = 0; k < 5; ++k)
        for (int j = 0; j < 5; ++j)
            for (int i = 0; i < 5; ++i)
                v[static_cast<size_t>(i) + 5ull * (static_cast<size_t>(j) + 5ull * k)] =
                    0.25 * i + 0.1 * j + 0.05 * k;

    const auto rv = refine_values(g, v, 2);
    const StructuredGrid rg = refine(g, 2);
    CHECK(rg.dims[0] == 9);
    CHECK(rv.size() == rg.node_count());

    // trilinear resampling of a linear field is exact
    for (int k = 0; k < 9; ++k)
        for (int j = 0; j < 9; ++j)
            for (int i = 0; i < 9; ++i)
            {
                // phi(x,y,z) = x + 0.4y + 0.2z (source spacing 0.25); refined
                // spacing 0.125 → value = 0.125i + 0.05j + 0.025k
                const double expect = 0.125 * i + 0.05 * j + 0.025 * k;
                const size_t idx = static_cast<size_t>(i) + 9ull * (static_cast<size_t>(j) + 9ull * k);
                CHECK(rv[idx] == doctest::Approx(expect).epsilon(1e-12));
            }
}

TEST_CASE("Mesh validation: boundary completeness")
{
    ModelStatus st;
    StructuredGrid g = make_structured_grid({0, 0, 0}, {1, 1, 1}, {5, 5, 5});
    BoundaryCondition in;
    in.name = "inlet";
    in.kind = BoundaryConditionKind::Inlet;
    std::array<const BoundaryCondition*, 6> faces{};
    faces[static_cast<size_t>(BoundaryId::XNeg)] = &in;
    // incomplete → rejected
    CHECK(!boundary_complete(g, faces, st));
    CHECK(!st.ok);

    ModelStatus st2;
    BoundaryCondition out;
    out.name = "outlet";
    out.kind = BoundaryConditionKind::Outlet;
    BoundaryCondition wall;
    wall.name = "wall";
    wall.kind = BoundaryConditionKind::Wall;
    for (int b = 0; b < 6; ++b)
        faces[static_cast<size_t>(b)] = (b == 0) ? &in : (b == 1 ? &out : &wall);
    CHECK(boundary_complete(g, faces, st2));
    CHECK(st2.ok);
}
