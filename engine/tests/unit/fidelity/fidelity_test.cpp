// fidelity_test.cpp — spec §45/§46: fidelity profiles.
#include <exd/engine/diagnostics/diagnostics.hpp>
#include <exd/engine/fidelity/profiles.hpp>
#include <exd/engine/mesh/generation.hpp>
#include <doctest/doctest.h>

using namespace exd::engine;
using namespace exd::engine::fidelity;

TEST_CASE("Fidelity profiles: REALTIME … HIGH_FIDELITY monotone strictness")
{
    const auto rt = profile(FidelityLevel::Realtime);
    const auto bal = profile(FidelityLevel::Balanced);
    const auto hf = profile(FidelityLevel::HighFidelity);

    CHECK(rt.solver_tolerance_scale > bal.solver_tolerance_scale);
    CHECK(hf.solver_tolerance_scale < bal.solver_tolerance_scale);
    CHECK(hf.coupling_tolerance < bal.coupling_tolerance);
    CHECK(hf.coupling_iterations > bal.coupling_iterations);
    CHECK(rt.simplified_physics == true);
    CHECK(bal.simplified_physics == false);
    CHECK(rt.precision == core::FieldPrecision::F32);
    CHECK(bal.precision == core::FieldPrecision::F64);
    CHECK(std::string(to_string(FidelityLevel::Accurate)) == "ACCURATE");
}

TEST_CASE("Diagnostics: CFL, conservation, stopwatch")
{
    using namespace exd::engine::diagnostics;
    const mesh::StructuredGrid g = mesh::make_structured_grid({0, 0, 0}, {1, 1, 1}, {11, 11, 11});
    const auto cfl = compute_cfl(g, 0.5, 0.01, 0.01);
    CHECK(cfl.cfl_advection == doctest::Approx(0.05).epsilon(1e-12));
    CHECK(cfl.cfl_diffusion == doctest::Approx(0.01).epsilon(1e-12));
    CHECK(cfl.stable());

    std::vector<double> before = {1.0, 2.0, 3.0};
    std::vector<double> after = {1.05, 1.95, 3.0};
    ConservationDiagnostics cons;
    cons.initial_total = field_total(before);
    cons.final_total = field_total(after);
    CHECK(cons.absolute_error() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(cons.conserved(1e-12));

    Stopwatch sw;
    sw.start();
    CHECK(sw.seconds() >= 0.0);
}
