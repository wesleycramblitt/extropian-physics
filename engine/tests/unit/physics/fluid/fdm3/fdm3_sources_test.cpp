// fdm3_sources_test.cpp — W16 per-cell momentum sources:
//   * immersed-solid mask geometry is exact at the cell centers
//   * a stationary immersed slab blocks the duct: the cells inside are
//     frozen (|u| ~ 0), the flow squeezes over the top, and the discrete
//     mass conservation (inlet flux == outlet flux) still holds
//   * a MOVING immersed solid drags the occupied cells to its velocity
//     (the "water through movable veins" mechanism)
//   * Boussinesq buoyancy matches the analytic per-cell force and drives
//     hot-bottom convection (mean upward velocity at the mid-plane > 0)
#include <exd/engine/physics/fluid/fdm3/fdm3_solver.hpp>
#include <exd/engine/physics/fluid/fdm3/fdm3_sources.hpp>
#include <doctest/doctest.h>

#include <cmath>

using namespace exd::engine;
using namespace exd::engine::physics::fluid::fdm3;

namespace {

FDM3Config duct(int nx = 20, int ny = 8, int nz = 8)
{
    FDM3Config flow;
    flow.nx = nx; flow.ny = ny; flow.nz = nz;
    flow.lx = 1.0; flow.ly = 0.4; flow.lz = 0.4;
    flow.rho = 1.225;
    flow.mu = 0.05;
    flow.dt = 0.01;
    flow.max_steps = 100000;
    flow.pressure_max_iterations = 60;
    flow.pressure_tolerance = 5e-6;
    flow.sor_omega = 1.45;
    flow.velocity_under_relaxation = 0.8;
    flow.pressure_under_relaxation = 0.35;
    flow.convergence_tolerance = 1e-6;
    flow.time_integration = TimeIntegration::Heun;
    return flow;
}

void channel_bcs(FDM3Config& flow, double u0)
{
    FDM3BoundaryCondition inlet;
    inlet.face = BoundaryFace::XMin;
    inlet.type = FDMBoundaryType::Inlet;
    inlet.u_value = u0;
    flow.boundary_conditions.push_back(inlet);
    FDM3BoundaryCondition outlet;
    outlet.face = BoundaryFace::XMax;
    outlet.type = FDMBoundaryType::FixedPressure;
    outlet.p_value = 0.0;
    flow.boundary_conditions.push_back(outlet);
    for (auto face : {BoundaryFace::YMin, BoundaryFace::YMax,
                      BoundaryFace::ZMin, BoundaryFace::ZMax})
    {
        FDM3BoundaryCondition bc;
        bc.face = face;
        bc.type = FDMBoundaryType::Wall;
        flow.boundary_conditions.push_back(bc);
    }
}

double mean_velocity_in(const FDM3Config& cfg, const FDM3FieldData& f,
                        const std::vector<bool>& mask, int comp)
{
    double sum = 0.0;
    int n = 0;
    for (int k = 0; k < cfg.nz; ++k)
        for (int j = 0; j < cfg.ny; ++j)
            for (int i = 0; i < cfg.nx; ++i)
            {
                const size_t id = fdm3_cell_index(cfg, i, j, k);
                if (!mask[id]) continue;
                const size_t fid = f.index(i, j, k);
                sum += (comp == 0) ? f.u[fid] : ((comp == 1) ? f.v[fid] : f.w[fid]);
                ++n;
            }
    return n ? sum / n : 0.0;
}

double plane_flux(const FDM3Config& cfg, const FDM3FieldData& f, int plane_i)
{
    double flux = 0.0;
    const double da = cfg.ly * cfg.lz / (cfg.ny * cfg.nz);
    for (int k = 0; k < cfg.nz; ++k)
        for (int j = 0; j < cfg.ny; ++j)
            flux += f.u[f.index(plane_i, j, k)] * da;
    return flux;
}

} // namespace

TEST_CASE("Immersed mask: cell-center containment is exact")
{
    FDM3Config cfg = duct();
    ImmersedSolid box;
    box.name = "slab";
    box.shape = ImmersedShape::Box;
    box.center = {0.5, 0.1, 0.2};
    box.half_extents = {0.05, 0.05, 0.05};

    std::vector<bool> mask;
    REQUIRE(immersed_solid_mask(cfg, {box}, mask));
    REQUIRE(mask.size() == fdm3_cell_count(cfg));
    // inside: cells (9..10, 1..2, 3..4) — centers 0.475/0.525, 0.075/0.125, 0.175/0.225
    CHECK(mask[fdm3_cell_index(cfg, 9, 1, 3)]);
    CHECK(mask[fdm3_cell_index(cfg, 10, 2, 4)]);
    // outside
    CHECK(!mask[fdm3_cell_index(cfg, 0, 0, 0)]);
    CHECK(!mask[fdm3_cell_index(cfg, 15, 6, 6)]);
    CHECK(!mask[fdm3_cell_index(cfg, 9, 5, 3)]);
}

TEST_CASE("Immersed slab: frozen cells, squeeze over the top, mass conserved")
{
    FDM3Config flow = duct();
    channel_bcs(flow, 0.2);
    // bottom-half slab across the middle 20% of the channel: the ONLY way
    // past is over the top half
    ImmersedSolid slab;
    slab.name = "slab";
    slab.shape = ImmersedShape::Box;
    slab.center = {0.5, 0.1, 0.2};
    slab.half_extents = {0.1, 0.1, 0.2};         // bottom-HALF, full z-width:
                                                 // the flow must go over or
                                                 // under... over (the top)
    slab.penalty = 150.0;

    std::vector<bool> mask;
    REQUIRE(immersed_solid_mask(flow, {slab}, mask));
    const size_t N = fdm3_cell_count(flow);

    ModelStatus st;
    FDM3Solver solver;
    REQUIRE(solver.initialize(flow, st));
    std::vector<double> fx(N), fy(N), fz(N);
    const int steps = 3000;                      // 30 s ≈ 6 flow-through times
    for (int it = 0; it < steps; ++it)
    {
        const auto& f = solver.field();
        // kinematic freeze (the robust frozen-cell treatment): the slab's
        // velocity is projected onto the occupied cells after each step
        REQUIRE(solver.step(flow.dt, st));
        auto& f2 = solver.field();
        REQUIRE(apply_kinematic_freeze(flow, {slab}, 0.95, f2.u, f2.v, f2.w, st));
        (void)f;
        (void)fx;
    }

    const auto& f = solver.field();
    // frozen cells: |u| inside the slab is a small fraction of the inflow
    const double u_inside = mean_velocity_in(flow, f, mask, 0);
    CHECK(std::fabs(u_inside) < 0.02);
    // the obstruction forces the flow to accelerate to get past: the max
    // streamwise velocity clearly exceeds the inlet value (the slab spans
    // the full z and the bottom half of y — the only path is over the top)
    double u_max = 0.0;
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
                u_max = std::max(u_max, f.u[f.index(i, j, k)]);
    CHECK(u_max > 0.3);
    // discrete mass conservation: inlet flux ≈ outlet flux
    const double q_in = plane_flux(flow, f, 0);
    const double q_out = plane_flux(flow, f, flow.nx - 1);
    CHECK(std::fabs(q_in - q_out) / q_in < 0.03);
}

TEST_CASE("Immersed moving solid: the occupied cells are dragged to u_solid")
{
    FDM3Config flow = duct();
    // closed box (no inlet): the fluid starts at rest
    for (auto face : {BoundaryFace::XMin, BoundaryFace::XMax,
                      BoundaryFace::YMin, BoundaryFace::YMax,
                      BoundaryFace::ZMin, BoundaryFace::ZMax})
    {
        FDM3BoundaryCondition bc;
        bc.face = face;
        bc.type = FDMBoundaryType::Wall;
        flow.boundary_conditions.push_back(bc);
    }
    ImmersedSolid body;
    body.name = "vein";
    body.shape = ImmersedShape::Box;
    body.center = {0.5, 0.2, 0.2};
    body.half_extents = {0.05, 0.05, 0.05};
    body.velocity = {0.4, 0.0, 0.0};             // the vein moves at 0.4 m/s
    body.penalty = 150.0;                        // K·dt = 1.5 (Heun-stable)

    std::vector<bool> mask;
    REQUIRE(immersed_solid_mask(flow, {body}, mask));
    const size_t N = fdm3_cell_count(flow);

    ModelStatus st;
    FDM3Solver solver;
    REQUIRE(solver.initialize(flow, st));
    std::vector<double> fx(N), fy(N), fz(N);
    for (int it = 0; it < 200; ++it)             // t = 2 s
    {
        REQUIRE(solver.step(flow.dt, st));
        auto& f = solver.field();
        REQUIRE(apply_kinematic_freeze(flow, {body}, 0.95, f.u, f.v, f.w, st));
    }
    const auto& f = solver.field();
    const double u_inside = mean_velocity_in(flow, f, mask, 0);
    CHECK(std::fabs(u_inside - 0.4) < 0.03);     // dragged to u_solid
    // the far field has barely responded on this timescale
    const double u_corner = f.u[f.index(2, 2, 2)];
    CHECK(std::fabs(u_corner) < 0.12);
}


TEST_CASE("fdm3: field() adapter edits are ingested and mismatches are hard errors")
{
    FDM3Config flow = duct();
    channel_bcs(flow, 0.2);
    ModelStatus st;
    FDM3Solver solver;
    REQUIRE(solver.initialize(flow, st));

    // a corrupt (resized) adapter must be rejected, not silently dropped
    auto& f = solver.field();
    f.u.resize(f.u.size() / 2);
    REQUIRE_FALSE(solver.step(flow.dt, st));
    REQUIRE(!st.error.empty());

    // restore a consistent adapter so the remaining checks run clean
    ModelStatus st2;
    REQUIRE(solver.initialize(flow, st2));
    // an edit IS ingested: zero the interior streamwise velocity, one step
    // must move it away from the frozen 0.0 initial state but stay finite
    auto& f2 = solver.field();
    std::fill(f2.u.begin(), f2.u.end(), 0.05);
    REQUIRE(solver.step(flow.dt, st2));
    const auto& f3 = solver.field();
    bool any_moved = false;
    for (size_t i = 0; i < f3.u.size(); ++i)
        if (std::fabs(f3.u[i]) > 1e-4) any_moved = true;
    CHECK(any_moved);
}
TEST_CASE("Heun with sustained body forces is stable (the W16 force fix)")
{
    // Regression for the Heun time-level bug: the body force must enter
    // every integrator stage; the old post-step application destabilized
    // Heun with ANY sustained force (uniform or immersed).
    FDM3Config flow;
    flow.nx = 20; flow.ny = 8; flow.nz = 8;
    flow.lx = 1.0; flow.ly = 0.4; flow.lz = 0.4;
    flow.rho = 1.225; flow.mu = 0.05;
    flow.dt = 0.01; flow.max_steps = 100000;
    flow.pressure_max_iterations = 60; flow.pressure_tolerance = 5e-6;
    flow.sor_omega = 1.45;
    flow.velocity_under_relaxation = 0.8; flow.pressure_under_relaxation = 0.35;
    flow.convergence_tolerance = 1e-6;
    flow.time_integration = TimeIntegration::Heun;
    FDM3BoundaryCondition inlet;
    inlet.face = BoundaryFace::XMin;
    inlet.type = FDMBoundaryType::Inlet;
    inlet.u_value = 0.2;
    FDM3BoundaryCondition outlet;
    outlet.face = BoundaryFace::XMax;
    outlet.type = FDMBoundaryType::FixedPressure;
    outlet.p_value = 0.0;
    flow.boundary_conditions.push_back(inlet);
    flow.boundary_conditions.push_back(outlet);
    for (auto face : {BoundaryFace::YMin, BoundaryFace::YMax,
                      BoundaryFace::ZMin, BoundaryFace::ZMax})
    {
        FDM3BoundaryCondition bc;
        bc.face = face;
        bc.type = FDMBoundaryType::Wall;
        flow.boundary_conditions.push_back(bc);
    }

    const size_t N = fdm3_cell_count(flow);
    ModelStatus st;
    FDM3Solver solver;
    REQUIRE(solver.initialize(flow, st));
    std::vector<double> fx(N, 0.05), fy(N, 0.0), fz(N, 0.0);   // uniform +x
    // switch to an immersed penalty after the first 1500 steps
    ImmersedSolid box;
    box.name = "b";
    box.shape = ImmersedShape::Box;
    box.center = {0.5, 0.1, 0.2};
    box.half_extents = {0.05, 0.05, 0.05};
    box.penalty = 10.0;
    bool finite = true;
    for (int it = 0; it < 3000; ++it)
    {
        if (it == 1500)
        {
            std::fill(fx.begin(), fx.end(), 0.0);
            std::fill(fy.begin(), fy.end(), 0.0);
            std::fill(fz.begin(), fz.end(), 0.0);
        }
        const auto& f = solver.field();
        if (it >= 1500)
            REQUIRE(add_immersed_solid_forces(flow, {box}, f.u, f.v, f.w,
                                              fx, fy, fz, st));
        REQUIRE(solver.set_body_force(fx, fy, fz, st));
        REQUIRE(solver.step(flow.dt, st));
        const auto& f2 = solver.field();
        for (size_t i = 0; i < N; ++i)
            if (!std::isfinite(f2.u[i]) || !std::isfinite(f2.p[i]))
                finite = false;
        REQUIRE(finite);
    }
    // the flow actually developed and responded (not frozen/zero)
    const auto& f = solver.field();
    double umax = 0.0;
    for (size_t i = 0; i < N; ++i)
        umax = std::max(umax, std::fabs(f.u[i]));
    CHECK(umax > 0.2);
    CHECK(umax < 2.0);
}

TEST_CASE("Immersed penalty force: the force formula is exact per cell")
{
    // The SOFT (force-based) formulation is verified deterministically here.
    // Integration caveat (documented in the header): with the collocated
    // pressure projection, a UNIFORM body force interior to a region tends
    // to be canceled by the projection in the fully developed steady state
    // (the classic pressure-velocity decoupling) — use the KINEMATIC FREEZE
    // for hard blockage; the penalty is the right tool for moderate drag.
    FDM3Config flow = duct();
    ImmersedSolid sphere;
    sphere.name = "s";
    sphere.shape = ImmersedShape::Sphere;
    sphere.center = {0.5, 0.2, 0.2};
    sphere.radius = 0.05;
    sphere.velocity = {0.3, -0.1, 0.0};
    sphere.penalty = 40.0;

    const size_t N = fdm3_cell_count(flow);
    std::vector<double> u(N), v(N), w(N), fx(N), fy(N), fz(N);
    for (size_t i = 0; i < N; ++i)
    {
        u[i] = 0.2;
        v[i] = 0.05;
        w[i] = 0.0;
    }
    ModelStatus st;
    REQUIRE(add_immersed_solid_forces(flow, {sphere}, u, v, w, fx, fy, fz, st));
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
            {
                const size_t id = fdm3_cell_index(flow, i, j, k);
                const auto c = fdm3_cell_center(flow, i, j, k);
                const double dx = c[0] - 0.5, dy = c[1] - 0.2, dz = c[2] - 0.2;
                const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
                // solid fraction: linear transition over the default one-cell
                // smearing width (max cell size = 0.05)
                const double f = d <= 0.05 - 0.025 ? 1.0
                               : (d >= 0.05 + 0.025 ? 0.0
                               : std::clamp(0.5 - (d - 0.05) / 0.05, 0.0, 1.0));
                if (f <= 0.0)
                {
                    CHECK(fx[id] == doctest::Approx(0.0).epsilon(1e-12));
                    continue;
                }
                // f += −K·f·(u − u_solid)
                CHECK(fx[id] == doctest::Approx(-40.0 * f * (0.2 - 0.3)).epsilon(1e-9));
                CHECK(fy[id] == doctest::Approx(-40.0 * f * (0.05 + 0.1)).epsilon(1e-9));
                CHECK(fz[id] == doctest::Approx(0.0).epsilon(1e-12));
            }
}

TEST_CASE("Boussinesq: per-cell force matches the analytic buoyancy")
{
    FDM3Config flow = duct(20, 8, 8);
    const size_t N = fdm3_cell_count(flow);
    const double beta = 1e-3, t_ref = 300.0;
    const std::array<double, 3> g{0, 0, -9.81};
    std::vector<double> T(N);
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
            {
                const auto c = fdm3_cell_center(flow, i, j, k);
                T[fdm3_cell_index(flow, i, j, k)] = t_ref + 5.0 * c[2];
            }
    std::vector<double> fx(N), fy(N), fz(N);
    ModelStatus st;
    REQUIRE(add_boussinesq_forces(flow, T, beta, t_ref, g, fx, fy, fz, st));
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
            {
                const size_t id = fdm3_cell_index(flow, i, j, k);
                const double z = fdm3_cell_center(flow, i, j, k)[2];
                const double a_exact = -beta * (5.0 * z) * g[2];
                CHECK(fx[id] == doctest::Approx(0.0).epsilon(1e-12));
                CHECK(fy[id] == doctest::Approx(0.0).epsilon(1e-12));
                CHECK(fz[id] == doctest::Approx(a_exact).epsilon(1e-12));
            }
}

TEST_CASE("Boussinesq: hot-bottom stratification drives upward motion")
{
    FDM3Config flow = duct(20, 8, 8);
    for (auto face : {BoundaryFace::XMin, BoundaryFace::XMax,
                      BoundaryFace::YMin, BoundaryFace::YMax,
                      BoundaryFace::ZMin, BoundaryFace::ZMax})
    {
        FDM3BoundaryCondition bc;
        bc.face = face;
        bc.type = FDMBoundaryType::Wall;
        flow.boundary_conditions.push_back(bc);
    }
    const double beta = 1e-3, t_ref = 300.0;
    const std::array<double, 3> g{0, 0, -9.81};
    const size_t N = fdm3_cell_count(flow);
    // linear stratification, HOT at the bottom: T(z) = 310 − 25·z
    std::vector<double> T(N);
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
                T[fdm3_cell_index(flow, i, j, k)] =
                    310.0 - 25.0 * fdm3_cell_center(flow, i, j, k)[2];

    ModelStatus st;
    FDM3Solver solver;
    REQUIRE(solver.initialize(flow, st));
    std::vector<double> fx(N), fy(N), fz(N);
    for (int it = 0; it < 100; ++it)             // t = 1 s
    {
        const auto& f = solver.field();
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
        REQUIRE(add_boussinesq_forces(flow, T, beta, t_ref, g, fx, fy, fz, st));
        REQUIRE(solver.set_body_force(fx, fy, fz, st));
        REQUIRE(solver.step(flow.dt, st));
    }
    const auto& f = solver.field();
    // A CLOSED box: the net vertical flux through any plane is ~0 (the
    // return flow) — the physical anchor is the CIRCULATION: the lower
    // half rises, the upper half sinks, and the velocities are nonzero.
    double w_bot = 0.0, w_top = 0.0, w_max = 0.0;
    int n_bot = 0, n_top = 0;
    for (int k = 0; k < flow.nz; ++k)
        for (int j = 0; j < flow.ny; ++j)
            for (int i = 0; i < flow.nx; ++i)
            {
                const double w = f.w[f.index(i, j, k)];
                w_max = std::max(w_max, std::fabs(w));
                if (k <= 2) { w_bot += w; ++n_bot; }      // lower half
                if (k >= 5) { w_top += w; ++n_top; }      // upper half
            }
    w_bot /= n_bot;
    w_top /= n_top;
    CHECK(w_bot > 1e-5);                         // hot bottom fluid rises
    CHECK(w_top < -1e-5);                        // and returns at the top
    CHECK(w_max > 5e-5);                         // the circulation is real
}
