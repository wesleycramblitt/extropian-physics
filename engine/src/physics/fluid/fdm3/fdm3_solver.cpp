#include <exd/engine/physics/fluid/fdm3/fdm3_solver.hpp>
#include "fdm3_internal.hpp"

#include <algorithm>
#include <cmath>

namespace exd::engine::physics::fluid::fdm3 {

struct FDM3Solver::GridImpl : FDM3Grid {};

FDM3Solver::FDM3Solver() = default;
FDM3Solver::~FDM3Solver() = default;

// ─────────────────────────────────────────────────────
// Configuration validation
// ─────────────────────────────────────────────────────

bool FDM3Config::validate(std::string& error, std::vector<std::string>& warnings) const {
    error.clear();
    warnings.clear();

    if (nx < 4 || ny < 4 || nz < 4) { error = "Grid must be at least 4x4x4"; return false; }
    if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
        error = "Domain dimensions must be positive"; return false;
    }
    if (rho <= 0.0) { error = "Density must be positive"; return false; }
    if (mu < 0.0) { error = "Viscosity must be non-negative"; return false; }
    if (dt <= 0.0) { error = "Time step must be positive"; return false; }
    if (max_steps < 1) { error = "Maximum steps must be >= 1"; return false; }
    if (pressure_max_iterations < 1) { error = "Pressure iterations must be >= 1"; return false; }
    if (convergence_window < 1) { error = "Convergence window must be >= 1"; return false; }
    if (sor_omega <= 1.0 || sor_omega >= 2.0) {
        error = "SOR omega must be in (1, 2)"; return false;
    }
    if (velocity_under_relaxation <= 0.0 || velocity_under_relaxation > 1.0) {
        error = "Velocity under-relaxation must be in (0, 1]"; return false;
    }
    if (pressure_under_relaxation <= 0.0 || pressure_under_relaxation > 1.0) {
        error = "Pressure under-relaxation must be in (0, 1]"; return false;
    }
    if (adaptive_dt && cfl_target <= 0.0) {
        error = "CFL target must be positive when adaptive dt is enabled"; return false;
    }

    if (boundary_conditions.empty()) {
        error = "At least one boundary condition required"; return false;
    }

    bool seen[6] = {false, false, false, false, false, false};
    for (const auto& bc : boundary_conditions) {
        const int f = static_cast<int>(bc.face);
        if (f < 0 || f >= 6) { error = "Invalid boundary face"; return false; }
        if (seen[f]) { error = "Duplicate boundary condition for a face"; return false; }
        seen[f] = true;

        if (bc.type == FDMBoundaryType::Periodic && bc.paired_face != natural_opposite(bc.face)) {
            warnings.push_back("Periodic boundary pair is not the natural opposite face; "
                               "only the natural opposite pairing is honored in v1");
        }
    }
    return true;
}

// ─────────────────────────────────────────────────────
// Utility / diagnostics
// ─────────────────────────────────────────────────────

double max_velocity(const FDM3Grid& g) {
    double vmax = 0.0;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                double mag = std::sqrt(g.u[id] * g.u[id] + g.v[id] * g.v[id] +
                                       g.w[id] * g.w[id]);
                if (mag > vmax) vmax = mag;
            }
    return vmax;
}

double max_divergence(const FDM3Grid& g) {
    double dmax = 0.0;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                double div = std::abs(spatial::divergence(g, i, j, k));
                if (div > dmax) dmax = div;
            }
    return dmax;
}

void compute_diagnostics(const FDM3Grid& g, const FDM3Config&, FDM3StepResult& r) {
    double ru = 0.0, rv = 0.0, rw = 0.0, rp = 0.0, vmax = 0.0, dmax = 0.0;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                size_t id = g.idx(i, j, k);
                ru = std::max(ru, std::abs(g.u[id] - g.u_old[id]));
                rv = std::max(rv, std::abs(g.v[id] - g.v_old[id]));
                rw = std::max(rw, std::abs(g.w[id] - g.w_old[id]));
                rp = std::max(rp, std::abs(g.p[id] - g.p_old[id]));
                double mag = std::sqrt(g.u[id] * g.u[id] + g.v[id] * g.v[id] +
                                       g.w[id] * g.w[id]);
                vmax = std::max(vmax, mag);
                dmax = std::max(dmax, std::abs(spatial::divergence(g, i, j, k)));
            }
    r.residual_u = ru;
    r.residual_v = rv;
    r.residual_w = rw;
    r.residual_p = rp;
    r.max_velocity = vmax;
    r.divergence = dmax;
}

// ─────────────────────────────────────────────────────
// Field extraction
// ─────────────────────────────────────────────────────

void extract_field(const FDM3Grid& g, const FDM3Config& config, FDM3FieldData& field) {
    field.nx = g.nx;
    field.ny = g.ny;
    field.nz = g.nz;

    const size_t cells = static_cast<size_t>(g.nx) * g.ny * g.nz;
    field.u.resize(cells);
    field.v.resize(cells);
    field.w.resize(cells);
    field.p.resize(cells);
    field.x.resize(g.nx);
    field.y.resize(g.ny);
    field.z.resize(g.nz);

    for (int i = 0; i < g.nx; ++i) field.x[i] = (i + 0.5) * config.dx();
    for (int j = 0; j < g.ny; ++j) field.y[j] = (j + 0.5) * config.dy();
    for (int k = 0; k < g.nz; ++k) field.z[k] = (k + 0.5) * config.dz();

    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                size_t id = field.index(i, j, k);
                size_t gid = g.idx(i + 1, j + 1, k + 1);
                field.u[id] = g.u[gid];
                field.v[id] = g.v[gid];
                field.w[id] = g.w[gid];
                field.p[id] = g.p[gid];
            }
}

// ─────────────────────────────────────────────────────
// Adaptive time step (CFL clamp with a floor)
// ─────────────────────────────────────────────────────

static void clamp_dt_from_cfl(const FDM3Grid& g, const FDM3Config& config, double& dt) {
    if (!config.adaptive_dt) return;
    const double hmin = std::min(g.dx, std::min(g.dy, g.dz));
    const double vmax = max_velocity(g);
    double dt_cfl = (vmax > 1e-12) ? config.cfl_target * hmin / vmax : config.dt;
    dt_cfl = std::min(dt_cfl, config.dt);
    const double floor = 1e-6 * hmin;
    dt = std::max(dt_cfl, floor);
}

// ─────────────────────────────────────────────────────
// FDM3Solver
// ─────────────────────────────────────────────────────

bool FDM3Solver::initialize(const FDM3Config& config, ModelStatus& status) {
    status.ok = true;
    status.error.clear();
    status.warnings.clear();

    std::string error;
    if (!config.validate(error, status.warnings)) {
        status.ok = false;
        status.error = error;
        return false;
    }

    config_ = config;
    grid_ = std::make_unique<GridImpl>();
    initialize_grid(*grid_, config_);
    apply_boundary_conditions(*grid_, config_);

    time_ = 0.0;
    step_count_ = 0;
    body_force_active_ = false;
    body_fx_.clear();
    body_fy_.clear();
    body_fz_.clear();

    extract_field(*grid_, config_, field_);
    last_step_ = FDM3StepResult{};
    return true;
}

bool FDM3Solver::step(double dt, ModelStatus& status) {
    status.ok = true;
    status.error.clear();
    if (!grid_) {
        status.ok = false;
        status.error = "Solver not initialized";
        return false;
    }

    double dt_eff = dt;
    clamp_dt_from_cfl(*grid_, config_, dt_eff);

    // 1. Save state for SIMPLE
    std::copy(grid_->u.begin(), grid_->u.end(), grid_->u_old.begin());
    std::copy(grid_->v.begin(), grid_->v.end(), grid_->v_old.begin());
    std::copy(grid_->w.begin(), grid_->w.end(), grid_->w_old.begin());
    std::copy(grid_->p.begin(), grid_->p.end(), grid_->p_old.begin());

    // 2. Predictor (explicit time integration on the momentum equations;
    //    the body force is constant during the step and enters EVERY stage
    //    of the integrator's RHS — the correct treatment for all methods;
    //    the old post-step application mismatched the time levels and made
    //    Heun unstable with sustained forces).
    apply_time_integration(*grid_, config_, dt_eff,
                           body_force_active_ ? &body_fx_ : nullptr,
                           body_force_active_ ? &body_fy_ : nullptr,
                           body_force_active_ ? &body_fz_ : nullptr);

    // 3. Pressure-velocity coupling (SIMPLE): Poisson solve for p' with a
    //    zero initial guess; periodic face ghosts are refreshed during the
    //    relaxation (see solve_pressure_poisson).
    compute_pressure_rhs(*grid_, config_, dt_eff);
    std::fill(grid_->p.begin(), grid_->p.end(), 0.0);
    solve_pressure_poisson(*grid_, config_);

    // Stash p' and restore the real pressure field.
    std::copy(grid_->p.begin(), grid_->p.end(), grid_->p_prime.begin());
    std::copy(grid_->p_old.begin(), grid_->p_old.end(), grid_->p.begin());

    // 4-5. Correct velocity (relative to _old) and update pressure.
    correct_velocity(*grid_, config_, dt_eff);

    // 6. Boundary conditions.
    apply_boundary_conditions(*grid_, config_);

    time_ += dt_eff;
    ++step_count_;

    // 7. Diagnostics.
    FDM3StepResult r;
    compute_diagnostics(*grid_, config_, r);
    r.time = time_;
    r.step = step_count_;
    const double hmin = std::min(grid_->dx, std::min(grid_->dy, grid_->dz));
    r.cfl = (hmin > 0.0) ? r.max_velocity * dt_eff / hmin : 0.0;
    last_step_ = r;

    extract_field(*grid_, config_, field_);
    return true;
}

FDM3FieldData& FDM3Solver::field() { return field_; }

bool FDM3Solver::set_body_force(std::span<const double> fx, std::span<const double> fy,
                                std::span<const double> fz, ModelStatus& status) {
    status.ok = true;
    status.error.clear();
    if (!grid_) {
        status.ok = false;
        status.error = "Solver not initialized";
        return false;
    }

    const size_t cells = static_cast<size_t>(config_.nx) * config_.ny * config_.nz;
    if (fx.empty() && fy.empty() && fz.empty()) {
        body_force_active_ = false;
        body_fx_.clear();
        body_fy_.clear();
        body_fz_.clear();
        return true;
    }
    if (fx.size() != cells || fy.size() != cells || fz.size() != cells) {
        status.ok = false;
        status.error = "Body force arrays must each contain nx*ny*nz values";
        return false;
    }

    body_fx_.assign(fx.begin(), fx.end());
    body_fy_.assign(fy.begin(), fy.end());
    body_fz_.assign(fz.begin(), fz.end());
    body_force_active_ = true;
    return true;
}

// ─────────────────────────────────────────────────────
// One-shot solve
// ─────────────────────────────────────────────────────

FDM3Result solve_fdm3(const FDM3Config& config) {
    FDM3Result result;
    FDM3Solver solver;
    ModelStatus status;

    if (!solver.initialize(config, status)) {
        result.valid = false;
        result.error = status.error;
        result.warnings = status.warnings;
        return result;
    }
    result.warnings = status.warnings;

    for (int step = 0; step < config.max_steps; ++step) {
        if (!solver.step(config.dt, status)) {
            result.valid = false;
            result.error = status.error;
            return result;
        }

        const FDM3StepResult& sr = solver.last_step();
        result.history.push_back(sr);

        // Windowed convergence check (mirrors the 2D solver).
        if (step > 0 && sr.step % config.convergence_window == 0) {
            const double max_res = std::max(std::max(sr.residual_u, sr.residual_v),
                                            std::max(sr.residual_w, sr.residual_p));
            if (max_res < config.convergence_tolerance) {
                result.converged = true;
                break;
            }
        }
    }

    result.valid = true;
    result.steps_taken = solver.step_count();
    result.final_time = solver.time();
    result.field = solver.field();
    if (!result.converged)
        result.warnings.push_back("Solver did not converge within max_steps");

    return result;
}

// ─────────────────────────────────────────────────────
// Stamping driver (Wave 7: real runs with field output)
// ─────────────────────────────────────────────────────

FDM3Result run_fdm3_simulation(const FDM3Config& config,
                               exd::engine::output::IFieldWriter* writer,
                               exd::engine::output::OutputScheduler* scheduler,
                               ModelStatus* status)
{
    FDM3Result result;

    FDM3Solver solver;
    {
        std::string error;
        if (!config.validate(error, result.warnings))
        {
            result.valid = false;
            result.error = error;
            if (status) *status = {false, error, result.warnings};
            return result;
        }
        ModelStatus init_status;
        if (!solver.initialize(config, init_status))
        {
            result.valid = false;
            result.error = init_status.error;
            result.warnings.insert(result.warnings.end(),
                                   init_status.warnings.begin(),
                                   init_status.warnings.end());
            if (status) *status = init_status;
            return result;
        }
    }

    bool failed = false;
    for (int step = 0; step < config.max_steps && !failed; ++step)
    {
        ModelStatus step_status;
        if (!solver.step(config.dt, step_status))
        {
            failed = true;
            result.error = step_status.error;
            break;
        }
        if (writer)
        {
            bool emit = false;
            if (scheduler)
            {
                scheduler->set_now(solver.time());
                emit = scheduler->should_emit(static_cast<uint64_t>(step));
            }
            else
            {
                emit = config.field_stamp_interval > 0
                       && (static_cast<uint64_t>(step) % config.field_stamp_interval == 0);
            }
            if (emit)
            {
                const auto& field = solver.field();
                exd::engine::output::FieldGeometry geo;
                geo.origin = {0.5 * config.dx(), 0.5 * config.dy(), 0.5 * config.dz()};
                geo.spacing = {config.dx(), config.dy(), config.dz()};
                geo.dims = {static_cast<uint32_t>(config.nx),
                            static_cast<uint32_t>(config.ny),
                            static_cast<uint32_t>(config.nz)};

                const size_t n = static_cast<size_t>(config.nx) * config.ny * config.nz;
                std::vector<float> vel(3 * n);
                std::vector<float> pres(n);
                for (size_t i = 0; i < n; ++i)
                {
                    vel[3 * i + 0] = static_cast<float>(field.u[i]);
                    vel[3 * i + 1] = static_cast<float>(field.v[i]);
                    vel[3 * i + 2] = static_cast<float>(field.w[i]);
                    pres[i] = static_cast<float>(field.p[i]);
                }
                if (!writer->begin_stamp(solver.time(), static_cast<uint64_t>(step)))
                {
                    result.warnings.push_back("fdm3 driver: field stamp begin failed");
                }
                else
                {
                    writer->write_vector_field("velocity", geo, vel);
                    writer->write_scalar_field("pressure", geo, pres);
                    writer->end_stamp();
                }
            }
        }
    }

    result.valid = !failed;
    result.field = solver.field();
    result.steps_taken = solver.step_count();
    result.final_time = solver.time();
    result.history = {solver.last_step()};
    if (status)
    {
        status->ok = result.valid;
        status->error = result.error;
        status->warnings = result.warnings;
    }
    return result;
}
} // namespace exd::engine::physics::fluid::fdm3
