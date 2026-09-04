#pragma once

// Public solver interface for the 3D incompressible FDM CFD solver (fdm3).
//
// FDM3Solver is a persistent-state stepper: initialize() sets up a grid and
// step() advances the SAME grid each call (unlike the legacy 2D step_fdm,
// which re-initializes the field every call).  solve_fdm3() is the
// one-shot convenience wrapper.

#include "fdm3_config.hpp"
#include "fdm3_result.hpp"
#include <exd/engine/output/field_writer.hpp>
#include <exd/engine/output/output_policy.hpp>
#include <exd/engine/core/model_status.hpp>

#include <memory>
#include <span>
#include <vector>

namespace exd::engine::physics::fluid::fdm3 {

/// Persistent-state 3D incompressible Navier-Stokes solver using finite
/// differences (collocated cell-centered variables, one ghost layer, SOR
/// pressure-projection).
class FDM3Solver
{
public:
    FDM3Solver();
    ~FDM3Solver();

    /// Allocate the grid and apply the initial conditions from `config`.
    /// On failure `status.ok` is false and `status.error` describes why.
    bool initialize(const FDM3Config& config, ModelStatus& status);

    /// Advance the solver by one time step.  `dt` is the requested step;
    /// when `config.adaptive_dt` is set the solver clamps it by the CFL
    /// limit (with a small floor).  The same grid is advanced each call.
    bool step(double dt, ModelStatus& status);

    /// Simulated time after the most recent completed step (s).
    double time() const { return time_; }

    /// Number of completed steps.
    int step_count() const { return step_count_; }

    /// The configuration this solver was initialized with.
    const FDM3Config& config() const { return config_; }

    /// Mutable field state for per-step external manipulation (immersed
    /// solids, actuator disks, HITL edits): apply the change, then `step`.
    /// Accessing the mutable adapter marks the state dirty; `step()` ingests
    /// it into the grid (a hard error on an adapter/grid size mismatch) and
    /// clears the flag — untouched adapters cost nothing.
    FDM3FieldData& field();

    /// Current field state (cell-centered values, updated each step).
    const FDM3FieldData& field() const { return field_; }

    /// Diagnostics of the most recent completed step.
    const FDM3StepResult& last_step() const { return last_step_; }

    /// Set per-cell body force (acceleration, m/s², one value per interior
    /// cell: nx*ny*nz ordered i + nx*(j + ny*k)).  The force is the
    /// acceleration experienced by the FLUID (i.e. force per unit mass).
    ///
    /// Sign convention: a positive value of, e.g., fx means the fluid is
    /// accelerated in +x.  A turbine disk decelerates the flow, so the force
    /// applied to the fluid is the NEGATED blade force per unit mass
    /// (fluid receives -F_blade/rho after dividing by cell mass).  This
    /// module only implements the convention; the coupling layer (Wave 4)
    /// is responsible for producing the correctly-signed force field.
    ///
    /// Passing empty spans clears the body force.
    bool set_body_force(std::span<const double> fx, std::span<const double> fy,
                        std::span<const double> fz, ModelStatus& status);

private:
    FDM3Config config_;
    FDM3FieldData field_;
    bool field_dirty_ = false;
    FDM3StepResult last_step_;
    double time_ = 0.0;
    int step_count_ = 0;

    // Internal grid is accessed through the private implementation header.
    struct GridImpl;
    std::unique_ptr<GridImpl> grid_;

    bool body_force_active_ = false;
    std::vector<double> body_fx_;
    std::vector<double> body_fy_;
    std::vector<double> body_fz_;
};

/// Solve 3D incompressible Navier-Stokes using finite differences.
///
/// Constructs a solver, initializes it, runs up to `config.max_steps`
/// (honoring adaptive dt), records a history entry every step, and returns
/// the result with a `converged` flag from the windowed residual checks.
FDM3Result solve_fdm3(const FDM3Config& config);
/// Driver: run the simulation AND stamp fields (velocity + pressure, exd-fld
/// convention) at a cadence. `writer`/`scheduler` are non-owning; when
/// `scheduler` is null the stamp interval falls back to
/// `config.field_stamp_interval`. Without a writer this is equivalent to
/// solve_fdm3() — solvers stay pure, drivers own the stamps.
FDM3Result run_fdm3_simulation(const FDM3Config& config,
                               exd::engine::output::IFieldWriter* writer = nullptr,
                               exd::engine::output::OutputScheduler* scheduler = nullptr,
                               ModelStatus* status = nullptr);


} // namespace exd::engine::physics::fluid::fdm3