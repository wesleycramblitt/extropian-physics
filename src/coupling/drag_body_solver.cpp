// drag_body_solver.cpp
// Partitioned FSI-lite driver (W12): 6-DOF rigid body <-> fdm3 through a
// Gaussian-smeared point force.  See the header for the model contract.

#include <exd/physics/coupling/drag_body_solver.hpp>
#include <exd/physics/coupling/field_sampler.hpp>
#include <exd/physics/mechanics/rigid_body.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace exd::physics::coupling {

namespace
{

constexpr double kBlobSupport = 4.0; // ±4ε cutoff; exp(-8) ~ 3.3e-4

bool fail(ModelStatus& status, DragBodyResult& result, const std::string& msg)
{
    status.ok = false;
    status.error = msg;
    result.ok = false;
    result.status = status;
    return false;
}

/// Gaussian-smeared momentum source with DISCRETE normalization Σw_i = 1
/// over the ±4ε support.  The exact sum is what makes the injected fluid
/// momentum equal the applied force.
/// `force_action` = the force applied TO THE FLUID (i.e. -F_drag on the body).
/// Out-of-domain cells are clipped, which breaks the identity -- the driver
/// warns when the body is within kBlobSupport·ε of a wall.
void smear_point_force(std::vector<double>& fx, std::vector<double>& fy,
                       std::vector<double>& fz, const fluid::fdm3::FDM3Config& flow,
                       const std::array<double, 3>& position,
                       const std::array<double, 3>& force_action, double eps,
                       double rho)
{
    const int nx = flow.nx, ny = flow.ny, nz = flow.nz;
    const double dx = flow.lx / static_cast<double>(nx);
    const double dy = flow.ly / static_cast<double>(ny);
    const double dz = flow.lz / static_cast<double>(nz);

    std::fill(fx.begin(), fx.end(), 0.0);
    std::fill(fy.begin(), fy.end(), 0.0);
    std::fill(fz.begin(), fz.end(), 0.0);

    // Pass 1: raw weights over the support box.
    double weight_sum = 0.0;
    for (int k = 0; k < nz; ++k)
    {
        const double z = (static_cast<double>(k) + 0.5) * dz;
        const double dzc = z - position[2];
        if (std::fabs(dzc) > kBlobSupport * eps) continue;
        for (int j = 0; j < ny; ++j)
        {
            const double y = (static_cast<double>(j) + 0.5) * dy;
            const double dyc = y - position[1];
            if (std::fabs(dyc) > kBlobSupport * eps) continue;
            for (int i = 0; i < nx; ++i)
            {
                const double x = (static_cast<double>(i) + 0.5) * dx;
                const double dxc = x - position[0];
                if (std::fabs(dxc) > kBlobSupport * eps) continue;
                const double r2 = dxc * dxc + dyc * dyc + dzc * dzc;
                weight_sum += std::exp(-r2 / (2.0 * eps * eps));
            }
        }
    }
    if (!(weight_sum > 0.0)) return; // body outside the grid: nothing to smear

    // Pass 2: per-cell acceleration = w_i · F_action / (ρ · V_cell).
    const double inv_mass = 1.0 / (rho * dx * dy * dz);
    for (int k = 0; k < nz; ++k)
    {
        const double z = (static_cast<double>(k) + 0.5) * dz;
        const double dzc = z - position[2];
        if (std::fabs(dzc) > kBlobSupport * eps) continue;
        for (int j = 0; j < ny; ++j)
        {
            const double y = (static_cast<double>(j) + 0.5) * dy;
            const double dyc = y - position[1];
            if (std::fabs(dyc) > kBlobSupport * eps) continue;
            for (int i = 0; i < nx; ++i)
            {
                const double x = (static_cast<double>(i) + 0.5) * dx;
                const double dxc = x - position[0];
                if (std::fabs(dxc) > kBlobSupport * eps) continue;
                const double r2 = dxc * dxc + dyc * dyc + dzc * dzc;
                const double w = std::exp(-r2 / (2.0 * eps * eps)) / weight_sum;
                const std::size_t idx = static_cast<std::size_t>(i) +
                                        static_cast<std::size_t>(nx) *
                                            (static_cast<std::size_t>(j) +
                                             static_cast<std::size_t>(ny) *
                                                 static_cast<std::size_t>(k));
                fx[idx] = w * force_action[0] * inv_mass;
                fy[idx] = w * force_action[1] * inv_mass;
                fz[idx] = w * force_action[2] * inv_mass;
            }
        }
    }
}

/// Sum of fluid momentum Σ ρ·u·dV from the solver field (cell-centered,
/// collocated: cell-center values are the cell averages).
std::array<double, 3> fluid_momentum(const fluid::fdm3::FDM3Solver& solver)
{
    const auto& f = solver.field();
    const fluid::fdm3::FDM3Config& cfg = solver.config();
    const double rho = cfg.rho;
    const double dV = (cfg.lx / static_cast<double>(cfg.nx)) *
                      (cfg.ly / static_cast<double>(cfg.ny)) *
                      (cfg.lz / static_cast<double>(cfg.nz));
    std::array<double, 3> mom{0.0, 0.0, 0.0};
    const double scale = rho * dV;
    for (std::size_t i = 0; i < f.u.size(); ++i)
    {
        mom[0] += f.u[i] * scale;
        mom[1] += f.v[i] * scale;
        mom[2] += f.w[i] * scale;
    }
    return mom;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Smearing (public; acceptance pins Σ ρ·a·dV = F exactly)
// ---------------------------------------------------------------------------

void apply_smeared_point_force(std::vector<double>& fx, std::vector<double>& fy,
                               std::vector<double>& fz,
                               const fluid::fdm3::FDM3Config& flow,
                               const std::array<double, 3>& position,
                               const std::array<double, 3>& force, double eps,
                               double rho)
{
    smear_point_force(fx, fy, fz, flow, position, force, eps, rho);
}

// ---------------------------------------------------------------------------
// Drag law
// ---------------------------------------------------------------------------

std::array<double, 3> body_drag_force(const std::array<double, 3>& relative_velocity,
                                      double rho, double drag_area)
{
    const double speed = std::sqrt(relative_velocity[0] * relative_velocity[0] +
                                   relative_velocity[1] * relative_velocity[1] +
                                   relative_velocity[2] * relative_velocity[2]);
    const double coeff = -0.5 * rho * drag_area * speed;
    return {coeff * relative_velocity[0],
            coeff * relative_velocity[1],
            coeff * relative_velocity[2]};
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

DragBodyResult simulate_drag_body(const DragBodyConfig& config, ModelStatus& status)
{
    DragBodyResult result;
    status = ModelStatus{};

    // ── Validation ──────────────────────────────────────────────────────
    if (!(config.mass > 0.0))
        return fail(status, result, "drag body: mass must be > 0"), result;
    if (!(config.drag_area >= 0.0))
        return fail(status, result, "drag body: drag_area must be >= 0"), result;
    for (int c = 0; c < 3; ++c)
    {
        if (!(config.inertia_principal[static_cast<size_t>(c)] > 0.0))
            return fail(status, result,
                        "drag body: inertia_principal must be positive"), result;
    }
    if (config.fluid_steps_per_exchange < 1)
        return fail(status, result,
                    "drag body: fluid_steps_per_exchange must be >= 1"), result;
    if (!(config.force_relaxation > 0.0 && config.force_relaxation <= 1.0))
        return fail(status, result,
                    "drag body: force_relaxation must be in (0, 1]"), result;
    if (!(config.smear_cells >= 0.5))
        return fail(status, result, "drag body: smear_cells must be >= 0.5"), result;
    if (!(config.sample_lead >= 0.0))
        return fail(status, result, "drag body: sample_lead must be >= 0"), result;
    if (config.flow_precondition_steps < 1 && config.flow_precondition_steps != 0)
        return fail(status, result, "drag body: flow_precondition_steps must be >= 0"), result;
    if (config.max_steps < 1)
        return fail(status, result, "drag body: max_steps must be >= 1"), result;
    if (config.flow.adaptive_dt)
        return fail(status, result,
                    "drag body: flow.adaptive_dt must be false (body window = "
                    "cadence · flow.dt)"), result;
    if (!(config.flow.dt > 0.0))
        return fail(status, result, "drag body: flow.dt must be > 0"), result;
    if (config.flow.nx < 1 || config.flow.ny < 1 || config.flow.nz < 1)
        return fail(status, result, "drag body: flow grid must be nonempty"), result;

    const double dx = config.flow.lx / static_cast<double>(config.flow.nx);
    const double dy = config.flow.ly / static_cast<double>(config.flow.ny);
    const double dz = config.flow.lz / static_cast<double>(config.flow.nz);
    const double hmin = std::min({dx, dy, dz});
    const double eps = config.smear_cells * hmin;

    // The ±4ε blob must fit inside the cell-center domain at the start so
    // that the discrete normalization (Σw_i = 1) is exact.
    const std::array<double, 3> cell_lo = {0.5 * dx, 0.5 * dy, 0.5 * dz};
    const std::array<double, 3> cell_hi = {config.flow.lx - 0.5 * dx,
                                           config.flow.ly - 0.5 * dy,
                                           config.flow.lz - 0.5 * dz};
    // The blob AND the upstream probe must fit inside the cell-center domain
    // at the start (probe direction fallback: upstream of the -z gravity
    // convention when at rest, i.e. +z).
    const double lead = config.sample_lead * eps;
    std::array<double, 3> start_probe = config.initial_position;
    {
        std::array<double, 3> v_rel = config.initial_velocity; // v_f = 0 at rest
        const double sp = std::sqrt(v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] +
                                    v_rel[2] * v_rel[2]);
        if (sp > 1e-9)
        {
            for (int c = 0; c < 3; ++c)
                start_probe[static_cast<size_t>(c)] -=
                    lead * v_rel[static_cast<size_t>(c)] / sp;
        }
        else
        {
            start_probe[2] += lead; // fallback: probe above the body
        }
    }
    for (int c = 0; c < 3; ++c)
    {
        const double lo = cell_lo[static_cast<size_t>(c)];
        const double hi = cell_hi[static_cast<size_t>(c)];
        const double pb = config.initial_position[static_cast<size_t>(c)];
        if (pb - kBlobSupport * eps < lo || pb + kBlobSupport * eps > hi)
            return fail(status, result,
                        "drag body: initial position must place the ±4ε smear "
                        "blob inside the fluid domain"), result;
        const double pp = start_probe[static_cast<size_t>(c)];
        if (pp < lo || pp > hi)
            return fail(status, result,
                        "drag body: initial position must keep the upstream "
                        "probe inside the fluid domain"), result;
    }

    // ── Setup ───────────────────────────────────────────────────────────
    fluid::fdm3::FDM3Solver solver;
    if (!solver.initialize(config.flow, status))
        return fail(status, result, "drag body: " + status.error), result;

    mechanics::RigidBodyConfig rb_cfg;
    rb_cfg.mass = config.mass;
    rb_cfg.inertia_principal = config.inertia_principal;
    auto dynamics = mechanics::make_rigid_body_dynamics(
        mechanics::RigidBodyIntegration::SymplecticEuler, rb_cfg);

    mechanics::RigidBodyState body;
    body.position = config.initial_position;
    body.linear_velocity = config.initial_velocity;
    // orientation = identity, angular velocity = 0 (no torque in v1)

    std::vector<double> fx(static_cast<size_t>(config.flow.nx) *
                           static_cast<size_t>(config.flow.ny) *
                           static_cast<size_t>(config.flow.nz), 0.0);
    std::vector<double> fy(fx.size(), 0.0);
    std::vector<double> fz(fx.size(), 0.0);

    const double window_dt = static_cast<double>(config.fluid_steps_per_exchange) *
                             config.flow.dt;
    const double rho = config.flow.rho;
    const double beta = config.force_relaxation;
    std::array<double, 3> applied_fluid = {0.0, 0.0, 0.0}; // β-relaxed (prev)

    bool warned_wall = false;
    bool warned_stab = false;
    bool warned_oob = false;
    const std::array<double, 3> gravity = config.gravity;

    result.probe_time.reserve(config.max_steps);
    result.probe_position.reserve(config.max_steps);
    result.probe_velocity.reserve(config.max_steps);
    result.probe_fluid_momentum.reserve(config.max_steps);

    std::array<double, 3> last_drag = {0.0, 0.0, 0.0};
    std::array<double, 3> last_vf = {0.0, 0.0, 0.0};

    // ── Precondition: establish the flow before the body is released ────
    {
        ModelStatus flow_status;
        for (uint64_t s = 0; s < config.flow_precondition_steps; ++s)
        {
            if (!solver.step(config.flow.dt, flow_status))
                return fail(status, result, "drag body: precondition step: " +
                                            flow_status.error), result;
        }
    }

    // ── Main loop: exchanges ────────────────────────────────────────────
    for (uint64_t exchange = 0; exchange < config.max_steps; ++exchange)
    {
        // 1. Sample the fluid velocity AHEAD of the body (upstream probe:
        //    excludes the body's own smeared reaction -- see the header).
        std::unique_ptr<IFlowField3D> adapter = make_fdm3_field_adapter(solver);
        {
            const double sp = std::sqrt(body.linear_velocity[0] *
                                            body.linear_velocity[0] +
                                        body.linear_velocity[1] *
                                            body.linear_velocity[1] +
                                        body.linear_velocity[2] *
                                            body.linear_velocity[2]);
            std::array<double, 3> probe = body.position;
            if (sp > 1e-9)
            {
                for (int c = 0; c < 3; ++c)
                    probe[static_cast<size_t>(c)] -=
                        lead * body.linear_velocity[static_cast<size_t>(c)] / sp;
            }
            else
            {
                probe[2] += lead; // at rest: probe above the body
            }
            std::array<double, 3> vel{0.0, 0.0, 0.0};
            double p = 0.0;
            if (!adapter->sample(probe, vel, p))
            {
                if (!warned_oob)
                {
                    status.warnings.push_back(
                        "drag body: upstream probe outside the fluid domain; "
                        "v_f = 0");
                    warned_oob = true;
                }
                vel = {0.0, 0.0, 0.0};
            }
            last_vf = vel;

            // 2. Drag on relative velocity; Newton's third law for the fluid.
            std::array<double, 3> v_rel = {
                body.linear_velocity[0] - vel[0],
                body.linear_velocity[1] - vel[1],
                body.linear_velocity[2] - vel[2],
            };
            last_drag = body_drag_force(v_rel, rho, config.drag_area);

            // 3. β-relaxed action on the fluid, then smear (Σw = 1).
            std::array<double, 3> action = {
                -last_drag[0], -last_drag[1], -last_drag[2],
            };
            std::array<double, 3> applied = {
                beta * action[0] + (1.0 - beta) * applied_fluid[0],
                beta * action[1] + (1.0 - beta) * applied_fluid[1],
                beta * action[2] + (1.0 - beta) * applied_fluid[2],
            };
            applied_fluid = applied;
            smear_point_force(fx, fy, fz, config.flow, body.position, applied,
                              eps, rho);
            if (!solver.set_body_force(fx, fy, fz, status))
                return fail(status, result,
                            "drag body: set_body_force: " + status.error), result;
        }

        // 4. Advance the body one window (loads held constant).
        mechanics::RigidBodyForces loads;
        for (int c = 0; c < 3; ++c)
        {
            loads.force[static_cast<size_t>(c)] =
                config.mass * gravity[static_cast<size_t>(c)] +
                last_drag[static_cast<size_t>(c)];
        }
        ModelStatus body_status;
        body = dynamics->advance(window_dt, loads, body, body_status);
        if (!body_status.ok)
            return fail(status, result, "drag body: body advance: " +
                                        body_status.error), result;

        // 5. Run the fluid window under the smeared source.
        ModelStatus flow_status;
        for (int s = 0; s < config.fluid_steps_per_exchange; ++s)
        {
            if (!solver.step(config.flow.dt, flow_status))
                return fail(status, result, "drag body: fluid step: " +
                                            flow_status.error), result;
        }

        // 6. Probes + guards.
        const double t = solver.time();
        result.probe_time.push_back(t);
        result.probe_position.push_back(body.position);
        result.probe_velocity.push_back(body.linear_velocity);
        result.probe_fluid_momentum.push_back(fluid_momentum(solver));

        const double speed = std::sqrt(body.linear_velocity[0] *
                                           body.linear_velocity[0] +
                                       body.linear_velocity[1] *
                                           body.linear_velocity[1] +
                                       body.linear_velocity[2] *
                                           body.linear_velocity[2]);
        if (!warned_stab && config.flow.dt >= eps / (3.0 * std::max(speed, 1e-3)))
        {
            status.warnings.push_back(
                "drag body: dt is too coarse to resolve the smeared force "
                "layer (dt >= eps/(3·v)); results degrade");
            warned_stab = true;
        }
        if (!warned_wall)
        {
            bool near_wall = false;
            for (int c = 0; c < 3; ++c)
            {
                const double lo = cell_lo[static_cast<size_t>(c)];
                const double hi = cell_hi[static_cast<size_t>(c)];
                const double p = body.position[static_cast<size_t>(c)];
                if (p - kBlobSupport * eps < lo || p + kBlobSupport * eps > hi)
                    near_wall = true;
            }
            if (near_wall)
            {
                status.warnings.push_back(
                    "drag body: body within ±4ε of a boundary; smear blob "
                    "clipped, momentum conservation degrades");
                warned_wall = true;
            }
        }
        {
            // Near-wall check for the upstream probe (explicit, once).
            const double sp = std::sqrt(body.linear_velocity[0] *
                                            body.linear_velocity[0] +
                                        body.linear_velocity[1] *
                                            body.linear_velocity[1] +
                                        body.linear_velocity[2] *
                                            body.linear_velocity[2]);
            std::array<double, 3> probe = body.position;
            if (sp > 1e-9)
            {
                for (int c = 0; c < 3; ++c)
                    probe[static_cast<size_t>(c)] -=
                        lead * body.linear_velocity[static_cast<size_t>(c)] / sp;
            }
            else
            {
                probe[2] += lead;
            }
            for (int c = 0; c < 3; ++c)
            {
                const double lo = cell_lo[static_cast<size_t>(c)];
                const double hi = cell_hi[static_cast<size_t>(c)];
                const double p = probe[static_cast<size_t>(c)];
                if (p < lo || p > hi)
                {
                    status.warnings.push_back(
                        "drag body: upstream probe outside the fluid domain; "
                        "v_f = 0");
                    warned_oob = true;
                    break;
                }
            }
        }
    }

    // ── Close out ───────────────────────────────────────────────────────
    const auto& last_v = result.probe_velocity.back();
    result.terminal_velocity = std::sqrt(last_v[0] * last_v[0] +
                                         last_v[1] * last_v[1] +
                                         last_v[2] * last_v[2]);
    result.drag_at_end = std::sqrt(last_drag[0] * last_drag[0] +
                                   last_drag[1] * last_drag[1] +
                                   last_drag[2] * last_drag[2]);
    result.sampled_fluid_velocity_final = last_vf;
    result.fluid_momentum = fluid_momentum(solver);
    result.exchanges = config.max_steps;
    result.ok = true;
    result.status = status;
    return result;
}

} // namespace exd::physics::coupling
