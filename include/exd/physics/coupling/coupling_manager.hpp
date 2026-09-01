#pragma once

#include <exd/physics/coupling/field_channels.hpp>
#include <exd/physics/coupling/surface_mapping.hpp>
#include <exd/physics/model_status.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace exd::physics::coupling {

// -----------------------------------------------------
// Coupling manager (Phase H-lite): the real exchange of
// field values between solver domains.  Domains register
// channel callbacks (read + sink) so the manager never
// owns solver state.  Links read a source channel at a
// probe point and write the (possibly relaxed) value into
// a target domain's sink.
//
// Explicit exchange applies each due link once (staggered
// operator splitting).  Implicit exchange sub-iterates the
// due links to a tolerance with per-link relaxation.
// -----------------------------------------------------

/// Directional coupling link: reads one field channel from the source domain
/// and writes into the target domain's sink at a set of probe points.
struct CouplingLink
{
    std::string id;                                  // unique link identifier
    std::string source_domain;
    std::string source_channel;                      // registered channel name (read side)
    std::string target_domain;
    std::string target_channel;                      // registered sink name (write side)
    bool vector_field = false;                       // false = scalar link, true = vector link
    InterpolationMode interpolation = InterpolationMode::Nearest;
    double interval = 0.0;                           // seconds between exchanges (0 = every step)
    double relaxation = 1.0;                         // under-relaxation in (0, 1]
    int sub_iterations = 1;                          // implicit sub-iterations per exchange (>= 1)
    std::vector<std::array<double, 3>> probe_points; // where to sample / write
};

class CouplingManager
{
public:
    /// Domain registration uses callbacks so the manager never owns solver
    /// state.  All callbacks are optional (nullptr default) except the read /
    /// write side actually used by a link.
    struct DomainHandle
    {
        std::string name;
        std::function<const IScalarField3D*(std::string_view)> scalar_channel = nullptr;
        std::function<const IVectorField3D*(std::string_view)> vector_channel = nullptr;
        // Sinks receive the value the manager computed for the probe point
        // (nearest-node semantics are implemented by the domain).
        std::function<bool(std::string_view, const std::array<double, 3>&, double)>
            scalar_write = nullptr;
        std::function<bool(std::string_view, const std::array<double, 3>&,
                           const std::array<double, 3>&)>
            vector_write = nullptr;
        std::function<double()> current_time = nullptr; // for interval gating (informational:
                                                        // the manager times links from the `t`
                                                        // passed to exchange*)
    };

    /// Register a domain.  Fails on an empty or duplicate name.
    bool register_domain(const DomainHandle& domain, exd::physics::ModelStatus& status);

    /// Add a link.  Validates: domains registered, source read channel exists,
    /// target write side exists, probe points non-empty, id unique,
    /// relaxation in (0,1], sub_iterations >= 1.
    bool add_link(const CouplingLink& link, exd::physics::ModelStatus& status);

    /// Explicit (staggered) exchange: applies every due link once.  Returns
    /// the number of links executed.  On the first hard failure processing
    /// stops, status.ok is false, and the count covers only completed links.
    int exchange(double t, exd::physics::ModelStatus& status);

    /// Implicit exchange: sub-iterates every due link until the max probe-point
    /// delta between consecutive sweeps is below `tolerance` or the link's
    /// sub_iterations are exhausted.  Returns the number of links executed.
    /// Non-converged links produce a warning; hard failures set status.ok=false.
    int exchange_implicit(double t, double tolerance, exd::physics::ModelStatus& status);

    /// Whether the most recent exchange_implicit converged for all due links
    /// (true when no links were due).
    bool last_implicit_converged() const { return last_implicit_converged_; }

    size_t link_count() const { return links_.size(); }
    size_t domain_count() const { return domains_.size(); }

private:
    struct LinkRuntime
    {
        CouplingLink link;
        double last_exchange = -std::numeric_limits<double>::infinity();
        // Per-probe relaxed value last written.  The manager has no read side
        // to a solver's current value, so the first exchange treats the
        // previous value as 0.0 and tracks from there.
        std::vector<double> last_written_scalar;
        std::vector<std::array<double, 3>> last_written_vector;
    };

    static bool due(const LinkRuntime& rt, double t);
    const DomainHandle* find_domain(const std::string& name) const;
    bool perform_exchange(LinkRuntime& rt, double t, exd::physics::ModelStatus& status);
    bool perform_implicit_exchange(LinkRuntime& rt, double t, double tolerance,
                                   exd::physics::ModelStatus& status,
                                   bool& converged_out);

    std::vector<DomainHandle> domains_;
    std::vector<LinkRuntime> links_;
    bool last_implicit_converged_ = true;
};

/// A solver domain registered with a coupled run: `step` advances the domain
/// by its own timestep and `handle` exposes its coupling channels.
struct CoupledDomainSpec
{
    std::string name;
    std::function<bool(double dt)> step; // advances the domain by dt (returns ok)
    double dt = 0.0;                     // per-domain timestep (> 0; multi-rate allowed)
    CouplingManager::DomainHandle handle; // channels / sinks for coupling
};

/// Macro-step driver over a set of coupled domains.
///
/// Staggered mode: each macro-step advances all domains (sub-cycling smaller
/// timesteps to cover the macro-step), then calls exchange().  Implicit mode:
/// advance, then exchange_implicit(t, tolerance).  The link fields inherit
/// the Config values when a link keeps its struct defaults (relaxation 1.0,
/// sub_iterations 1), so one Config drives the whole run unless a link
/// overrides.
class CoupledSimulation
{
public:
    struct Config
    {
        bool implicit = false;
        double relaxation = 1.0;
        double tolerance = 1e-8;
        int max_sub_iterations = 50;
    };

    struct RunReport
    {
        bool ok = false;
        exd::physics::ModelStatus status;
        double total_time = 0.0;      // scheduler macro-clock when the run stopped
        size_t steps = 0;             // macro-steps taken
        size_t total_exchanges = 0;   // links executed over the whole run
        bool converged = true;        // false when any implicit exchange failed to converge
    };

    CoupledSimulation() = default;
    explicit CoupledSimulation(const Config& config);
    void set_config(const Config& config);

    bool add_domain(const CoupledDomainSpec& spec, exd::physics::ModelStatus& status);
    bool add_link(const CouplingLink& link, exd::physics::ModelStatus& status);

    /// Run to `time_horizon`.  The macro-step is max(dt over domains); every
    /// domain is advanced by its own dt repeatedly until it has covered the
    /// macro-step (v1 sub-cycling), then the due links are exchanged.
    RunReport run(double time_horizon, exd::physics::ModelStatus& status);

private:
    Config config_;
    std::vector<CoupledDomainSpec> domains_;
    CouplingManager manager_;
};

} // namespace exd::physics::coupling