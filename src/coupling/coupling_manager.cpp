// Coupling manager (Phase H-lite): the real coupling exchange.  A manager
// owns a set of registered domains (channel callbacks) and directional links.
// Each link reads a source channel at its probe points and writes the
// (possibly relaxed) value into a target domain's sink.
//
// Operator-splitting coupling with under-relaxation: the manager tracks the
// last value it wrote per probe and writes
//     relaxed = (1 - w) * last_written + w * sampled
// so the explicit exchange can be under-relaxed and the implicit exchange is
// a sub-iterated fixed point over those relaxed writes.

#include <exd/physics/coupling/coupling_manager.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>

namespace exd::physics::coupling
{

namespace
{

// Interval-gating tolerance so boundary queries like t = last + interval are
// treated as due despite floating-point subtraction.
constexpr double kIntervalEpsilon = 1e-12;

} // anonymous namespace

// -- CouplingManager ----------------------------------------------

bool CouplingManager::register_domain(const DomainHandle& domain,
                                      exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    if (domain.name.empty())
    {
        status.ok = false;
        status.error = "register_domain: domain name must not be empty";
        return false;
    }
    if (find_domain(domain.name) != nullptr)
    {
        status.ok = false;
        status.error = "register_domain: domain '" + domain.name + "' already registered";
        return false;
    }

    domains_.push_back(domain);
    return true;
}

bool CouplingManager::add_link(const CouplingLink& link, exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    auto fail = [&status](const std::string& message)
    {
        status.ok = false;
        status.error = message;
        return false;
    };

    if (link.id.empty())
    {
        return fail("add_link: link id must not be empty");
    }
    for (const auto& rt : links_)
    {
        if (rt.link.id == link.id)
        {
            return fail("add_link: duplicate link id '" + link.id + "'");
        }
    }

    const DomainHandle* source = find_domain(link.source_domain);
    if (source == nullptr)
    {
        return fail("add_link: source domain '" + link.source_domain + "' is not registered");
    }
    const DomainHandle* target = find_domain(link.target_domain);
    if (target == nullptr)
    {
        return fail("add_link: target domain '" + link.target_domain + "' is not registered");
    }

    if (link.probe_points.empty())
    {
        return fail("add_link: link '" + link.id + "' has no probe points");
    }
    if (!(link.relaxation > 0.0) || link.relaxation > 1.0)
    {
        return fail("add_link: link '" + link.id + "' relaxation must be in (0, 1]");
    }
    if (link.sub_iterations < 1)
    {
        return fail("add_link: link '" + link.id + "' sub_iterations must be >= 1");
    }

    if (!link.vector_field)
    {
        if (!source->scalar_channel)
        {
            return fail("add_link: source domain '" + link.source_domain +
                        "' has no scalar channel lookup");
        }
        if (source->scalar_channel(link.source_channel) == nullptr)
        {
            return fail("add_link: source scalar channel '" + link.source_channel +
                        "' not found on domain '" + link.source_domain + "'");
        }
        if (!target->scalar_write)
        {
            return fail("add_link: target domain '" + link.target_domain +
                        "' has no scalar sink");
        }
    }
    else
    {
        if (!source->vector_channel)
        {
            return fail("add_link: source domain '" + link.source_domain +
                        "' has no vector channel lookup");
        }
        if (source->vector_channel(link.source_channel) == nullptr)
        {
            return fail("add_link: source vector channel '" + link.source_channel +
                        "' not found on domain '" + link.source_domain + "'");
        }
        if (!target->vector_write)
        {
            return fail("add_link: target domain '" + link.target_domain +
                        "' has no vector sink");
        }
    }

    LinkRuntime rt;
    rt.link = link;
    rt.last_written_scalar.assign(link.probe_points.size(), 0.0);
    rt.last_written_vector.assign(link.probe_points.size(),
                                  std::array<double, 3>{0.0, 0.0, 0.0});
    links_.push_back(std::move(rt));
    return true;
}

bool CouplingManager::due(const LinkRuntime& rt, double t)
{
    if (rt.link.interval <= 0.0)
    {
        return true; // every step
    }
    return (t - rt.last_exchange) >= rt.link.interval - kIntervalEpsilon;
}

const CouplingManager::DomainHandle* CouplingManager::find_domain(const std::string& name) const
{
    for (const auto& domain : domains_)
    {
        if (domain.name == name)
        {
            return &domain;
        }
    }
    return nullptr;
}

bool CouplingManager::gather_exchange(LinkRuntime& rt,
                                        exd::physics::ModelStatus& status)
{
    const CouplingLink& link = rt.link;
    const DomainHandle* source = find_domain(link.source_domain);
    const DomainHandle* target = find_domain(link.target_domain);
    if (source == nullptr || target == nullptr)
    {
        status.ok = false;
        status.error = "exchange: link '" + link.id +
                       "' references a domain that is no longer registered";
        return false;
    }

    rt.pending_sampled.assign(link.probe_points.size(), 0.0);
    rt.pending_previous.assign(link.probe_points.size(), 0.0);

    if (!link.vector_field)
    {
        if (!source->scalar_channel || !target->scalar_write)
        {
            status.ok = false;
            status.error = "exchange: link '" + link.id +
                           "' scalar read/write side is no longer available";
            return false;
        }
        const IScalarField3D* channel = source->scalar_channel(link.source_channel);
        if (channel == nullptr)
        {
            status.ok = false;
            status.error = "exchange: link '" + link.id + "' source scalar channel '" +
                           link.source_channel + "' not found";
            return false;
        }
        for (std::size_t i = 0; i < link.probe_points.size(); ++i)
        {
            if (!channel->sample(link.probe_points[i], rt.pending_sampled[i]))
            {
                status.warnings.push_back("exchange: link '" + link.id +
                                          "' source sample out of bounds at probe " +
                                          std::to_string(i) + "; probe skipped");
                rt.pending_sampled[i] = 0.0;
            }
            rt.pending_previous[i] = rt.last_written_scalar[i];
            if (target->scalar_read)
            {
                double current = rt.pending_previous[i];
                if (!target->scalar_read(link.target_channel, link.probe_points[i], current))
                {
                    status.warnings.push_back("exchange: link '" + link.id +
                                              "' target read-back failed at probe " +
                                              std::to_string(i) + "; using last written");
                }
                else
                {
                    rt.pending_previous[i] = current;
                }
            }
        }
    }
    else
    {
        if (!source->vector_channel || !target->vector_write)
        {
            status.ok = false;
            status.error = "exchange: link '" + link.id +
                           "' vector read/write side is no longer available";
            return false;
        }
        const IVectorField3D* channel = source->vector_channel(link.source_channel);
        if (channel == nullptr)
        {
            status.ok = false;
            status.error = "exchange: link '" + link.id + "' source vector channel '" +
                           link.source_channel + "' not found";
            return false;
        }
        rt.pending_sampled_vec.assign(link.probe_points.size(),
                                      std::array<double, 3>{0.0, 0.0, 0.0});
        for (std::size_t i = 0; i < link.probe_points.size(); ++i)
        {
            if (!channel->sample(link.probe_points[i], rt.pending_sampled_vec[i]))
            {
                status.warnings.push_back("exchange: link '" + link.id +
                                          "' source sample out of bounds at probe " +
                                          std::to_string(i) + "; probe skipped");
                rt.pending_sampled_vec[i] = {0.0, 0.0, 0.0};
            }
        }
    }
    return true;
}

bool CouplingManager::apply_exchange(LinkRuntime& rt,
                                     exd::physics::ModelStatus& status)
{
    const CouplingLink& link = rt.link;
    const DomainHandle* target = find_domain(link.target_domain);
    if (target == nullptr)
    {
        status.ok = false;
        status.error = "exchange: link '" + link.id +
                       "' references a domain that is no longer registered";
        return false;
    }
    const double w = link.relaxation;

    if (!link.vector_field)
    {
        for (std::size_t i = 0; i < link.probe_points.size(); ++i)
        {
            const double relaxed =
                (1.0 - w) * rt.pending_previous[i] + w * rt.pending_sampled[i];
            if (!target->scalar_write(link.target_channel, link.probe_points[i], relaxed))
            {
                status.ok = false;
                status.error = "exchange: link '" + link.id + "' target scalar write failed "
                               "at probe " + std::to_string(i);
                return false;
            }
            rt.last_written_scalar[i] = relaxed;
        }
    }
    else
    {
        for (std::size_t i = 0; i < link.probe_points.size(); ++i)
        {
            std::array<double, 3> relaxed{};
            for (int c = 0; c < 3; ++c)
            {
                const double prev = rt.last_written_vector[i][c];
                relaxed[c] = (1.0 - w) * prev + w * rt.pending_sampled_vec[i][c];
            }
            if (!target->vector_write(link.target_channel, link.probe_points[i], relaxed))
            {
                status.ok = false;
                status.error = "exchange: link '" + link.id + "' target vector write failed "
                               "at probe " + std::to_string(i);
                return false;
            }
            rt.last_written_vector[i] = relaxed;
        }
    }
    return true;
}

bool CouplingManager::perform_implicit_exchange(LinkRuntime& rt, double t, double tolerance,
                                                exd::physics::ModelStatus& status,
                                                bool& converged_out)
{
    const CouplingLink& link = rt.link;
    const DomainHandle* source = find_domain(link.source_domain);
    const DomainHandle* target = find_domain(link.target_domain);
    if (source == nullptr || target == nullptr)
    {
        status.ok = false;
        status.error = "exchange_implicit: link '" + link.id +
                       "' references a domain that is no longer registered";
        return false;
    }

    const double w = link.relaxation;
    const int max_iterations = std::max(1, link.sub_iterations);
    bool converged = false;

    if (!link.vector_field)
    {
        if (!source->scalar_channel || !target->scalar_write)
        {
            status.ok = false;
            status.error = "exchange_implicit: link '" + link.id +
                           "' scalar read/write side is no longer available";
            return false;
        }
        const IScalarField3D* channel = source->scalar_channel(link.source_channel);
        if (channel == nullptr)
        {
            status.ok = false;
            status.error = "exchange_implicit: link '" + link.id +
                           "' source scalar channel '" + link.source_channel + "' not found";
            return false;
        }

        for (int iter = 0; iter < max_iterations; ++iter)
        {
            double max_delta = 0.0;
            std::size_t skipped = 0;
            for (std::size_t i = 0; i < link.probe_points.size(); ++i)
            {
                double sampled = 0.0;
                if (!channel->sample(link.probe_points[i], sampled))
                {
                    ++skipped;
                    continue;
                }
                double previous = rt.last_written_scalar[i];
                if (target->scalar_read)
                {
                    double current = previous;
                    if (!target->scalar_read(link.target_channel, link.probe_points[i], current))
                    {
                        status.warnings.push_back("exchange_implicit: link '" + link.id +
                                                  "' target read-back failed at probe " +
                                                  std::to_string(i) + "; using last written");
                    }
                    else
                    {
                        previous = current;
                    }
                }
                const double relaxed = (1.0 - w) * previous + w * sampled;
                if (!target->scalar_write(link.target_channel, link.probe_points[i], relaxed))
                {
                    status.ok = false;
                    status.error = "exchange_implicit: link '" + link.id +
                                   "' target scalar write failed at probe " + std::to_string(i);
                    return false;
                }
                rt.last_written_scalar[i] = relaxed;
                max_delta = std::max(max_delta, std::fabs(relaxed - previous));
            }
            if (skipped > 0 && iter == 0)
            {
                status.warnings.push_back("exchange_implicit: link '" + link.id + "' " +
                                          std::to_string(skipped) +
                                          " probe(s) out of bounds; skipped");
            }
            if (max_delta < tolerance)
            {
                converged = true;
                break;
            }
        }
    }
    else
    {
        if (!source->vector_channel || !target->vector_write)
        {
            status.ok = false;
            status.error = "exchange_implicit: link '" + link.id +
                           "' vector read/write side is no longer available";
            return false;
        }
        const IVectorField3D* channel = source->vector_channel(link.source_channel);
        if (channel == nullptr)
        {
            status.ok = false;
            status.error = "exchange_implicit: link '" + link.id +
                           "' source vector channel '" + link.source_channel + "' not found";
            return false;
        }

        for (int iter = 0; iter < max_iterations; ++iter)
        {
            double max_delta = 0.0;
            std::size_t skipped = 0;
            for (std::size_t i = 0; i < link.probe_points.size(); ++i)
            {
                std::array<double, 3> sampled{0.0, 0.0, 0.0};
                if (!channel->sample(link.probe_points[i], sampled))
                {
                    ++skipped;
                    continue;
                }
                const std::array<double, 3>& previous = rt.last_written_vector[i];
                std::array<double, 3> relaxed;
                for (int c = 0; c < 3; ++c)
                {
                    relaxed[c] = (1.0 - w) * previous[c] + w * sampled[c];
                }
                if (!target->vector_write(link.target_channel, link.probe_points[i], relaxed))
                {
                    status.ok = false;
                    status.error = "exchange_implicit: link '" + link.id +
                                   "' target vector write failed at probe " + std::to_string(i);
                    return false;
                }
                rt.last_written_vector[i] = relaxed;
                max_delta = std::max(max_delta,
                                     std::fabs(relaxed[0] - previous[0]) +
                                         std::fabs(relaxed[1] - previous[1]) +
                                         std::fabs(relaxed[2] - previous[2]));
            }
            if (skipped > 0 && iter == 0)
            {
                status.warnings.push_back("exchange_implicit: link '" + link.id + "' " +
                                          std::to_string(skipped) +
                                          " probe(s) out of bounds; skipped");
            }
            if (max_delta < tolerance)
            {
                converged = true;
                break;
            }
        }
    }

    if (!converged)
    {
        status.warnings.push_back("exchange_implicit: link '" + link.id +
                                  "' did not converge within " + std::to_string(max_iterations) +
                                  " sub-iteration(s)");
    }
    converged_out = converged;

    rt.last_exchange = t;
    return true;
}

int CouplingManager::exchange(double t, exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    // JACOBI ORDER (W11): gather every due link's samples + read-backs
    // against the pre-exchange state, THEN apply all writes.  The two-slab
    // CHT acceptance test pins the symmetric interface fixed point; a
    // sequential per-link order would fold the write order into the result.
    int executed = 0;
    for (auto& rt : links_)
    {
        if (!due(rt, t))
        {
            continue;
        }
        if (!gather_exchange(rt, status))
        {
            status.ok = false;
            return executed;
        }
        ++executed;
    }
    for (auto& rt : links_)
    {
        if (!due(rt, t))
        {
            continue;
        }
        if (!apply_exchange(rt, status))
        {
            status.ok = false;
            return executed;
        }
        rt.last_exchange = t;
    }
    return executed;
}
int CouplingManager::exchange_implicit(double t, double tolerance,
                                       exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    int executed = 0;
    bool all_converged = true;
    for (auto& rt : links_)
    {
        if (!due(rt, t))
        {
            continue;
        }
        bool converged = true;
        if (!perform_implicit_exchange(rt, t, tolerance, status, converged))
        {
            status.ok = false;
            last_implicit_converged_ = false;
            return executed;
        }
        all_converged = all_converged && converged;
        ++executed;
    }
    last_implicit_converged_ = all_converged;
    return executed;
}

// -- CoupledSimulation --------------------------------------------

CoupledSimulation::CoupledSimulation(const Config& config)
    : config_(config)
{
}

void CoupledSimulation::set_config(const Config& config)
{
    config_ = config;
}

bool CoupledSimulation::add_domain(const CoupledDomainSpec& spec,
                                   exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    if (spec.name.empty())
    {
        status.ok = false;
        status.error = "CoupledSimulation::add_domain: domain name must not be empty";
        return false;
    }
    if (spec.dt <= 0.0)
    {
        status.ok = false;
        status.error = "CoupledSimulation::add_domain: domain '" + spec.name +
                       "' dt must be > 0";
        return false;
    }
    if (!spec.step)
    {
        status.ok = false;
        status.error = "CoupledSimulation::add_domain: domain '" + spec.name +
                       "' step must be callable";
        return false;
    }
    for (const auto& existing : domains_)
    {
        if (existing.name == spec.name)
        {
            status.ok = false;
            status.error = "CoupledSimulation::add_domain: duplicate domain '" +
                           spec.name + "'";
            return false;
        }
    }

    // Links resolve domains by name, so the registered handle must carry the
    // same name the domain is added under.
    if (!spec.handle.name.empty() && spec.handle.name != spec.name)
    {
        status.ok = false;
        status.error = "CoupledSimulation::add_domain: handle name '" +
                       spec.handle.name + "' does not match domain name '" + spec.name + "'";
        return false;
    }

    CoupledDomainSpec stored = spec;
    stored.handle.name = spec.name;
    if (!manager_.register_domain(stored.handle, status))
    {
        return false;
    }

    domains_.push_back(std::move(stored));
    return true;
}

bool CoupledSimulation::add_link(const CouplingLink& link, exd::physics::ModelStatus& status)
{
    status = exd::physics::ModelStatus{};

    // Config supplies the coupling defaults; a link keeps its own override
    // when it differs from the CouplingLink struct defaults.
    CouplingLink effective = link;
    if (effective.relaxation == 1.0)
    {
        effective.relaxation = config_.relaxation;
    }
    if (effective.sub_iterations == 1)
    {
        effective.sub_iterations = config_.max_sub_iterations;
    }

    return manager_.add_link(effective, status);
}

CoupledSimulation::RunReport CoupledSimulation::run(double time_horizon,
                                                    exd::physics::ModelStatus& status)
{
    RunReport report;
    status = exd::physics::ModelStatus{};

    if (domains_.empty() || manager_.domain_count() == 0)
    {
        status.ok = false;
        status.error = "CoupledSimulation::run: no domains registered";
        report.status = status;
        return report;
    }
    if (time_horizon < 0.0)
    {
        status.ok = false;
        status.error = "CoupledSimulation::run: time_horizon must be >= 0";
        report.status = status;
        return report;
    }

    double macro_dt = 0.0;
    for (const auto& domain : domains_)
    {
        if (!(domain.dt > 0.0))
        {
            status.ok = false;
            status.error = "CoupledSimulation::run: domain '" + domain.name +
                           "' dt must be > 0";
            report.status = status;
            return report;
        }
        macro_dt = std::max(macro_dt, domain.dt);
    }

    // v1 scheduler: the macro-step is max(dt over domains).  Every macro-step
    // advances each domain by repeated calls of the domain's own dt until it
    // has covered the macro-step (sub-cycling keeps multi-rate domains on the
    // same clock; exact for integer dt ratios).  The due links are then
    // exchanged once at the end of the macro-step.
    std::vector<double> pending(domains_.size(), 0.0);

    double t = 0.0;
    const size_t macro_plan = static_cast<size_t>(
        std::ceil(time_horizon / macro_dt));
    size_t steps = 0;
    size_t total_exchanges = 0;
    bool converged = true;

    for (size_t m = 0; m < macro_plan; ++m)
    {
        const double step = std::min(macro_dt, time_horizon - t);

        for (std::size_t i = 0; i < domains_.size(); ++i)
        {
            pending[i] += step;
            while (pending[i] >= domains_[i].dt - 1e-12 * domains_[i].dt)
            {
                if (!domains_[i].step(domains_[i].dt))
                {
                    status.ok = false;
                    status.error = "CoupledSimulation::run: domain '" +
                                   domains_[i].name + "' step failed";
                    report.status = status;
                    return report;
                }
                pending[i] -= domains_[i].dt;
            }
        }

        const double exchange_t = t + step;
        exd::physics::ModelStatus exchange_status;
        int executed = 0;
        if (config_.implicit)
        {
            executed = manager_.exchange_implicit(exchange_t, config_.tolerance,
                                                  exchange_status);
            if (!exchange_status.ok)
            {
                status = exchange_status;
                report.status = status;
                return report;
            }
            if (!manager_.last_implicit_converged())
            {
                converged = false;
            }
        }
        else
        {
            executed = manager_.exchange(exchange_t, exchange_status);
            if (!exchange_status.ok)
            {
                status = exchange_status;
                report.status = status;
                return report;
            }
        }
        total_exchanges += static_cast<size_t>(executed);

        t += step;
        ++steps;
    }

    report.ok = true;
    report.status = status;
    report.total_time = t;
    report.steps = steps;
    report.total_exchanges = total_exchanges;
    report.converged = converged;
    return report;
}

} // namespace exd::physics::coupling