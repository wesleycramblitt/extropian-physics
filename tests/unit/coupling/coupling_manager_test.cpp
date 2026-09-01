// Coupling-manager tests: explicit/implicit exchange, interval gating,
// under-relaxation, validation failures, and a coupled linear ODE sanity
// system (staggered vs implicit agreement plus multi-rate).

#include <exd/physics/coupling/coupling_manager.hpp>
#include <exd/physics/model_status.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <vector>

using namespace exd::physics;
using namespace exd::physics::coupling;

namespace
{

// -- Fakes --------------------------------------------------------

// Constant scalar field used as a fake source channel.
class ConstScalarField final : public IScalarField3D
{
public:
    explicit ConstScalarField(double value = 0.0)
        : value_(value)
    {
    }

    bool sample(const std::array<double, 3>&, double& value_out) const override
    {
        value_out = value_;
        return true;
    }

private:
    double value_;
};

// -- Coupled 2-ODE linear system ----------------------------------
//
//   dx/dt = -2x + y,   dy/dt = 2x - 3y
//
// The matrix [[-2,1],[2,-3]] has eigenvalues -1 and -4 with eigenvectors
// [1,1] and [1,-2].  With x(0) = x0, y(0) = y0:
//
//   x(t) = c1 e^(-t) + c2 e^(-4t)
//   y(t) = c1 e^(-t) - 2 c2 e^(-4t)
//   c1 = (2 x0 + y0) / 3,  c2 = (x0 - y0) / 3
struct OdeState
{
    double x = 2.0;
    double y = 1.0;
    double y_at_a = 0.0; // y value written into domain A (source B.y)
    double x_at_b = 0.0; // x value written into domain B (source A.x)
};

// Scalar channel whose value tracks a live state field.
class LiveScalarField final : public IScalarField3D
{
public:
    const double* value = nullptr;

    bool sample(const std::array<double, 3>&, double& value_out) const override
    {
        value_out = *value;
        return true;
    }
};

struct OdeFixture
{
    OdeState state;
    LiveScalarField x_field;
    LiveScalarField y_field;
    CoupledDomainSpec domain_a;
    CoupledDomainSpec domain_b;
    CouplingLink link_x_to_b;
    CouplingLink link_y_to_a;

    explicit OdeFixture(double dt_a, double dt_b)
    {
        x_field.value = &state.x;
        y_field.value = &state.y;

        domain_a.name = "A";
        domain_a.dt = dt_a;
        domain_a.step = [this](double dt) -> bool {
            state.x += dt * (-2.0 * state.x + state.y_at_a);
            return true;
        };
        domain_a.handle.name = "A";
        domain_a.handle.scalar_channel = [this](std::string_view name) -> const IScalarField3D* {
            return name == "x" ? &x_field : nullptr;
        };
        domain_a.handle.scalar_write =
            [this](std::string_view, const std::array<double, 3>&, double value) -> bool {
            state.y_at_a = value;
            return true;
        };

        domain_b.name = "B";
        domain_b.dt = dt_b;
        domain_b.step = [this](double dt) -> bool {
            state.y += dt * (2.0 * state.x_at_b - 3.0 * state.y);
            return true;
        };
        domain_b.handle.name = "B";
        domain_b.handle.scalar_channel = [this](std::string_view name) -> const IScalarField3D* {
            return name == "y" ? &y_field : nullptr;
        };
        domain_b.handle.scalar_write =
            [this](std::string_view, const std::array<double, 3>&, double value) -> bool {
            state.x_at_b = value;
            return true;
        };

        // A.x -> B probe (B receives the current x).
        link_x_to_b.id = "A.x->B";
        link_x_to_b.source_domain = "A";
        link_x_to_b.source_channel = "x";
        link_x_to_b.target_domain = "B";
        link_x_to_b.target_channel = "probe";
        link_x_to_b.probe_points = {{0.0, 0.0, 0.0}};

        // B.y -> A probe (A receives the current y).
        link_y_to_a.id = "B.y->A";
        link_y_to_a.source_domain = "B";
        link_y_to_a.source_channel = "y";
        link_y_to_a.target_domain = "A";
        link_y_to_a.target_channel = "probe";
        link_y_to_a.probe_points = {{0.0, 0.0, 0.0}};
    }

    /// Analytic solution at time t (derived above for x0=2, y0=1).
    static std::array<double, 2> analytic(double t)
    {
        const double x0 = 2.0;
        const double y0 = 1.0;
        const double c1 = (2.0 * x0 + y0) / 3.0;
        const double c2 = (x0 - y0) / 3.0;
        return {c1 * std::exp(-t) + c2 * std::exp(-4.0 * t),
                c1 * std::exp(-t) - 2.0 * c2 * std::exp(-4.0 * t)};
    }
};

// Register the two ODE domains and both links on a simulation.
void add_ode_system(CoupledSimulation& simulation, OdeFixture& fixture,
                    ModelStatus& status)
{
    REQUIRE(simulation.add_domain(fixture.domain_a, status));
    REQUIRE(status.ok);
    REQUIRE(simulation.add_domain(fixture.domain_b, status));
    REQUIRE(status.ok);
    REQUIRE(simulation.add_link(fixture.link_x_to_b, status));
    REQUIRE(status.ok);
    REQUIRE(simulation.add_link(fixture.link_y_to_a, status));
    REQUIRE(status.ok);
}

} // anonymous namespace

// -- Explicit exchange --------------------------------------------

TEST_CASE("explicit exchange applies sampled values")
{
    ConstScalarField source(42.0);
    std::vector<double> written;

    CouplingManager manager;
    ModelStatus status;

    CouplingManager::DomainHandle src;
    src.name = "src";
    src.scalar_channel = [&source](std::string_view name) -> const IScalarField3D* {
        return name == "field" ? &source : nullptr;
    };

    CouplingManager::DomainHandle dst;
    dst.name = "dst";
    dst.scalar_write = [&written](std::string_view, const std::array<double, 3>&, double value) -> bool {
        written.push_back(value);
        return true;
    };

    REQUIRE(manager.register_domain(src, status));
    REQUIRE(status.ok);
    REQUIRE(manager.register_domain(dst, status));
    REQUIRE(status.ok);

    CouplingLink link;
    link.id = "constant";
    link.source_domain = "src";
    link.source_channel = "field";
    link.target_domain = "dst";
    link.target_channel = "sink";
    link.interval = 0.0;
    link.relaxation = 1.0;
    link.probe_points = {{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}};

    REQUIRE(manager.add_link(link, status));
    REQUIRE(status.ok);

    const int executed = manager.exchange(0.0, status);
    REQUIRE(status.ok);
    CHECK(executed == 1);
    REQUIRE(written.size() == 2);
    CHECK(written[0] == doctest::Approx(42.0));
    CHECK(written[1] == doctest::Approx(42.0));
}

TEST_CASE("vector link transfers all three components")
{
    class ConstVectorField final : public IVectorField3D
    {
    public:
        std::array<double, 3> value = {0.0, 0.0, 0.0};

        bool sample(const std::array<double, 3>&,
                    std::array<double, 3>& value_out) const override
        {
            value_out = value;
            return true;
        }
    };

    ConstVectorField source;
    source.value = {1.0, 2.0, 3.0};
    std::vector<std::array<double, 3>> written;

    CouplingManager manager;
    ModelStatus status;

    CouplingManager::DomainHandle src;
    src.name = "vsrc";
    src.vector_channel = [&source](std::string_view name) -> const IVectorField3D* {
        return name == "E" ? &source : nullptr;
    };

    CouplingManager::DomainHandle dst;
    dst.name = "vdst";
    dst.vector_write = [&written](std::string_view, const std::array<double, 3>&,
                                  const std::array<double, 3>& value) -> bool {
        written.push_back(value);
        return true;
    };

    REQUIRE(manager.register_domain(src, status));
    REQUIRE(manager.register_domain(dst, status));

    CouplingLink link;
    link.id = "vector";
    link.source_domain = "vsrc";
    link.source_channel = "E";
    link.target_domain = "vdst";
    link.target_channel = "sink";
    link.vector_field = true;
    link.probe_points = {{0.0, 0.0, 0.0}};

    REQUIRE(manager.add_link(link, status));
    REQUIRE(status.ok);

    const int executed = manager.exchange(0.0, status);
    REQUIRE(status.ok);
    CHECK(executed == 1);
    REQUIRE(written.size() == 1);
    CHECK(written[0][0] == doctest::Approx(1.0));
    CHECK(written[0][1] == doctest::Approx(2.0));
    CHECK(written[0][2] == doctest::Approx(3.0));
}

// -- Interval gating ----------------------------------------------

TEST_CASE("interval gating")
{
    ConstScalarField source(7.0);
    std::vector<double> written;

    CouplingManager manager;
    ModelStatus status;

    CouplingManager::DomainHandle src;
    src.name = "src";
    src.scalar_channel = [&source](std::string_view) -> const IScalarField3D* {
        return &source;
    };

    CouplingManager::DomainHandle dst;
    dst.name = "dst";
    dst.scalar_write = [&written](std::string_view, const std::array<double, 3>&, double value) -> bool {
        written.push_back(value);
        return true;
    };

    REQUIRE(manager.register_domain(src, status));
    REQUIRE(manager.register_domain(dst, status));

    CouplingLink link;
    link.id = "gated";
    link.source_domain = "src";
    link.source_channel = "field";
    link.target_domain = "dst";
    link.target_channel = "sink";
    link.interval = 0.5;
    link.probe_points = {{0.0, 0.0, 0.0}};

    REQUIRE(manager.add_link(link, status));
    REQUIRE(status.ok);

    CHECK(manager.exchange(0.0, status) == 1); // due (first exchange)
    REQUIRE(status.ok);
    CHECK(manager.exchange(0.2, status) == 0); // 0.2 < 0.5
    REQUIRE(status.ok);
    CHECK(manager.exchange(0.5, status) == 1); // 0.5 - 0.0 >= 0.5
    REQUIRE(status.ok);

    CHECK(written.size() == 2);
}

// -- Under-relaxation ---------------------------------------------

TEST_CASE("relaxation accumulates per probe")
{
    ConstScalarField source(100.0);
    std::vector<double> written;

    CouplingManager manager;
    ModelStatus status;

    CouplingManager::DomainHandle src;
    src.name = "src";
    src.scalar_channel = [&source](std::string_view) -> const IScalarField3D* {
        return &source;
    };

    CouplingManager::DomainHandle dst;
    dst.name = "dst";
    dst.scalar_write = [&written](std::string_view, const std::array<double, 3>&, double value) -> bool {
        written.push_back(value);
        return true;
    };

    REQUIRE(manager.register_domain(src, status));
    REQUIRE(manager.register_domain(dst, status));

    CouplingLink link;
    link.id = "relaxed";
    link.source_domain = "src";
    link.source_channel = "field";
    link.target_domain = "dst";
    link.target_channel = "sink";
    link.interval = 0.0;
    link.relaxation = 0.25;
    link.probe_points = {{0.0, 0.0, 0.0}};

    REQUIRE(manager.add_link(link, status));
    REQUIRE(status.ok);

    const int first = manager.exchange(0.0, status);
    REQUIRE(status.ok);
    const int second = manager.exchange(0.1, status);
    REQUIRE(status.ok);

    CHECK(first == 1);
    CHECK(second == 1);
    REQUIRE(written.size() == 2);
    // relaxed = (1 - w) * prev + w * sampled, starting from the manager's
    // per-probe baseline (0.0 first exchange).
    CHECK(written[0] == doctest::Approx(25.0));
    CHECK(written[1] == doctest::Approx(43.75));
}

// -- Validation failures ------------------------------------------

TEST_CASE("validation failures")
{
    ConstScalarField source(1.0);
    std::vector<double> written;

    CouplingManager manager;
    ModelStatus status;

    CouplingManager::DomainHandle src;
    src.name = "src";
    src.scalar_channel = [&source](std::string_view name) -> const IScalarField3D* {
        return name == "field" ? &source : nullptr;
    };

    CouplingManager::DomainHandle dst;
    dst.name = "dst";
    dst.scalar_write = [&written](std::string_view, const std::array<double, 3>&, double value) -> bool {
        written.push_back(value);
        return true;
    };

    REQUIRE(manager.register_domain(src, status));
    REQUIRE(status.ok);
    REQUIRE(manager.register_domain(dst, status));
    REQUIRE(status.ok);

    auto make_link = []() -> CouplingLink {
        CouplingLink link;
        link.id = "default";
        link.source_domain = "src";
        link.source_channel = "field";
        link.target_domain = "dst";
        link.target_channel = "sink";
        link.probe_points = {{0.0, 0.0, 0.0}};
        return link;
    };

    SUBCASE("unknown source domain")
    {
        auto link = make_link();
        link.source_domain = "ghost";
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("unknown target domain")
    {
        auto link = make_link();
        link.target_domain = "ghost";
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("unknown source channel")
    {
        auto link = make_link();
        link.source_channel = "missing";
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("missing write side")
    {
        CouplingManager::DomainHandle sinkless;
        sinkless.name = "sinkless";
        REQUIRE(manager.register_domain(sinkless, status));
        REQUIRE(status.ok);

        auto link = make_link();
        link.target_domain = "sinkless";
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("empty probe set")
    {
        auto link = make_link();
        link.probe_points.clear();
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("duplicate id")
    {
        auto link = make_link();
        REQUIRE(manager.add_link(link, status));
        REQUIRE(status.ok);
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("relaxation of zero")
    {
        auto link = make_link();
        link.relaxation = 0.0;
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("relaxation above one")
    {
        auto link = make_link();
        link.relaxation = 1.5;
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }

    SUBCASE("sub_iterations of zero")
    {
        auto link = make_link();
        link.sub_iterations = 0;
        CHECK_FALSE(manager.add_link(link, status));
        CHECK_FALSE(status.ok);
        CHECK_FALSE(status.error.empty());
    }
}

// -- Coupled linear ODE system ------------------------------------

TEST_CASE("staggered vs implicit agreement on a linear system")
{
    SUBCASE("staggered mode matches the analytic solution within 1%")
    {
        OdeFixture fixture(1e-4, 1e-4);
        CoupledSimulation::Config config;
        config.implicit = false;

        CoupledSimulation simulation(config);
        ModelStatus status;
        add_ode_system(simulation, fixture, status);

        const auto report = simulation.run(2.0, status);
        REQUIRE(report.ok);
        const auto exact = OdeFixture::analytic(2.0);
        CHECK(fixture.state.x == doctest::Approx(exact[0]).epsilon(0.01));
        CHECK(fixture.state.y == doctest::Approx(exact[1]).epsilon(0.01));
    }

    SUBCASE("implicit mode converges and matches the analytic within 1e-3")
    {
        OdeFixture fixture(1e-4, 1e-4);
        CoupledSimulation::Config config;
        config.implicit = true;
        config.relaxation = 0.5;
        config.tolerance = 1e-8;

        CoupledSimulation simulation(config);
        ModelStatus status;
        add_ode_system(simulation, fixture, status);

        const auto report = simulation.run(2.0, status);
        REQUIRE(report.ok);
        CHECK(report.converged);
        const auto exact = OdeFixture::analytic(2.0);
        CHECK(fixture.state.x == doctest::Approx(exact[0]).epsilon(1e-3));
        CHECK(fixture.state.y == doctest::Approx(exact[1]).epsilon(1e-3));
    }
}

TEST_CASE("multi-rate: slower domain still converges within 2%")
{
    // Domain B steps with twice the timestep of domain A.
    OdeFixture fixture(1e-4, 2e-4);
    CoupledSimulation::Config config;
    config.implicit = false;

    CoupledSimulation simulation(config);
    ModelStatus status;
    add_ode_system(simulation, fixture, status);

    const auto report = simulation.run(2.0, status);
    REQUIRE(report.ok);
    const auto exact = OdeFixture::analytic(2.0);
    CHECK(fixture.state.x == doctest::Approx(exact[0]).epsilon(0.02));
    CHECK(fixture.state.y == doctest::Approx(exact[1]).epsilon(0.02));
}