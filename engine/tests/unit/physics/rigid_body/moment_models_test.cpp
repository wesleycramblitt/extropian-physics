#include <exd/engine/physics/rigid_body/moment_model.hpp>
#include <exd/engine/physics/rigid_body/status.hpp>

#include <doctest/doctest.h>

using namespace exd::engine::physics::rigid_body;

TEST_CASE("Constant moment model: returns the configured torque independent of state")
{
    ConstantMomentConfig config;
    config.torque = 40.0;

    auto model = make_constant_moment(config);
    REQUIRE(model);
    CHECK(model->name() == "constant");

    RotationalState state;
    ModelStatus status;

    state.omega = 0.0;
    CHECK(model->moment(state, status) == doctest::Approx(40.0).epsilon(1e-12));

    state.omega = 100.0;
    CHECK(model->moment(state, status) == doctest::Approx(40.0).epsilon(1e-12));

    CHECK(status.ok);
    CHECK(status.error.empty());
}

TEST_CASE("Linear moment model: T = k*omega + offset at several speeds")
{
    LinearMomentConfig config;
    config.k = 2.0;
    config.offset = 3.0;

    auto model = make_linear_moment(config);
    REQUIRE(model);
    CHECK(model->name() == "linear");

    RotationalState state;
    ModelStatus status;

    state.omega = 0.0;
    CHECK(model->moment(state, status) == doctest::Approx(3.0).epsilon(1e-12));

    state.omega = 5.0;
    CHECK(model->moment(state, status) == doctest::Approx(13.0).epsilon(1e-12));

    state.omega = -2.0;
    CHECK(model->moment(state, status) == doctest::Approx(-1.0).epsilon(1e-12));

    CHECK(status.ok);
}

TEST_CASE("Curve moment model: piecewise-linear interpolation between table points")
{
    CurveMomentConfig config;
    config.omega_pts = {0.0, 5.0, 10.0};
    config.torque_pts = {2.0, 1.0, 0.5};

    auto model = make_curve_moment(config);
    REQUIRE(model);
    CHECK(model->name() == "curve");

    RotationalState state;
    ModelStatus status;

    state.omega = 2.5; // halfway on the segment {0,2.0} -> {5,1.0}
    CHECK(model->moment(state, status) == doctest::Approx(1.5).epsilon(1e-12));

    state.omega = 7.5; // halfway on the segment {5,1.0} -> {10,0.5}
    CHECK(model->moment(state, status) == doctest::Approx(0.75).epsilon(1e-12));

    CHECK(status.ok);
    CHECK(status.error.empty());
}

TEST_CASE("Curve moment model: clamps flat beyond the table ends")
{
    CurveMomentConfig config;
    config.omega_pts = {0.0, 5.0, 10.0};
    config.torque_pts = {2.0, 1.0, 0.5};

    auto model = make_curve_moment(config);
    REQUIRE(model);

    RotationalState state;
    ModelStatus status;

    state.omega = -3.0; // below the first point
    CHECK(model->moment(state, status) == doctest::Approx(2.0).epsilon(1e-12));

    state.omega = 42.0; // above the last point
    CHECK(model->moment(state, status) == doctest::Approx(0.5).epsilon(1e-12));

    CHECK(status.ok);
}

TEST_CASE("Curve moment model: empty config is a documented no-load model")
{
    CurveMomentConfig config; // both tables empty

    auto model = make_curve_moment(config);
    REQUIRE(model);

    RotationalState state;
    state.omega = 8.0;

    ModelStatus status;
    CHECK(model->moment(state, status) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(status.ok);
    CHECK(status.error.empty());
}

TEST_CASE("Curve moment model: non-increasing omega table is a config error")
{
    CurveMomentConfig config;
    config.omega_pts = {5.0, 0.0, 10.0}; // not strictly increasing
    config.torque_pts = {1.0, 2.0, 0.5};

    auto model = make_curve_moment(config);
    REQUIRE(model);

    RotationalState state;
    state.omega = 2.0;

    ModelStatus status;
    const double result = model->moment(state, status);

    CHECK(result == doctest::Approx(0.0).epsilon(1e-12));
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
}

TEST_CASE("Curve moment model: size-mismatched tables are a config error")
{
    CurveMomentConfig config;
    config.omega_pts = {0.0, 5.0, 10.0};
    config.torque_pts = {2.0, 1.0}; // one point short

    auto model = make_curve_moment(config);
    REQUIRE(model);

    RotationalState state;
    state.omega = 2.0;

    ModelStatus status;
    const double result = model->moment(state, status);

    CHECK(result == doctest::Approx(0.0).epsilon(1e-12));
    CHECK_FALSE(status.ok);
    CHECK_FALSE(status.error.empty());
}