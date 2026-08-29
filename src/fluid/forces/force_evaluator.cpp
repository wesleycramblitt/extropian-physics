// force_evaluator.cpp
// Dispatcher factory: constructs the concrete IForceEvaluator implementation
// selected by ForceEvaluatorParams::type. Concrete evaluator factories are
// defined in their own translation units (repo factory pattern).

#include <exd/physics/fluid/forces/force_evaluator.hpp>

#include <memory>

namespace exd::physics::fluid::forces
{

// Concrete evaluator factories (defined in their own translation units).
std::unique_ptr<IForceEvaluator> make_pressure_integration_evaluator(
    const PressureIntegrationConfig& config);
std::unique_ptr<IForceEvaluator> make_momentum_balance_evaluator(
    const MomentumBalanceConfig& config, const PolarDatabase* polars);
std::unique_ptr<IForceEvaluator> make_table_lookup_evaluator(
    const TableLookupConfig& config, const PolarDatabase* polars);

std::unique_ptr<IForceEvaluator> make_force_evaluator(const ForceEvaluatorParams& params,
                                                      mechanics::ModelStatus& /*status*/)
{
    switch (params.type)
    {
        case ForceEvaluatorType::PressureIntegration:
            return make_pressure_integration_evaluator(params.pressure);
        case ForceEvaluatorType::TableLookup:
            return make_table_lookup_evaluator(params.table, params.polars);
        case ForceEvaluatorType::MomentumBalance:
        default:
            return make_momentum_balance_evaluator(params.momentum, params.polars);
    }
}

} // namespace exd::physics::fluid::forces