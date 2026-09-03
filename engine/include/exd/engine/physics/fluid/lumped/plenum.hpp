#pragma once

// ---------------------------------------------------------------------------
// Lumped plenum (Greitzer-style surge cell) -- isentropic plenum mapping.
//
// State: (p_plenum, mdot_duct)
//   dmdot_dt = (compressor(mdot) - (p - p_ambient)) / I    [I = L/A inertance]
//   dp_dt     = (a^2 / volume) * (mdot - throttle(p))
//
// The plenum sound speed uses the ISENTROPIC mapping from the ambient
// reference state:
//   T_p = T_ambient * (p / p_ambient) ^ ((gamma - 1) / gamma)
//
// Assumption (documented): the lumped cell stores compression work as
// pressure alone and the evolution is fast enough that heat transfer and
// entropy generation inside the cell are neglected -- the compressor's own
// losses live in the compressor characteristic. This is the classical
// Greitzer surge-cell approximation and is appropriate when the cell volume
// is large compared with the duct.
//
// gamma and R are taken from the IEos (air defaults gamma = 1.4,
// R = 287.05 J/(kg.K)). Errors are reported through ModelStatus and never
// thrown.
// ---------------------------------------------------------------------------

#include <exd/engine/core/model_status.hpp>
#include <exd/engine/numerics/integrators.hpp>
#include <exd/engine/physics/thermo/eos.hpp>

#include <functional>
#include <string>
#include <vector>

namespace exd::engine::physics::fluid::lumped {

struct PlenumModelConfig
{
    double volume = 0.5;         // plenum volume (m^3, > 0)
    double duct_area = 0.05;     // duct cross-section (m^2, > 0)
    double duct_length = 2.0;    // duct length (m, > 0)
    double p_ambient = 101325.0; // ambient/reference pressure (Pa, > 0)
    double T_ambient = 288.15;   // ambient/reference temperature (K, > 0)
};

/// Validate the plenum config. Fatal problems return false with `error`;
/// non-fatal observations live in `warnings`.
bool validate_plenum_config(const PlenumModelConfig& config,
                            std::string& error,
                            std::vector<std::string>& warnings);

struct PlenumState
{
    double p_plenum = 101325.0; // plenum pressure (Pa), > 0
    double mdot_duct = 0.0;     // duct mass flow (kg/s), finite
};

// Compressor: pressure RISE (Pa) vs duct flow (kg/s).
// Throttle: flow (kg/s) vs plenum pressure (Pa).
using PlenumCharacteristic = std::function<double(double)>;

struct PlenumDerivative
{
    bool ok = false;
    exd::engine::core::ModelStatus status;

    double dp_dt = 0.0;            // plenum pressure rate (Pa/s)
    double dmdot_dt = 0.0;         // duct flow rate ((kg/s)/s)
    double speed_of_sound_sq = 0.0;// isentropic sound speed^2 in the plenum (m^2/s^2)
    double inertia = 0.0;          // duct inertance I = L/A (1/m)
};

/// Evaluate the lumped-plenum state derivative at `state`.
/// `compressor` maps duct flow to pressure RISE over ambient; `throttle`
/// maps plenum pressure to flow. Errors are mirrored into both the return
/// value and `status`.
PlenumDerivative plenum_derivative(const PlenumModelConfig& config,
                                   const PlenumState& state,
                                   const PlenumCharacteristic& compressor,
                                   const PlenumCharacteristic& throttle,
                                   const exd::engine::physics::thermo::IEos& eos,
                                   exd::engine::core::ModelStatus& status);

/// Advance the state by `dt` with `integration` (default RK4). The state
/// vector is {p_plenum, mdot_duct}. On failure the input state is returned
/// unchanged with status.ok = false.
PlenumState step_plenum(double dt,
                        const PlenumModelConfig& config,
                        const PlenumState& state,
                        const PlenumCharacteristic& compressor,
                        const PlenumCharacteristic& throttle,
                        const exd::engine::physics::thermo::IEos& eos,
                        const exd::engine::numerics::IntegratorConfig& integration,
                        exd::engine::core::ModelStatus& status);

} // namespace exd::engine::physics::fluid::lumped