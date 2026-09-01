#pragma once

// ---------------------------------------------------------------------------
// 0D well-stirred constant-volume reactor with generic mass-action reactions.
//
//   dy_i/dt = sum_r (nu_p - nu_r)_i * k_r * prod_j c_j^nu_r_j
//
// with mass-action power-law rates.  Optional Arrhenius temperature scaling:
// when activation_energy > 0, k_r = pre_exponential * exp(-Ea/(R*T)) with
// R = 8.314 J/(mol*K); otherwise k_r = rate_constant.
//
// Integrated with the shared solver::integrate_step() integrator; species are
// clamped to zero after each step (a small negative excursion is clamped with
// a warning, a large one is a hard error).  Deterministic and exception-free.
//
// Lumped/0D modules are exempt from the per-domain channel rule (same as the
// engine and compressor modules).
// ---------------------------------------------------------------------------

#include <exd/physics/model_status.hpp>
#include <exd/physics/solver/integrators.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace exd::physics::chemistry {

/// A chemical species with an initial concentration.
struct SpeciesSpec
{
    std::string name;                    // diagnostic label only
    double initial_concentration = 1.0;  // mol/m^3, >= 0
};

/// One mass-action reaction (generalized Arrhenius optional).
struct ReactionSpec
{
    std::vector<int32_t> reactant_indices; // species indices
    std::vector<double> reactant_stoich;   // positive stoichimetric coefficients
    std::vector<int32_t> product_indices;
    std::vector<double> product_stoich;

    double rate_constant = 0.0;        // mass-action k, >= 0 (used when activation_energy == 0)
    double activation_energy = 0.0;    // Ea in J/mol, >= 0 (0 disables Arrhenius)
    double pre_exponential = 0.0;      // A for k = A*exp(-Ea/(R*T)) when Ea > 0, >= 0
};

/// Configuration for a batch reactor run.
struct ChemistryConfig
{
    std::vector<SpeciesSpec> species;  // >= 1
    std::vector<ReactionSpec> reactions; // may be empty (inert species)

    double temperature = 300.0;        // K, > 0 (Arrhenius reference)
    double dt = 1e-3;                  // s, > 0
    double end_time = 1.0;             // s, > 0

    solver::IntegratorConfig integration; // default RK4
};

/// Validate the reactor config.  Fatal problems return false and fill
/// `error`; non-fatal observations are appended to `warnings`.
bool validate_chemistry_config(const ChemistryConfig& config,
                               std::string& error,
                               std::vector<std::string>& warnings);

/// Result of a batch reactor run.
struct ChemistryResult
{
    bool ok = false;
    ModelStatus status;

    std::vector<double> final_concentrations;      // mol/m^3
    std::vector<std::vector<double>> history;      // concentration snapshots
    std::vector<double> time_history;              // matching times (s)

    double max_total_mass_error = 0.0; // max over steps of |sum(dc)|/(1+sum|c_old|)
};

/// Integrate the reactor.  Deterministic: identical configs give bit-identical
/// results.  Failures surface through `status` and the result.ok flag.
ChemistryResult solve_chemistry(const ChemistryConfig& config, ModelStatus& status);

} // namespace exd::physics::chemistry