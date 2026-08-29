// 3D FDTD (Yee lattice) time-domain Maxwell solver.
// PEC box with a soft gaussian plane-wave source; step-driven, matching the
// house style of src/fluid/fdm/fdm_solver.cpp (config + in-place field +
// step result + valid/error status).
#include <exd/physics/electrical/fdtd.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace exd::physics::electrical {

namespace {

// ── Constants ──

constexpr double EPS0 = 8.8541878128e-12;  // vacuum permittivity (F/m)
constexpr double MU0 = 1.25663706212e-6;   // vacuum permeability (H/m)

// ── Grid helpers ──

size_t index_of(const std::array<int32_t, 3>& dims, int i, int j, int k) {
    const size_t nx = static_cast<size_t>(dims[0]);
    const size_t ny = static_cast<size_t>(dims[1]);
    return static_cast<size_t>(i)
         + nx * (static_cast<size_t>(j) + ny * static_cast<size_t>(k));
}

size_t cell_count(const std::array<int32_t, 3>& dims) {
    return static_cast<size_t>(dims[0])
         * static_cast<size_t>(dims[1])
         * static_cast<size_t>(dims[2]);
}

/// Read a neighbor cell; out-of-domain reads return zero, which supplies the
/// implicit (tangential E = 0) wall from the update formulas.
double read_cell(const std::vector<double>& v,
                 const std::array<int32_t, 3>& dims, int i, int j, int k) {
    if (i < 0 || i >= dims[0] || j < 0 || j >= dims[1] || k < 0 || k >= dims[2])
        return 0.0;
    return v[index_of(dims, i, j, k)];
}

// ── Validation ──

bool validate_config(const FdtdConfig& config, std::string& error) {
    for (int axis = 0; axis < 3; ++axis) {
        if (config.dims[axis] < 3) {
            error = "dims must be >= 3 on every axis";
            return false;
        }
        if (!(config.spacing[axis] > 0.0)) {
            error = "spacing must be positive on every axis";
            return false;
        }
    }
    if (!(config.eps_r > 0.0)) {
        error = "eps_r must be positive";
        return false;
    }
    if (!(config.mu_r > 0.0)) {
        error = "mu_r must be positive";
        return false;
    }
    if (!(config.courant_factor > 0.0 && config.courant_factor < 1.0)) {
        error = "courant_factor must be in (0, 1)";
        return false;
    }
    if (config.source_plane_index < 0 || config.source_plane_index >= config.dims[0]) {
        error = "source_plane_index must be in [0, dims.x)";
        return false;
    }
    return true;
}

/// dt is the configured value when fixed; otherwise courant_factor times the
/// 3D CFL limit 1/(c·sqrt(1/dx² + 1/dy² + 1/dz²)).
double time_step(const FdtdConfig& config) {
    if (config.dt > 0.0)
        return config.dt;
    const double c = 1.0 / std::sqrt(EPS0 * MU0);
    const double dx = config.spacing[0];
    const double dy = config.spacing[1];
    const double dz = config.spacing[2];
    const double inv_cfl = c * std::sqrt(1.0 / (dx * dx)
                                       + 1.0 / (dy * dy)
                                       + 1.0 / (dz * dz));
    return config.courant_factor / inv_cfl;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────
// Public API: init_fdtd_field
// ─────────────────────────────────────────────────────

bool init_fdtd_field(const FdtdConfig& config, FdtdField& field,
                     exd::physics::ModelStatus& status) {
    status.ok = true;
    status.error.clear();

    std::string error;
    if (!validate_config(config, error)) {
        status.ok = false;
        status.error = "invalid config: " + error;
        return false;
    }

    const size_t n = cell_count(config.dims);
    field.dims = config.dims;
    field.ex.assign(n, 0.0);
    field.ey.assign(n, 0.0);
    field.ez.assign(n, 0.0);
    field.hx.assign(n, 0.0);
    field.hy.assign(n, 0.0);
    field.hz.assign(n, 0.0);
    field.step = 0;
    field.t = 0.0;
    field.valid = true;
    return true;
}

// ─────────────────────────────────────────────────────
// Public API: step_fdtd
// ─────────────────────────────────────────────────────

bool step_fdtd(const FdtdConfig& config, FdtdField& field,
               FdtdStepResult& step_result, exd::physics::ModelStatus& status) {
    status.ok = true;
    status.error.clear();

    std::string error;
    if (!validate_config(config, error)) {
        status.ok = false;
        status.error = "invalid config: " + error;
        return false;
    }
    if (!field.valid) {
        status.ok = false;
        status.error = "field is not initialized";
        return false;
    }
    if (field.dims != config.dims) {
        status.ok = false;
        status.error = "field dims do not match config dims";
        return false;
    }
    const size_t n = cell_count(config.dims);
    if (field.ex.size() != n || field.ey.size() != n || field.ez.size() != n ||
        field.hx.size() != n || field.hy.size() != n || field.hz.size() != n) {
        status.ok = false;
        status.error = "field arrays are not sized to the config dims";
        return false;
    }

    const double dt = time_step(config);
    const double dx = config.spacing[0];
    const double dy = config.spacing[1];
    const double dz = config.spacing[2];
    const double eps = EPS0 * config.eps_r;
    const double mu = MU0 * config.mu_r;
    const double ht_h = dt / mu;
    const double ht_e = dt / eps;
    const int nx = config.dims[0];
    const int ny = config.dims[1];
    const int nz = config.dims[2];

    // H at n + 1/2: H -= (dt/mu) · curl(E)  (Faraday's law carries the minus)
    //   dHx = (dEz/dy - dEy/dz),  dHy = (dEx/dz - dEz/dx),  dHz = (dEy/dx - dEx/dy)
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const size_t id = index_of(config.dims, i, j, k);
                field.hx[id] -= ht_h * (
                    (read_cell(field.ez, config.dims, i, j + 1, k) - field.ez[id]) / dy -
                    (read_cell(field.ey, config.dims, i, j, k + 1) - field.ey[id]) / dz);
                field.hy[id] -= ht_h * (
                    (read_cell(field.ex, config.dims, i, j, k + 1) - field.ex[id]) / dz -
                    (read_cell(field.ez, config.dims, i + 1, j, k) - field.ez[id]) / dx);
                field.hz[id] -= ht_h * (
                    (read_cell(field.ey, config.dims, i + 1, j, k) - field.ey[id]) / dx -
                    (read_cell(field.ex, config.dims, i, j + 1, k) - field.ex[id]) / dy);
            }
        }
    }

    // Soft additive plane-wave source: gaussian gauss(n) into Ez at the
    // source x-plane before the E update, so the update sees the injected
    // field (additive, so ordering only shifts the injection by a half step).
    {
        const double arg = (static_cast<double>(field.step) - config.source_t0)
                         / config.source_sigma;
        const double gauss = std::isfinite(arg)
                           ? config.source_amplitude * std::exp(-(arg * arg))
                           : 0.0;
        const int i0 = config.source_plane_index;
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                field.ez[index_of(config.dims, i0, j, k)] += gauss;
            }
        }
    }

    // E at n + 1: E += (dt/eps) · curl(H)
    //   dEx = (dHz/dy - dHy/dz),  dEy = (dHx/dz - dHz/dx),  dEz = (dHy/dx - dHx/dy)
    // Out-of-range H reads are zero, matching the wall convention used in
    // the H update above (reads across the domain edge count as zero).
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const size_t id = index_of(config.dims, i, j, k);
                field.ex[id] += ht_e * (
                    (field.hz[id] - read_cell(field.hz, config.dims, i, j - 1, k)) / dy -
                    (field.hy[id] - read_cell(field.hy, config.dims, i, j, k - 1)) / dz);
                field.ey[id] += ht_e * (
                    (field.hx[id] - read_cell(field.hx, config.dims, i, j, k - 1)) / dz -
                    (field.hz[id] - read_cell(field.hz, config.dims, i - 1, j, k)) / dx);
                field.ez[id] += ht_e * (
                    (field.hy[id] - read_cell(field.hy, config.dims, i - 1, j, k)) / dx -
                    (field.hx[id] - read_cell(field.hx, config.dims, i, j - 1, k)) / dy);
            }
        }
    }

    // Diagnostics: peak |E|, peak |H|, total EM energy over all cells.
    double max_e = 0.0;
    double max_h = 0.0;
    double energy = 0.0;
    const double dv = dx * dy * dz;
    for (size_t i = 0; i < n; ++i) {
        max_e = std::max(max_e, std::max(std::abs(field.ex[i]),
                        std::max(std::abs(field.ey[i]), std::abs(field.ez[i]))));
        max_h = std::max(max_h, std::max(std::abs(field.hx[i]),
                        std::max(std::abs(field.hy[i]), std::abs(field.hz[i]))));
        if (config.record_energy) {
            energy += (eps * (field.ex[i] * field.ex[i]
                            + field.ey[i] * field.ey[i]
                            + field.ez[i] * field.ez[i])
                     + mu * (field.hx[i] * field.hx[i]
                           + field.hy[i] * field.hy[i]
                           + field.hz[i] * field.hz[i])) * dv;
        }
    }

    step_result.max_e = max_e;
    step_result.max_h = max_h;
    step_result.energy = energy;

    field.t += dt;
    field.step += 1;
    return true;
}

// ─────────────────────────────────────────────────────
// Public API: run_fdtd
// ─────────────────────────────────────────────────────

FdtdResult run_fdtd(const FdtdConfig& config) {
    FdtdResult result;
    exd::physics::ModelStatus status;

    FdtdField field;
    if (!init_fdtd_field(config, field, status)) {
        result.valid = false;
        result.error = status.error;
        return result;
    }

    for (int step = 0; step < config.max_steps; ++step) {
        FdtdStepResult step_result;
        if (!step_fdtd(config, field, step_result, status)) {
            result.valid = false;
            result.error = status.error;
            result.field = std::move(field);  // partial field on failure
            return result;
        }
        result.history.push_back(step_result);
    }

    result.valid = true;
    result.field = std::move(field);
    result.warnings = status.warnings;
    return result;
}

} // namespace exd::physics::electrical