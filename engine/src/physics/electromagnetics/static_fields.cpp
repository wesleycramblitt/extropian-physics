// Static field solvers (electrostatic φ and magnetostatic A_z) on a regular
// 3D grid using SOR for the discrete Poisson equation, following the FDM
// pressure-Poisson machinery.  Solves ∇²φ = 0 (E = −∇φ) or
// ∇²A_z = −μ₀·J_z (B = ∇×A) on a box centered on the origin with Dirichlet
// faces and optional internal electrode/current boxes.

#include <exd/engine/physics/electromagnetics/static_fields.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace exd::engine::physics::electromagnetics
{

namespace
{

// Vacuum permeability (H/m).
constexpr double MU0 = 1.25663706212e-6;

// Flat node index along the header convention idx(i,j,k) = i + nx·(j + ny·k).
inline std::size_t grid_index(int i, int j, int k, int nx, int ny)
{
    return static_cast<std::size_t>(i + nx * (j + ny * k));
}

// Finite-difference first derivative of `u` along one axis, central in the
// interior and one-sided on the boundary nodes.  `ia`/`length` are the node
// coordinate and axis length; `stride` is the flat offset between consecutive
// nodes along that axis and `h` the corresponding spacing.
inline double derivative(const std::vector<double>& u, std::size_t id,
                         int ia, int length, std::size_t stride, double h)
{
    if (ia == 0)
        return (u[id + stride] - u[id]) / h;
    if (ia == length - 1)
        return (u[id] - u[id - stride]) / h;
    return (u[id + stride] - u[id - stride]) / (2.0 * h);
}

// Is a world-space node position, expressed as `node[a] = origin[a] +
// index[a] * spacing[a]`, inside [center − half_extents, center + half_extents]?
bool node_in_box(const std::array<int32_t, 3>& index,
                 const std::array<double, 3>& origin,
                 const std::array<double, 3>& spacing,
                 const std::array<double, 3>& center,
                 const std::array<double, 3>& half_extents)
{
    for (int a = 0; a < 3; ++a)
    {
        const double pos = origin[a] + static_cast<double>(index[a]) * spacing[a];
        if (std::fabs(pos - center[a]) > half_extents[a])
            return false;
    }
    return true;
}

} // anonymous namespace

// ── Public solver ─────────────────────────────────────────────────

StaticFieldResult solve_static_field(const StaticFieldConfig& config)
{
    StaticFieldResult result;

    // ── Config validation ─────────────────────────────────────────
    for (int32_t d : config.dims)
        if (d < 2)
        {
            result.error = "static field: dims must be >= 2 per axis";
            return result;
        }
    for (double s : config.spacing)
        if (s <= 0.0)
        {
            result.error = "static field: spacing must be positive";
            return result;
        }
    if (config.sor_omega <= 0.0 || config.sor_omega >= 2.0)
    {
        result.error = "static field: sor_omega must be in (0, 2)";
        return result;
    }
    if (config.tolerance <= 0.0)
    {
        result.error = "static field: tolerance must be positive";
        return result;
    }
    if (config.max_iterations <= 0)
    {
        result.error = "static field: max_iterations must be positive";
        return result;
    }

    // ── Grid ──────────────────────────────────────────────────────
    // Node-centered grid with `dims` interior nodes per axis.  The box is
    // centered on the origin: origin[a] = −spacing[a]·(dims[a]−1)/2, so the
    // domain spans [−spacing·(dims−1)/2, +spacing·(dims−1)/2] per axis.
    const int nx = config.dims[0];
    const int ny = config.dims[1];
    const int nz = config.dims[2];
    const double hx = config.spacing[0];
    const double hy = config.spacing[1];
    const double hz = config.spacing[2];
    const std::array<double, 3> origin = {
        -hx * (nx - 1) / 2.0,
        -hy * (ny - 1) / 2.0,
        -hz * (nz - 1) / 2.0,
    };

    const std::size_t node_count = static_cast<std::size_t>(nx) *
                                   static_cast<std::size_t>(ny) *
                                   static_cast<std::size_t>(nz);

    std::vector<double> phi(node_count, 0.0);   // φ (V) or A_z (T·m)
    std::vector<double> rhs(node_count, 0.0);   // −μ₀·J_z for magnetostatic
    std::vector<bool> dirichlet(node_count, false);

    // ── Dirichlet faces ───────────────────────────────────────────
    // All nodes on the six outer faces get face_values in order
    // {x−, x+, y−, y+, z−, z+}; a corner node may be assigned by several
    // faces, later faces simply overwrite earlier assignments.
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t id = grid_index(i, j, k, nx, ny);
                if (i == 0) { phi[id] = config.face_values[0]; dirichlet[id] = true; }
                if (i == nx - 1) { phi[id] = config.face_values[1]; dirichlet[id] = true; }
                if (j == 0) { phi[id] = config.face_values[2]; dirichlet[id] = true; }
                if (j == ny - 1) { phi[id] = config.face_values[3]; dirichlet[id] = true; }
                if (k == 0) { phi[id] = config.face_values[4]; dirichlet[id] = true; }
                if (k == nz - 1) { phi[id] = config.face_values[5]; dirichlet[id] = true; }
            }

    // ── Dirichlet patches (electrodes / conductors) ───────────────
    // A node becomes a patch Dirichlet node when its cell — the node's volume
    // of influence, [position − spacing/2, position + spacing/2] per axis —
    // overlaps the patch box.  This keeps zero-thickness conductors (thinner
    // than one cell, e.g. the capacitor plates) from silently vanishing
    // between nodes.  Patches override the face values assigned above.
    for (const auto& patch : config.patches)
    {
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const std::array<int32_t, 3> index = {i, j, k};
                    bool inside = true;
                    for (int a = 0; a < 3; ++a)
                    {
                        const double pos = origin[a] +
                                           static_cast<double>(index[a]) * config.spacing[a];
                        if (std::fabs(pos - patch.center[a]) >
                            patch.half_extents[a] + config.spacing[a] / 2.0)
                            inside = false;
                    }
                    if (inside)
                    {
                        const std::size_t id = grid_index(i, j, k, nx, ny);
                        phi[id] = patch.value;
                        dirichlet[id] = true;
                    }
                }
    }

    // ── RHS (magnetostatic current boxes) ─────────────────────────
    // Each node whose position falls inside a current box carries
    // b = −μ₀·jz.  A current box whose full width is exactly one cell is
    // therefore captured by a single node, so the discrete source integral
    // reproduces I = jz·(2·half_extents_x)·(2·half_extents_y).
    if (config.mode == StaticFieldMode::Magnetostatic)
    {
        for (const auto& box : config.currents)
        {
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        const std::array<int32_t, 3> index = {i, j, k};
                        if (node_in_box(index, origin, config.spacing,
                                        box.center, box.half_extents))
                        {
                            const std::size_t id = grid_index(i, j, k, nx, ny);
                            rhs[id] += -MU0 * box.jz;
                        }
                    }
        }
    }

    // ── SOR iteration ─────────────────────────────────────────────
    // Discrete 7-point Poisson: at interior nodes that are not Dirichlet,
    //   phi_new = ((phi[i−1]+phi[i+1])/hx²
    //              + (phi[j−1]+phi[j+1])/hy²
    //              + (phi[k−1]+phi[k+1])/hz² − b) / (2/hx² + 2/hy² + 2/hz²)
    // with SOR relaxation phi += sor_omega·(phi_new − phi).  The sweep
    // residual is the maximum pre-update change |phi_new − phi|; the sweep
    // is converged when it drops below tolerance.  Dirichlet nodes are never
    // updated.  If convergence is not reached inside max_iterations the best
    // iterate is still returned with ok=true and a warning.
    const double hx2 = hx * hx;
    const double hy2 = hy * hy;
    const double hz2 = hz * hz;
    const double denom = 2.0 / hx2 + 2.0 / hy2 + 2.0 / hz2;
    const std::size_t nxs = static_cast<std::size_t>(nx);
    const std::size_t nxy = nxs * static_cast<std::size_t>(ny);

    int iterations = 0;
    double residual = 0.0;
    for (; iterations < config.max_iterations; ++iterations)
    {
        double sweep_max = 0.0;
        for (int k = 1; k < nz - 1; ++k)
            for (int j = 1; j < ny - 1; ++j)
                for (int i = 1; i < nx - 1; ++i)
                {
                    const std::size_t id = grid_index(i, j, k, nx, ny);
                    if (dirichlet[id])
                        continue;

                    const std::size_t ixm = id - 1;
                    const std::size_t ixp = id + 1;
                    const std::size_t jym = id - static_cast<std::size_t>(nx);
                    const std::size_t jyp = id + static_cast<std::size_t>(nx);
                    const std::size_t kzm = id - nxy;   // nxy = nx·ny
                    const std::size_t kzp = id + nxy;

                    const double phi_new = ((phi[ixm] + phi[ixp]) / hx2 +
                                            (phi[jym] + phi[jyp]) / hy2 +
                                            (phi[kzm] + phi[kzp]) / hz2 -
                                            rhs[id]) / denom;
                    const double delta = phi_new - phi[id];
                    sweep_max = std::max(sweep_max, std::fabs(delta));
                    phi[id] += config.sor_omega * delta;
                }
        residual = sweep_max;
        if (residual < config.tolerance)
            break;
    }
    if (iterations >= config.max_iterations)
        result.warnings.push_back("static field: SOR did not converge in " +
                                  std::to_string(config.max_iterations) +
                                  " iterations");

    // ── Derived field (central differences, one-sided on boundary) ─
    // Electrostatic:  E = −∇φ  →  Ex = −∂φ/∂x,  Ey = −∂φ/∂y,  Ez = −∂φ/∂z.
    // Magnetostatic:  B = ∇×A with A = A_z·ẑ
    //     Bx = ∂A_z/∂y,  By = −∂A_z/∂x,  Bz = 0.
    // The result grids share the same dims/spacing as the potential; boundary
    // cells use one-sided differences (Δ/h instead of central 2Δ/(2h)).
    std::vector<double> field(3 * node_count, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t id = grid_index(i, j, k, nx, ny);
                const double dphi_dx = derivative(phi, id, i, nx, 1, hx);
                const double dphi_dy = derivative(phi, id, j, ny, nxs, hy);
                const double dphi_dz = derivative(phi, id, k, nz, nxy, hz);

                if (config.mode == StaticFieldMode::Electrostatic)
                {
                    field[3 * id + 0] = -dphi_dx;
                    field[3 * id + 1] = -dphi_dy;
                    field[3 * id + 2] = -dphi_dz;
                }
                else
                {
                    field[3 * id + 0] = dphi_dy;    // ∂A_z/∂y
                    field[3 * id + 1] = -dphi_dx;   // −∂A_z/∂x
                    field[3 * id + 2] = 0.0;        // B_z ≡ 0 for A = A_z·ẑ
                }
            }

    // ── Result ────────────────────────────────────────────────────
    // NOTE: the discrete Poisson with RHS −μ₀·J and Dirichlet 0 on the box
    // faces approximates the infinite-domain wire only when the box is large;
    // the accompanying tests use a box large enough that the far field is
    // within tolerance of μ₀·I/(2πr).
    result.ok = true;
    result.iterations = iterations;
    result.residual = residual;
    result.potential = exd::engine::coupling::StructuredScalarGrid{origin, config.spacing,
                                                      config.dims, std::move(phi)};
    result.field_vector = exd::engine::coupling::StructuredVectorGrid{origin, config.spacing,
                                                         config.dims, std::move(field)};
    return result;
}

} // namespace exd::engine::physics::electromagnetics