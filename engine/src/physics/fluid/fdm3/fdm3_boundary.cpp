#include "fdm3_internal.hpp"
#include <algorithm>

namespace exd::engine::physics::fluid::fdm3 {

namespace {

// ─── face geometry helpers ──────────────────────────────────────────

struct FaceInfo {
    int axis;   // 0 = X-normal face, 1 = Y-normal, 2 = Z-normal
    int side;   // 0 = min boundary (ghost index 0), 1 = max (index n+1)
};

FaceInfo face_info(BoundaryFace f) {
    switch (f) {
    case BoundaryFace::XMin: return {0, 0};
    case BoundaryFace::XMax: return {0, 1};
    case BoundaryFace::YMin: return {1, 0};
    case BoundaryFace::YMax: return {1, 1};
    case BoundaryFace::ZMin: return {2, 0};
    case BoundaryFace::ZMax: return {2, 1};
    }
    return {0, 0};  // unreachable; keeps compilers quiet
}

int grid_dim(const FDM3Grid& g, int axis) {
    if (axis == 0) return g.nx;
    if (axis == 1) return g.ny;
    return g.nz;
}

// Ghost cell index along `axis` for a face on `side`.
int ghost_coord(int n, int side) { return side == 0 ? 0 : n + 1; }
// First interior cell index along `axis` for a face on `side`.
int inner_coord(int n, int side) { return side == 0 ? 1 : n; }
// Interior index whose value maps onto the ghost cell (periodic pairing).
int periodic_source_coord(int n, int side) { return side == 0 ? n : 1; }

// Velocity component reference by axis (0=u, 1=v, 2=w).
double& comp_ref(std::vector<double>& u, std::vector<double>& v,
                 std::vector<double>& w, int axis, size_t id) {
    if (axis == 0) return u[id];
    if (axis == 1) return v[id];
    return w[id];
}

} // anonymous namespace

// ────────────────────────────────────────────────────────────────────
// Boundary conditions on all six faces.
//
// Applied generically: per face we derive the normal axis (the velocity
// component whose direction is perpendicular to the face) and the two
// tangential axes, so the exact same per-type logic runs on every face
// without per-axis copy-paste.
// ────────────────────────────────────────────────────────────────────

void apply_boundary_conditions(FDM3Grid& g, const FDM3Config& config) {
    for (const auto& bc : config.boundary_conditions) {
        const FaceInfo fi = face_info(bc.face);
        const int axis = fi.axis;
        const int side = fi.side;
        const int t1 = (axis + 1) % 3;
        const int t2 = (axis + 2) % 3;
        const int n_axis = grid_dim(g, axis);
        const int n1 = grid_dim(g, t1);
        const int n2 = grid_dim(g, t2);

        const int gc = ghost_coord(n_axis, side);
        const int ic = inner_coord(n_axis, side);
        const int pc = periodic_source_coord(n_axis, side);

        for (int b = 1; b <= n2; ++b) {
            for (int a = 1; a <= n1; ++a) {
                int coord[3] = {0, 0, 0};
                coord[axis] = gc; coord[t1] = a; coord[t2] = b;
                const int i = coord[0], j = coord[1], k = coord[2];

                int icoord[3] = {0, 0, 0};
                icoord[axis] = ic; icoord[t1] = a; icoord[t2] = b;
                const int ii = icoord[0], jj = icoord[1], kk = icoord[2];

                int pcoord[3] = {0, 0, 0};
                pcoord[axis] = pc; pcoord[t1] = a; pcoord[t2] = b;
                const int pi = pcoord[0], pj = pcoord[1], pk = pcoord[2];

                const size_t gid = g.idx(i, j, k);
                const size_t iid = g.idx(ii, jj, kk);
                const size_t pid = g.idx(pi, pj, pk);

                switch (bc.type) {

                case FDMBoundaryType::Inlet:
                    // Specified velocity components; pressure zero-gradient.
                    g.u[gid] = bc.u_value;
                    g.v[gid] = bc.v_value;
                    g.w[gid] = bc.w_value;
                    g.p[gid] = g.p[iid];
                    break;

                case FDMBoundaryType::Wall:
                    // No-slip collocated wall.  The normal component is
                    // reflected (ghost = -interior) so the value at the cell
                    // face vanishes; the tangential components are likewise
                    // reflected with sign flip.  This enforces u = 0 at the
                    // face to linear-interpolation order — the standard
                    // collocated approximation.
                    comp_ref(g.u, g.v, g.w, axis, gid) = -comp_ref(g.u, g.v, g.w, axis, iid);
                    comp_ref(g.u, g.v, g.w, t1, gid) = -comp_ref(g.u, g.v, g.w, t1, iid);
                    comp_ref(g.u, g.v, g.w, t2, gid) = -comp_ref(g.u, g.v, g.w, t2, iid);
                    g.p[gid] = g.p[iid];
                    break;

                case FDMBoundaryType::Outlet:
                    // Zero-gradient everything.
                    g.u[gid] = g.u[iid];
                    g.v[gid] = g.v[iid];
                    g.w[gid] = g.w[iid];
                    g.p[gid] = g.p[iid];
                    break;

                case FDMBoundaryType::Symmetry:
                    // Normal component reflected (zero normal velocity at
                    // the face); tangential components copied.
                    comp_ref(g.u, g.v, g.w, axis, gid) = -comp_ref(g.u, g.v, g.w, axis, iid);
                    comp_ref(g.u, g.v, g.w, t1, gid) = comp_ref(g.u, g.v, g.w, t1, iid);
                    comp_ref(g.u, g.v, g.w, t2, gid) = comp_ref(g.u, g.v, g.w, t2, iid);
                    g.p[gid] = g.p[iid];
                    break;

                case FDMBoundaryType::Periodic:
                    // Ghost maps onto the interior cell adjacent to the
                    // opposite face along the same axis (natural opposite
                    // pair, e.g. XMin <-> XMax).
                    g.u[gid] = g.u[pid];
                    g.v[gid] = g.v[pid];
                    g.w[gid] = g.w[pid];
                    g.p[gid] = g.p[pid];
                    break;

                case FDMBoundaryType::FixedPressure:
                    // Dirichlet pressure: the cell-face value is pinned to
                    // p_value, so the ghost value is 2*p_value - interior.
                    // Velocity ghosts are zero-gradient (free).
                    g.p[gid] = 2.0 * bc.p_value - g.p[iid];
                    g.u[gid] = g.u[iid];
                    g.v[gid] = g.v[iid];
                    g.w[gid] = g.w[iid];
                    break;
                }
            }
        }
    }
}

void update_periodic_field_ghosts(FDM3Grid& g, const FDM3Config& config,
                                  std::vector<double>& field) {
    for (const auto& bc : config.boundary_conditions) {
        if (bc.type != FDMBoundaryType::Periodic)
            continue;
        const FaceInfo fi = face_info(bc.face);
        const int axis = fi.axis;
        const int side = fi.side;
        const int t1 = (axis + 1) % 3;
        const int t2 = (axis + 2) % 3;
        const int n1 = grid_dim(g, t1);
        const int n2 = grid_dim(g, t2);

        const int gc = ghost_coord(grid_dim(g, axis), side);
        const int pc = periodic_source_coord(grid_dim(g, axis), side);

        for (int b = 1; b <= n2; ++b) {
            for (int a = 1; a <= n1; ++a) {
                int gcoord[3] = {0, 0, 0};
                gcoord[axis] = gc; gcoord[t1] = a; gcoord[t2] = b;
                int pcoord[3] = {0, 0, 0};
                pcoord[axis] = pc; pcoord[t1] = a; pcoord[t2] = b;
                field[g.idx(gcoord[0], gcoord[1], gcoord[2])] =
                    field[g.idx(pcoord[0], pcoord[1], pcoord[2])];
            }
        }
    }
}

} // namespace exd::engine::physics::fluid::fdm3