#include "fdm_internal.hpp"
#include <algorithm>

namespace exd::engine::physics::fluid::fdm {

void apply_boundary_conditions(FDMGrid& g, const FDMConfig& config) {
    const int nx = g.nx;
    const int ny = g.ny;

    for (const auto& bc : config.boundary_conditions) {
        switch (bc.edge) {

        case BoundaryEdge::Left:
            switch (bc.type) {
            case FDMBoundaryType::Inlet:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(0, j)] = bc.u_value;
                    g.v[g.idx(0, j)] = bc.v_value;
                    g.p[g.idx(0, j)] = g.p[g.idx(1, j)];
                }
                break;
            case FDMBoundaryType::Wall:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(0, j)] = -g.u[g.idx(1, j)];  // linear interp to zero
                    g.v[g.idx(0, j)] = 0.0;
                    g.p[g.idx(0, j)] = g.p[g.idx(1, j)];
                }
                break;
            case FDMBoundaryType::Outlet:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(0, j)] = g.u[g.idx(1, j)];
                    g.v[g.idx(0, j)] = g.v[g.idx(1, j)];
                    g.p[g.idx(0, j)] = g.p[g.idx(1, j)];
                }
                break;
            case FDMBoundaryType::Symmetry:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(0, j)] = 0.0;
                    g.v[g.idx(0, j)] = g.v[g.idx(1, j)];
                    g.p[g.idx(0, j)] = g.p[g.idx(1, j)];
                }
                break;
            case FDMBoundaryType::Periodic:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(0, j)] = g.u[g.idx(nx, j)];
                    g.v[g.idx(0, j)] = g.v[g.idx(nx, j)];
                    g.p[g.idx(0, j)] = g.p[g.idx(nx, j)];
                }
                break;
            case FDMBoundaryType::FixedPressure:
                for (int j = 1; j <= ny; ++j) {
                    g.p[g.idx(0, j)] = bc.p_value;
                    g.u[g.idx(0, j)] = g.u[g.idx(1, j)];
                    g.v[g.idx(0, j)] = g.v[g.idx(1, j)];
                }
                break;
            }
            break;

        case BoundaryEdge::Right:
            switch (bc.type) {
            case FDMBoundaryType::Inlet:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(nx+1, j)] = bc.u_value;
                    g.v[g.idx(nx+1, j)] = bc.v_value;
                    g.p[g.idx(nx+1, j)] = g.p[g.idx(nx, j)];
                }
                break;
            case FDMBoundaryType::Wall:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(nx+1, j)] = -g.u[g.idx(nx, j)];
                    g.v[g.idx(nx+1, j)] = 0.0;
                    g.p[g.idx(nx+1, j)] = g.p[g.idx(nx, j)];
                }
                break;
            case FDMBoundaryType::Outlet:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(nx+1, j)] = g.u[g.idx(nx, j)];
                    g.v[g.idx(nx+1, j)] = g.v[g.idx(nx, j)];
                    g.p[g.idx(nx+1, j)] = g.p[g.idx(nx, j)];
                }
                break;
            case FDMBoundaryType::Symmetry:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(nx+1, j)] = 0.0;
                    g.v[g.idx(nx+1, j)] = g.v[g.idx(nx, j)];
                    g.p[g.idx(nx+1, j)] = g.p[g.idx(nx, j)];
                }
                break;
            case FDMBoundaryType::Periodic:
                for (int j = 1; j <= ny; ++j) {
                    g.u[g.idx(nx+1, j)] = g.u[g.idx(1, j)];
                    g.v[g.idx(nx+1, j)] = g.v[g.idx(1, j)];
                    g.p[g.idx(nx+1, j)] = g.p[g.idx(1, j)];
                }
                break;
            case FDMBoundaryType::FixedPressure:
                for (int j = 1; j <= ny; ++j) {
                    g.p[g.idx(nx+1, j)] = bc.p_value;
                    g.u[g.idx(nx+1, j)] = g.u[g.idx(nx, j)];
                    g.v[g.idx(nx+1, j)] = g.v[g.idx(nx, j)];
                }
                break;
            }
            break;

        case BoundaryEdge::Bottom:
            switch (bc.type) {
            case FDMBoundaryType::Inlet:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, 0)] = bc.u_value;
                    g.v[g.idx(i, 0)] = bc.v_value;
                    g.p[g.idx(i, 0)] = g.p[g.idx(i, 1)];
                }
                break;
            case FDMBoundaryType::Wall:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, 0)] = 0.0;
                    g.v[g.idx(i, 0)] = -g.v[g.idx(i, 1)];
                    g.p[g.idx(i, 0)] = g.p[g.idx(i, 1)];
                }
                break;
            case FDMBoundaryType::Outlet:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, 0)] = g.u[g.idx(i, 1)];
                    g.v[g.idx(i, 0)] = g.v[g.idx(i, 1)];
                    g.p[g.idx(i, 0)] = g.p[g.idx(i, 1)];
                }
                break;
            case FDMBoundaryType::Symmetry:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, 0)] = g.u[g.idx(i, 1)];
                    g.v[g.idx(i, 0)] = 0.0;
                    g.p[g.idx(i, 0)] = g.p[g.idx(i, 1)];
                }
                break;
            case FDMBoundaryType::Periodic:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, 0)] = g.u[g.idx(i, ny)];
                    g.v[g.idx(i, 0)] = g.v[g.idx(i, ny)];
                    g.p[g.idx(i, 0)] = g.p[g.idx(i, ny)];
                }
                break;
            case FDMBoundaryType::FixedPressure:
                for (int i = 1; i <= nx; ++i) {
                    g.p[g.idx(i, 0)] = bc.p_value;
                    g.u[g.idx(i, 0)] = g.u[g.idx(i, 1)];
                    g.v[g.idx(i, 0)] = g.v[g.idx(i, 1)];
                }
                break;
            }
            break;

        case BoundaryEdge::Top:
            switch (bc.type) {
            case FDMBoundaryType::Inlet:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, ny+1)] = bc.u_value;
                    g.v[g.idx(i, ny+1)] = bc.v_value;
                    g.p[g.idx(i, ny+1)] = g.p[g.idx(i, ny)];
                }
                break;
            case FDMBoundaryType::Wall:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, ny+1)] = 0.0;
                    g.v[g.idx(i, ny+1)] = -g.v[g.idx(i, ny)];
                    g.p[g.idx(i, ny+1)] = g.p[g.idx(i, ny)];
                }
                break;
            case FDMBoundaryType::Outlet:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, ny+1)] = g.u[g.idx(i, ny)];
                    g.v[g.idx(i, ny+1)] = g.v[g.idx(i, ny)];
                    g.p[g.idx(i, ny+1)] = g.p[g.idx(i, ny)];
                }
                break;
            case FDMBoundaryType::Symmetry:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, ny+1)] = g.u[g.idx(i, ny)];
                    g.v[g.idx(i, ny+1)] = 0.0;
                    g.p[g.idx(i, ny+1)] = g.p[g.idx(i, ny)];
                }
                break;
            case FDMBoundaryType::Periodic:
                for (int i = 1; i <= nx; ++i) {
                    g.u[g.idx(i, ny+1)] = g.u[g.idx(i, 1)];
                    g.v[g.idx(i, ny+1)] = g.v[g.idx(i, 1)];
                    g.p[g.idx(i, ny+1)] = g.p[g.idx(i, 1)];
                }
                break;
            case FDMBoundaryType::FixedPressure:
                for (int i = 1; i <= nx; ++i) {
                    g.p[g.idx(i, ny+1)] = bc.p_value;
                    g.u[g.idx(i, ny+1)] = g.u[g.idx(i, ny)];
                    g.v[g.idx(i, ny+1)] = g.v[g.idx(i, ny)];
                }
                break;
            }
            break;
        }
    }
}

} // namespace exd::engine::physics::fluid::fdm
