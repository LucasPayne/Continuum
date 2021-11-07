#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

std::vector<CentripetalBlockEntry> SurfaceNavierStokesSolver::compute_centripetal_block_coefficients()
{
    auto coeffs = std::vector<CentripetalBlockEntry>();
    
    // For each basis trial function psi^r ...
    //------------------------------------------------------------
    // For each psi^r (based on a vertex)
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        coeffs.emplace_back(v, v, normal[v]);
    }
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        coeffs.emplace_back(edge, edge, normal[edge]);
    }

    return coeffs;
}
