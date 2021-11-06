#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

std::vector<CentripetalBlockEntry> SurfaceNavierStokesSolver::compute_centripetal_block_coefficients()
{
    auto coeffs = std::vector<CentripetalBlockEntry>();
    
    // For each basis trial function psi^p ...
    //------------------------------------------------------------
    // For each psi^p (based on a vertex)
    Vertex last_vertex;
    for (auto v : geom.mesh.vertices()) {
        last_vertex = v;
    }
    for (auto v : geom.mesh.vertices()) {
        if (v == last_vertex) continue; // skip the last pressure node
        auto v_pos = geom.position[v];

        // For each triangle.
        auto start = v.halfedge(); // If v is a boundary vertex, this should correspond to a triangle and be on the boundary.
        auto he = start;
        do {
            auto tri = he.face();
	    assert(!tri.null());
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            vec3 K1 = v_pos - vpp_pos;
            vec3 K2 = vp_pos - v_pos;
            vec3 K3 = vpp_pos - vp_pos;

            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            vec3 tri_normal = triangle_normal[tri];
            double tri_area = geom.triangle_area(tri);
            
            P2Element elements[6] = {vp,vpp,v, edge_110, edge_011, edge_101};
            double weights[6] = {
                -1./120.,
                -1./120.,
                1./60.,
                1./30.,
                1./15.,
                1./15.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                coeffs.emplace_back(v, elements[i], 2*tri_area * weights[i] * tri_normal);
            }

            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    return coeffs;
}
