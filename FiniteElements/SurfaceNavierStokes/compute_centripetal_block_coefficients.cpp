#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

std::vector<CentripetalBlockEntry> SurfaceNavierStokesSolver::compute_centripetal_block_coefficients()
{
    auto coeffs = std::vector<CentripetalBlockEntry>();
    
    // For each basis trial function psi^r ...
    //------------------------------------------------------------
    // For each psi^r (based on a vertex)
    for (auto v : geom.mesh.vertices()) {
        vec3 v_pos = eigen_to_vec3(geom.position[v]);

        // For each triangle.
        auto start = v.halfedge(); // If v is a boundary vertex, this should correspond to a triangle and be on the boundary.
        auto he = start;
        do {
            auto tri = he.face();
	    assert(!tri.null());
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            vec3 vp_pos = eigen_to_vec3(geom.position[vp]);
            vec3 vpp_pos = eigen_to_vec3(geom.position[vpp]);
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

            #if 0
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
            #endif

            #if 0
            double weights[6*6] = {
// 0, 200:
    1./420., 1./2520., -1./630., -1./630., -1./210., -1./315.,
// 0., 020:
    1./2520., 1./420., -1./630., -1./630., -1./315., -1./210.,
// 0., 002:
    -1./630., -1./630., 1./84., -1./630., 1./210., 1./210.,
// 0., 110:
    -1./630., -1./630., -1./630., 4./315., 4./315., 4./315.,
// 0., 011:
    -1./210., -1./315., 1./210., 4./315., 4./105., 2./105.,
// 0., 101:
    -1./315., -1./210., 1./210., 4./315., 2./105., 4./105.,
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
	        vec3 integral = vec3(0,0,0);
                for (int j = 0; j < 6; j++) {
                    integral += weights[6*i + j] * normal[elements[j]];
                }
	        coeffs.emplace_back(v, elements[i], 2*tri_area * integral);
            }
            #endif


            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    return coeffs;
}
