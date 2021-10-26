#include "NavierStokes/NavierStokesSolver.h"

std::vector<BottomLeftEntry> NavierStokesSolver::compute_gateaux_matrix_bottom_left()
{
    auto coeffs = std::vector<BottomLeftEntry>();
    
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
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            
            double R = 0.5/geom.triangle_area(tri);

            vec2 val = vec2(0,0);
            val = (-1./6.) * K3.perp();
            if (!v.on_boundary()) {
                coeffs.emplace_back(v, P2Element(v), 0, val.x());
                coeffs.emplace_back(v, P2Element(v), 1, val.y());
            }

            val = (1./6.) * K3.perp();
            if (!edge_110.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_110), 0, val.x());
                coeffs.emplace_back(v, P2Element(edge_110), 1, val.y());
            }

            val = (-1./6.) * (K2 - K1).perp();
            if (!edge_011.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_011), 0, val.x());
                coeffs.emplace_back(v, P2Element(edge_011), 1, val.y());
            }

            // Integrate psi^p at v against phi^u_101.
            val = (-1./6.) * (K1 - K2).perp();
            if (!edge_101.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_101), 0, val.x());
                coeffs.emplace_back(v, P2Element(edge_101), 1, val.y());
            }
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    return coeffs;
}
