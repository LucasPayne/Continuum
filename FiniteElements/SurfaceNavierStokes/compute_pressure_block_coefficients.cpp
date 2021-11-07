#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

std::vector<PressureBlockEntry> SurfaceNavierStokesSolver::compute_pressure_block_coefficients()
{
    auto coeffs = std::vector<PressureBlockEntry>();
    
    // For each basis trial function psi^p ...
    //------------------------------------------------------------
    // For each psi^p (based on a vertex)
    Vertex last_vertex;
    for (auto v : geom.mesh.vertices()) {
        last_vertex = v;
    }
    for (auto v : geom.mesh.vertices()) {
        if (v == last_vertex) continue; // skip the last pressure node
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

            // Gram-Schmidt to get a basis on this triangle.
            vec3 tri_basis_1 = (vp_pos - v_pos).normalized();
            vec3 tri_basis_2 = vpp_pos - v_pos;
            tri_basis_2 -= tri_basis_1 * vec3::dot(tri_basis_2, tri_basis_1);
            tri_basis_2 = tri_basis_2.normalized();
            
            // Rotate a vector 90 degrees anticlockwise on the triangle plane.
            auto perp = [&](vec3 vec)->vec3 {
                float a = vec3::dot(vec, tri_basis_1);
                float b = vec3::dot(vec, tri_basis_2);
                return -tri_basis_2*a + tri_basis_1*b;
            };
            
            double R = 0.5/geom.triangle_area(tri);
            mat3x3 tri_proj = triangle_projection_matrix[tri];

            vec3 val = vec3(0,0,0);
            val = (-1./6.) * perp(K3);
            if (!v.on_boundary()) {
                coeffs.emplace_back(v, P2Element(v), tri_proj * val);
            }

            val = (1./6.) * perp(K3);
            if (!edge_110.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_110), tri_proj * val);
            }

            val = (-1./6.) * perp(K2 - K1);
            if (!edge_011.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_011), tri_proj * val);
            }

            // Integrate psi^p at v against phi^u_101.
            val = (-1./6.) * perp(K1 - K2);
            if (!edge_101.on_boundary()) {
                coeffs.emplace_back(v, P2Element(edge_101), tri_proj * val);
            }
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    return coeffs;
}
