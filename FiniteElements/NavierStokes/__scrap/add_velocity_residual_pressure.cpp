#include "NavierStokes/NavierStokesSolver.h"

void NavierStokesSolver::add_velocity_residual_pressure(P2Attachment<vec2> &velocity_residual)
{
    /*--------------------------------------------------------------------------------
        Pressure term.
    -p * div(psi^u)
    --------------------------------------------------------------------------------*/

    // For each velocity trial function on a vertex node.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;

        vec2 integral = vec2(0., 0.);

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            Edge edge_101 = he.edge();
            Edge edge_110 = he.next().edge();
            Edge edge_011 = he.next().next().edge();
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double C = kinematic_viscosity * 1.0/(2.0*geom.triangle_area(tri));

            const Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            // 002, 002: K1**2/2 + K1*K2 + K2**2/2
            // 002, 020: K1*K2/6 + K2**2/6
            // 002, 200: K1**2/6 + K1*K2/6
            // 002, 011: -2*K1*K2/3 - 2*K2**2/3
            // 002, 110: 0
            // 002, 101: -2*K1**2/3 - 2*K1*K2/3
            const double element_weights[6] = {
                /* v        */ vec2::dot(K1, K1)/2 + vec2::dot(K1,K2) + vec2::dot(K2, K2)/2,
                /* vp       */ vec2::dot(K1, K1)/6 + vec2::dot(K1, K2)/6,
                /* vpp      */ vec2::dot(K1, K2)/6 + vec2::dot(K2, K2)/6,
                /* edge_110 */ 0,
                /* edge_011 */ -2*vec2::dot(K1, K2)/3 - 2*vec2::dot(K2, K2)/3,
                /* edge_101 */ -2*vec2::dot(K1, K1)/3 - 2*vec2::dot(K1, K2)/3,
            };

            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += C * element_weights[i] * velocity[elements[i]];
            }

            he = he.twin().next();
        } while (he != start);

        velocity_residual[v] += integral;
    }
    
    // For each velocity trial function on an edge node.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        vec2 integral = vec2(0., 0.);

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.tip();
            auto vpp = he.vertex();
            Edge edge_110 = he.edge();
            Edge edge_011 = he.next().edge();
            Edge edge_101 = he.next().next().edge();
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double C = kinematic_viscosity * 1.0/(2.0*geom.triangle_area(tri));

            const Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            // 110, 002: 0
            // 110, 020: 2*K1*K2/3
            // 110, 200: 2*K1*K2/3
            // 110, 011: -4*K1**2/3 - 4*K1*K2/3
            // 110, 110: 4*K1**2/3 + 4*K1*K2/3 + 4*K2**2/3
            // 110, 101: -4*K1*K2/3 - 4*K2**2/3
            const double element_weights[6] = {
                /* v        */ 0,
                /* vp       */ 2*vec2::dot(K1, K2)/3,
                /* vpp      */ 2*vec2::dot(K1, K2)/3,
                /* edge_110 */ 4*vec2::dot(K1, K1)/3 + 4*vec2::dot(K1, K2)/3 + 4*vec2::dot(K2, K2)/3,
                /* edge_011 */ -4*vec2::dot(K1, K1)/3 - 4*vec2::dot(K1, K2)/3,
                /* edge_101 */ -4*vec2::dot(K1, K2)/3 - 4*vec2::dot(K2, K2)/3,
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += C * element_weights[i] * velocity[elements[i]];
            }
        }

        velocity_residual[edge] += integral;
    }
}
