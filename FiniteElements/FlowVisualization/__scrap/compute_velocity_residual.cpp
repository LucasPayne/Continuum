#include "NavierStokes/NavierStokesSolver.h"
#include "NavierStokes/core.h"

void NavierStokesSolver::compute_velocity_residual(SparseMatrix &gateaux_matrix, P2Attachment<vec2> &velocity_residual)
{
    // Make sure that the residual begins at zero.
    for (auto v : geom.mesh.vertices()) {
        velocity_residual[v] = vec2(0,0);
    }
    for (auto e : geom.mesh.edges()) {
        velocity_residual[e] = vec2(0,0);
    }

    double inv_dt = 1./m_current_time_step_dt;

    // For each velocity trial function on a vertex node.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

        vec2 integral = vec2(0., 0.);

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            Edge edge_101 = he.edge();
            Edge edge_110 = he.next().edge();
            Edge edge_011 = he.next().next().edge();
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double tri_area = geom.triangle_area(tri);
            P2Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            /*--------------------------------------------------------------------------------
                Time-step update term.
            dot((u - u_prev)/dt, psi^u)
            --------------------------------------------------------------------------------*/
            element_weights = {
                1./60., -1./360., -1./360., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * (velocity[elements[i]] - velocity_prev[elements[i]]);
            }
            /*--------------------------------------------------------------------------------
                Viscosity term.
            kinematic_viscosity * grad(u):grad(psi^u)
            --------------------------------------------------------------------------------*/
            // 002, 002: K1**2/2 + K1*K2 + K2**2/2
            // 002, 020: K1*K2/6 + K2**2/6
            // 002, 200: K1**2/6 + K1*K2/6
            // 002, 011: -2*K1*K2/3 - 2*K2**2/3
            // 002, 110: 0
            // 002, 101: -2*K1**2/3 - 2*K1*K2/3
            element_weights = {
                /* v        */ vec2::dot(K1, K1)/2 + vec2::dot(K1,K2) + vec2::dot(K2, K2)/2,
                /* vp       */ vec2::dot(K1, K1)/6 + vec2::dot(K1, K2)/6,
                /* vpp      */ vec2::dot(K1, K2)/6 + vec2::dot(K2, K2)/6,
                /* edge_110 */ 0,
                /* edge_011 */ -2*vec2::dot(K1, K2)/3 - 2*vec2::dot(K2, K2)/3,
                /* edge_101 */ -2*vec2::dot(K1, K1)/3 - 2*vec2::dot(K1, K2)/3,
            };
            double viscosity_term_coefficient = m_kinematic_viscosity * 1.0/(2.0*tri_area);
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += viscosity_term_coefficient * element_weights[i] * velocity[elements[i]];
            }
            /*--------------------------------------------------------------------------------
                Pressure term.
            -p * div(psi^u)
            --------------------------------------------------------------------------------*/
            integral += pressure[v] * (K1/6 + K2/6);

            /*--------------------------------------------------------------------------------
                Source term.
            -dot(source_function, psi^u).
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            element_weights = {
                1./60., -1./360., -1./360., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += -2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }

            /*--------------------------------------------------------------------------------
                Advection term.
            dot(dot(u, grad(u)), psi^u)
            --------------------------------------------------------------------------------*/

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
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            Edge edge_110 = he.edge();
            Edge edge_011 = he.next().edge();
            Edge edge_101 = he.next().next().edge();
            double tri_area = geom.triangle_area(tri);
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            P2Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            /*--------------------------------------------------------------------------------
                Time-step update term.
            dot((u - u_prev)/dt, psi^u)
            --------------------------------------------------------------------------------*/
            element_weights = {
                -1./90., 0, 0, 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * (velocity[elements[i]] - velocity_prev[elements[i]]);
            }
            /*--------------------------------------------------------------------------------
                Viscosity term.
            kinematic_viscosity * grad(u):grad(psi^u)
            --------------------------------------------------------------------------------*/
            // 110, 002: 0
            // 110, 020: 2*K1*K2/3
            // 110, 200: 2*K1*K2/3
            // 110, 011: -4*K1**2/3 - 4*K1*K2/3
            // 110, 110: 4*K1**2/3 + 4*K1*K2/3 + 4*K2**2/3
            // 110, 101: -4*K1*K2/3 - 4*K2**2/3
            element_weights = {
                /* v        */ 0,
                /* vp       */ 2*vec2::dot(K1, K2)/3,
                /* vpp      */ 2*vec2::dot(K1, K2)/3,
                /* edge_110 */ 4*vec2::dot(K1, K1)/3 + 4*vec2::dot(K1, K2)/3 + 4*vec2::dot(K2, K2)/3,
                /* edge_011 */ -4*vec2::dot(K1, K1)/3 - 4*vec2::dot(K1, K2)/3,
                /* edge_101 */ -4*vec2::dot(K1, K2)/3 - 4*vec2::dot(K2, K2)/3,
            };
            double viscosity_term_coefficient = m_kinematic_viscosity * 1.0/(2.0*tri_area);
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += viscosity_term_coefficient * element_weights[i] * velocity[elements[i]];
            }
            
            /*--------------------------------------------------------------------------------
                Pressure term.
            -p * div(psi^u)
            --------------------------------------------------------------------------------*/
            integral += pressure[v] * (-K1/6 - K2/6);
            if (!vp.on_boundary()) {
                integral += pressure[vp] * (-K1/6 - K2/3);
            }
            if (!vpp.on_boundary()) {
                integral += pressure[vpp] * (-K1/3 - K2/6);
            }

            /*--------------------------------------------------------------------------------
                Source term.
            -dot(source_function, psi^u).
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            element_weights = {
                -1./90., 0, 0, 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += -2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }
        }
        velocity_residual[edge] += integral;
    }
}
