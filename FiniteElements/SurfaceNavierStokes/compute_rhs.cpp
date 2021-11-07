#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

Eigen::VectorXd SurfaceNavierStokesSolver::compute_rhs()
{
    auto rhs = Eigen::VectorXd(m_system_N);
    for (int i = 0; i < m_system_N; i++) rhs[i] = 0.;

    double inv_dt = 1./m_current_time_step_dt;
    
    // For each velocity trial function on a vertex node.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

        vec3 integral = vec3(0., 0., 0.);

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
            double tri_area = geom.triangle_area(tri);
            P2Element elements[6] = {
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            double element_weights[6] = {
                -1./360., -1./360., 1./60., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * velocity[elements[i]];
            }

            /*--------------------------------------------------------------------------------
                Source term.
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }

            he = he.twin().next();
        } while (he != start);

        rhs[3*velocity_node_indices[v] + 0] += integral.x();
        rhs[3*velocity_node_indices[v] + 1] += integral.y();
        rhs[3*velocity_node_indices[v] + 2] += integral.z();
    }
    
    // For each velocity trial function on an edge node.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        vec3 integral = vec3(0., 0., 0.);

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

            P2Element elements[6] = {
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            /*--------------------------------------------------------------------------------
                Time-step update term.
            dot(-u_prev/dt, psi^u)
            --------------------------------------------------------------------------------*/
            double element_weights[6] = {
                0, 0, -1./90., 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * velocity[elements[i]];
            }

            /*--------------------------------------------------------------------------------
                Source term.
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }
        }
        rhs[3*velocity_node_indices[edge] + 0] += integral.x();
        rhs[3*velocity_node_indices[edge] + 1] += integral.y();
        rhs[3*velocity_node_indices[edge] + 2] += integral.z();
    }
    return rhs;
}


