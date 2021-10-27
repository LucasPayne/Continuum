#include "NavierStokes/NavierStokesSolver.h"
#include "NavierStokes/core.h"

#define ADVECTION 1

void NavierStokesSolver::add_nonlinear_velocity_residual(P2Attachment<vec2> &velocity_residual)
{
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
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            /*--------------------------------------------------------------------------------
                Time-step update term.
            dot(-u_prev/dt, psi^u)
            --------------------------------------------------------------------------------*/
            element_weights = {
                -1./360., -1./360., 1./60., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * (-velocity_prev[elements[i]]);
            }

            /*--------------------------------------------------------------------------------
                Source term.
            -dot(source_function, psi^u).
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            element_weights = {
                -1./360., -1./360., 1./60., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += -2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }

#if ADVECTION
            /*--------------------------------------------------------------------------------
                Advection term.
            u * dot(grad(u), psi^u)
            --------------------------------------------------------------------------------*/
            vec2 advection_weights[6*6] = {
//generated================================================================================
-K1/140 - K2/140, 11*K1/2520 + 11*K2/2520, -K1/280 - K2/280, -2*K1/315 - 2*K2/315, -K1/126 - K2/126, -4*K1/315 - 4*K2/315,
11*K1/2520 + 11*K2/2520, -K1/140 - K2/140, -K1/280 - K2/280, -2*K1/315 - 2*K2/315, -4*K1/315 - 4*K2/315, -K1/126 - K2/126,
-K1/280 - K2/280, -K1/280 - K2/280, 13*K1/420 + 13*K2/420, K1/210 + K2/210, 2*K1/105 + 2*K2/105, 2*K1/105 + 2*K2/105,
-2*K1/315 - 2*K2/315, -2*K1/315 - 2*K2/315, K1/210 + K2/210, -4*K1/105 - 4*K2/105, 2*K1/315 + 2*K2/315, 2*K1/315 + 2*K2/315,
-K1/126 - K2/126, -4*K1/315 - 4*K2/315, 2*K1/105 + 2*K2/105, 2*K1/315 + 2*K2/315, 4*K1/63 + 4*K2/63, 2*K1/63 + 2*K2/63,
-4*K1/315 - 4*K2/315, -K1/126 - K2/126, 2*K1/105 + 2*K2/105, 2*K1/315 + 2*K2/315, 2*K1/63 + 2*K2/63, 4*K1/63 + 4*K2/63,
//=========================================================================================
            };
            for (int i = 0; i < 6; i++) {
                for (int k = 0; k < 6; k++) {
                    double term = vec2::dot(velocity[elements[k]], advection_weights[6*i + k]);
                    integral += velocity[elements[i]] * term;
                }
            }
#endif



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
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            /*--------------------------------------------------------------------------------
                Time-step update term.
            dot(-u_prev/dt, psi^u)
            --------------------------------------------------------------------------------*/
            element_weights = {
                0, 0, -1./90., 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += 2 * tri_area * element_weights[i] * inv_dt * (-velocity_prev[elements[i]]);
            }

            /*--------------------------------------------------------------------------------
                Source term.
            -dot(source_function, psi^u).
                Approximate integration by samples.
            --------------------------------------------------------------------------------*/
            element_weights = {
                0, 0, -1./90., 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                integral += -2 * tri_area * element_weights[i] * source_samples_P2[elements[i]];
            }

#if ADVECTION
            /*--------------------------------------------------------------------------------
                Advection term.
            u * dot(grad(u), psi^u)
            --------------------------------------------------------------------------------*/
            vec2 advection_weights[6*6] = {
//generated================================================================================
-K1/105 - K2/21, 2*K1/315 + 2*K2/315, -K1/630 + 2*K2/315, 4*K1/315 - 2*K2/105, 2*K1/105 + 2*K2/315, 2*K1/315 - 2*K2/105,
2*K1/315 + 2*K2/315, -K1/21 - K2/105, 2*K1/315 - K2/630, -2*K1/105 + 4*K2/315, -2*K1/105 + 2*K2/315, 2*K1/315 + 2*K2/105,
-K1/630 + 2*K2/315, 2*K1/315 - K2/630, -K1/105 - K2/105, 2*K1/105 + 2*K2/105, 4*K1/315 + 2*K2/315, 2*K1/315 + 4*K2/315,
4*K1/315 - 2*K2/105, -2*K1/105 + 4*K2/315, 2*K1/105 + 2*K2/105, -16*K1/105 - 16*K2/105, -8*K1/105 - 16*K2/315, -16*K1/315 - 8*K2/105,
2*K1/105 + 2*K2/315, -2*K1/105 + 2*K2/315, 4*K1/315 + 2*K2/315, -8*K1/105 - 16*K2/315, -16*K1/105 - 16*K2/315, -16*K1/315 - 16*K2/315,
2*K1/315 - 2*K2/105, 2*K1/315 + 2*K2/105, 2*K1/315 + 4*K2/315, -16*K1/315 - 8*K2/105, -16*K1/315 - 16*K2/315, -16*K1/315 - 16*K2/105,
//=========================================================================================
            };
            for (int i = 0; i < 6; i++) {
                for (int k = 0; k < 6; k++) {
                    double term = vec2::dot(velocity[elements[k]], advection_weights[6*i + k]);
                    integral += velocity[elements[i]] * term;
                }
            }
#endif
        }
        velocity_residual[edge] += integral;
    }
}


void NavierStokesSolver::add_nonlinear_term_matrix_top_left(std::vector<TopLeftEntry> &coeffs)
{
    // For each velocity trial basis function psi^u
    //------------------------------------------------------------
    // For each psi^u based on a vertex.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

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
            Edge edge_110 = he.next().edge();
            Edge edge_011 = he.next().next().edge();
            Edge edge_101 = he.edge();
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double tri_area = geom.triangle_area(tri);

            /*--------------------------------------------------------------------------------
                Advection term.
            --------------------------------------------------------------------------------*/
            P2Element elements[6] = {
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            vec2 terms[6*6] = {
//generated================================================================================
-K1/280, 11*K2/2520, K1/280 + K2/280, K1/630 - 2*K2/315, -K1/630 - K2/126, 2*K2/315,
11*K1/2520, -K2/280, K1/280 + K2/280, -2*K1/315 + K2/630, 2*K1/315, -K1/126 - K2/630,
-K1/140, -K2/140, -13*K1/420 - 13*K2/420, K1/105 + K2/105, -K1/105 + 4*K2/105, 4*K1/105 - K2/105,
-K1/126, -K2/126, -K1/210 - K2/210, -2*K1/105 - 2*K2/105, 2*K1/105 + 4*K2/315, 4*K1/315 + 2*K2/105,
-2*K1/315, -4*K2/315, -2*K1/105 - 2*K2/105, -4*K1/315 - 2*K2/315, 4*K1/315 + 2*K2/63, 8*K1/315 + 2*K2/315,
-4*K1/315, -2*K2/315, -2*K1/105 - 2*K2/105, -2*K1/315 - 4*K2/315, 2*K1/315 + 8*K2/315, 2*K1/63 + 4*K2/315,
//=========================================================================================
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                vec2 total_x = vec2(0., 0.);
                vec2 total_y = vec2(0., 0.);
                for (int k = 0; k < 6; k++) {
                    total_x += velocity[elements[k]].x() * terms[6*i + k];
                    total_y += velocity[elements[k]].y() * terms[6*i + k];
                }
	        coeffs.emplace_back(P2Element(v), 0, P2Element(elements[i]), 0, total_x.x());
	        coeffs.emplace_back(P2Element(v), 0, P2Element(elements[i]), 1, total_x.y());
	        coeffs.emplace_back(P2Element(v), 1, P2Element(elements[i]), 0, total_y.x());
	        coeffs.emplace_back(P2Element(v), 1, P2Element(elements[i]), 1, total_y.y());
            }


            he = he.twin().next();
        } while (he != start);
    }
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            auto edge_110 = he.edge();
            auto edge_011 = he.next().edge();
            auto edge_101 = he.next().next().edge();
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double tri_area = geom.triangle_area(tri);

            /*--------------------------------------------------------------------------------
                Advection term.
            --------------------------------------------------------------------------------*/
            P2Element elements[6] = {
                vp, vpp, v, edge_110, edge_011, edge_101
            };
            vec2 terms[6*6] = {
//generated================================================================================
2*K1/105, -4*K2/315, 2*K1/315 + 2*K2/315, -4*K1/315 + 2*K2/105, 4*K1/315 + 2*K2/315, -8*K1/315 - 2*K2/105,
-4*K1/315, 2*K2/105, 2*K1/315 + 2*K2/315, 2*K1/105 - 4*K2/315, -2*K1/105 - 8*K2/315, 2*K1/315 + 4*K2/315,
-K1/126, -K2/126, -K1/210 - K2/210, -2*K1/105 - 2*K2/105, 2*K1/105 + 4*K2/315, 4*K1/315 + 2*K2/105,
4*K1/63, 4*K2/63, 4*K1/105 + 4*K2/105, 16*K1/105 + 16*K2/105, -16*K1/105 - 32*K2/315, -32*K1/315 - 16*K2/105,
2*K1/315, 2*K2/63, -2*K1/315 - 2*K2/315, 8*K1/105 + 16*K2/315, -8*K1/105 - 8*K2/315, -16*K2/315,
2*K1/63, 2*K2/315, -2*K1/315 - 2*K2/315, 16*K1/315 + 8*K2/105, -16*K1/315, -8*K1/315 - 8*K2/105,
//=========================================================================================
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                vec2 total_x = vec2(0., 0.);
                vec2 total_y = vec2(0., 0.);
                for (int k = 0; k < 6; k++) {
                    total_x += velocity[elements[k]].x() * terms[6*i + k];
                    total_y += velocity[elements[k]].y() * terms[6*i + k];
                }
	        coeffs.emplace_back(P2Element(edge), 0, P2Element(elements[i]), 0, total_x.x());
	        coeffs.emplace_back(P2Element(edge), 0, P2Element(elements[i]), 1, total_x.y());
	        coeffs.emplace_back(P2Element(edge), 1, P2Element(elements[i]), 0, total_y.x());
	        coeffs.emplace_back(P2Element(edge), 1, P2Element(elements[i]), 1, total_y.y());
            }
        }
    }
}
