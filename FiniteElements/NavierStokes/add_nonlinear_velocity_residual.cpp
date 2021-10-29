#include "NavierStokes/NavierStokesSolver.h"
#include "NavierStokes/core.h"



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
            vec2 K1 = vec2_extract(v_pos - vpp_pos).perp();
            vec2 K2 = vec2_extract(vp_pos - v_pos).perp();
            vec2 K3 = vec2_extract(vpp_pos - vp_pos).perp(); //NOTE: perp!
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

if (m_use_advection) {
            /*--------------------------------------------------------------------------------
                Explicit advection term.
                This is approximated by taking a sample u_sample = u at psi^u.
            dot(u_sample, grad(u)) psi^(u,s)
            --------------------------------------------------------------------------------*/
            // vertex
            // 002,200: -A*K1**2/30 - B*K1*K2/30
            // 002,020: -A*K1*K2/30 - B*K2**2/30
            // 002,002: -A*K1**2/15 - A*K1*K2/15 - B*K1*K2/15 - B*K2**2/15
            // 002,110: -A*K1**2/30 - A*K1*K2/30 - B*K1*K2/30 - B*K2**2/30
            // 002,011: A*K1**2/30 + A*K1*K2/10 + B*K1*K2/30 + B*K2**2/10
            // 002,101: A*K1**2/10 + A*K1*K2/30 + B*K1*K2/10 + B*K2**2/30
            // midpoint
            // 110,200: A*K1**2/10 + B*K1*K2/10
            // 110,020: A*K1*K2/10 + B*K2**2/10
            // 110,002: A*K1**2/30 + A*K1*K2/30 + B*K1*K2/30 + B*K2**2/30
            // 110,110: 4*A*K1**2/15 + 4*A*K1*K2/15 + 4*B*K1*K2/15 + 4*B*K2**2/15
            // 110,011: -4*A*K1**2/15 - 2*A*K1*K2/15 - 4*B*K1*K2/15 - 2*B*K2**2/15
            // 110,101: -2*A*K1**2/15 - 4*A*K1*K2/15 - 2*B*K1*K2/15 - 4*B*K2**2/15
            // vec2 u_sample = velocity_prev[v];
            // double A = u_sample.x();
            // double B = u_sample.y();
            // float grad_weights[6] = {
            //     -A*K1*K1/30 - B*K1*K2/30,
            //     -A*K1*K2/30 - B*K2*K2/30,
            //     -A*K1*K1/15 - A*K1*K2/15 - B*K1*K2/15 - B*K2*K2/15,
            //     -A*K1*K1/30 - A*K1*K2/30 - B*K1*K2/30 - B*K2*K2/30,
            //     A*K1*K1/30 + A*K1*K2/10 + B*K1*K2/30 + B*K2*K2/10,
            //     A*K1*K1/10 + A*K1*K2/30 + B*K1*K2/10 + B*K2*K2/30
            // };
            // for (int i = 0; i < 6; i++) {
            //     vec2 u_val = velocity_prev[elements[i]];
            //     integral += (1/(2*tri_area)) * u_val * grad_weights[i];
            // }

            
            vec2 grad_weights[6*3] = {
                -K1/72 - K2/72, -K1/60 - K2/60, -K1/360 - K2/360,
                -K1/72 - K2/72, -K1/360 - K2/360, -K1/60 - K2/60,
                K1/20 + K2/20, K1/120 + K2/120, K1/120 + K2/120,
                K1/90 + K2/90, -K1/45 - K2/45, -K1/45 - K2/45,
                K1/15 + K2/15, K1/90 + K2/90, K1/45 + K2/45,
                K1/15 + K2/15, K1/45 + K2/45, K1/90 + K2/90,
            };
            Vertex vertices[3] = {v, vp, vpp};
            for (int i = 0; i < 6; i++) {
                for (int K = 0; K < 3; K++) {
                    integral += velocity_prev[elements[i]] * vec2::dot(velocity_prev[vertices[K]], grad_weights[3*i + K]);
                }
            }
            



} // endif m_use_advection
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
            vec2 K1 = vec2_extract(v_pos - vpp_pos).perp();
            vec2 K2 = vec2_extract(vp_pos - v_pos).perp();
            vec2 K3 = vec2_extract(vpp_pos - vp_pos).perp(); //NOTE: perp!

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

            /*--------------------------------------------------------------------------------
                Explicit advection term.
                This is approximated by taking a sample u_sample = u at psi^u.
            dot(u_sample, grad(u)) psi^(u,s)
            --------------------------------------------------------------------------------*/
if (m_use_advection) {
            // vec2 u_sample = velocity_prev[v];
            // double A = u_sample.x();
            // double B = u_sample.y();
            // float grad_weights[6] = {
            //     A*K1*K1/10 + B*K1*K2/10,
            //     A*K1*K2/10 + B*K2*K2/10,
            //     A*K1*K1/30 + A*K1*K2/30 + B*K1*K2/30 + B*K2*K2/30,
            //     4*A*K1*K1/15 + 4*A*K1*K2/15 + 4*B*K1*K2/15 + 4*B*K2*K2/15,
            //     -4*A*K1*K1/15 - 2*A*K1*K2/15 - 4*B*K1*K2/15 - 2*B*K2*K2/15,
            //     -2*A*K1*K1/15 - 4*A*K1*K2/15 - 2*B*K1*K2/15 - 4*B*K2*K2/15
            // };
            // for (int i = 0; i < 6; i++) {
            //     vec2 u_val = velocity_prev[elements[i]];
            //     integral += (1/(2*tri_area)) * u_val * grad_weights[i];
            // }

            vec2 grad_weights[6*3] = {
                K1/90, -K2/15, K1/45,
                K2/90, K2/45, -K1/15,
                vec2(0,0), K1/90 + K2/45, K1/45 + K2/90,
                -2*K1/45 - 2*K2/45, -4*K1/45 - 2*K2/15, -2*K1/15 - 4*K2/45,
                -4*K1/45 - 2*K2/45, -2*K1/45 - 2*K2/45, -2*K1/15 - 2*K2/45,
                -2*K1/45 - 4*K2/45, -2*K1/45 - 2*K2/15, -2*K1/45 - 2*K2/45,
            };
            Vertex vertices[3] = {v, vp, vpp};
            for (int i = 0; i < 6; i++) {
                for (int K = 0; K < 3; K++) {
                    integral += velocity_prev[elements[i]] * vec2::dot(velocity_prev[vertices[K]], grad_weights[3*i + K]);
                }
            }
            



} // endif m_use_advection
        }
        velocity_residual[edge] += integral;
    }
}


