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

#if 0
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
#endif
            



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

#if 0
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
#endif
            



} // endif m_use_advection
        }
        velocity_residual[edge] += integral;
    }
}


void NavierStokesSolver::explicit_advection()
{
    auto explicit_advection_vector_x_proj = Eigen::VectorXd(num_velocity_variation_nodes());
    auto explicit_advection_vector_y_proj = Eigen::VectorXd(num_velocity_variation_nodes());

    // For each velocity trial function on a vertex node.
    int vertex_index = 0;
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

        explicit_advection_vector_x_proj[vertex_index] = integral.x();
        explicit_advection_vector_y_proj[vertex_index] = integral.y();
        vertex_index += 1;
    }
    
    // For each velocity trial function on an edge node.
    int edge_index = 0;
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
                    vec2 val = velocity_prev[elements[i]] * vec2::dot(velocity_prev[vertices[K]], grad_weights[3*i + K]);
                    // printf("%.6g\n", m_current_time_step_dt);
                    // std::cout << velocity_prev[elements[i]] << "\n";
                    // std::cout << "::" << grad_weights[3*i+K] << "\n";
                    // printf("%.6f\n", vec2::dot(velocity_prev[vertices[K]], grad_weights[3*i + K]));
                    // std::cout << "::" << grad_weights[3*i+K] << "\n";
                    integral += val;
                }
            }
} // endif m_use_advection
        }
        explicit_advection_vector_x_proj[geom.mesh.num_interior_vertices() + edge_index] = integral.x();
        explicit_advection_vector_y_proj[geom.mesh.num_interior_vertices() + edge_index] = integral.y();
        edge_index += 1;
    }
    // getchar();

    // Project into velocity space.
    SparseMatrix gramian = gramian_matrix_P2_0();
    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double> > linear_solver;
    linear_solver.compute(gramian);
    Eigen::VectorXd explicit_advection_vector_x = linear_solver.solve(explicit_advection_vector_x_proj);
    Eigen::VectorXd explicit_advection_vector_y = linear_solver.solve(explicit_advection_vector_y_proj);

    // Add this advection term to the velocity on the mesh.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        velocity[v] += m_current_time_step_dt * vec2(explicit_advection_vector_x[counter], explicit_advection_vector_y[counter]);
        counter += 1;
    }
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) continue;
        velocity[e] += m_current_time_step_dt * vec2(explicit_advection_vector_x[counter], explicit_advection_vector_y[counter]);
        counter += 1;
    }
}

// Gramian projection matrix for P2_0 (zero on the boundary).
SparseMatrix NavierStokesSolver::gramian_matrix_P2_0()
{
    std::vector<EigenTriplet> coefficients;
    int N = num_velocity_variation_nodes();
    auto add_entry = [&](int i, int j, double value) {
        printf("%d %d in %dx%d\n", i, j, N, N);
        assert(i >= 0 && j >= 0 && i <= N-1 && j <= N-1);
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^us on a vertex.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto start = v.halfedge();
        auto he = start;
        int v_index = velocity_node_indices[v];
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = velocity_node_indices[vp];
            int vpp_index = velocity_node_indices[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            int edge_110_index = velocity_node_indices[edge_110];
            int edge_011_index = velocity_node_indices[edge_011];
            int edge_101_index = velocity_node_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

	    printf("vertex --- v\n");
            add_entry(v_index, v_index, (1./60.)*R);
            if (!vp.on_boundary()) {
                printf("vertex --- vp\n");
                add_entry(v_index, vp_index, (-1./360.)*R);
            }
            if (!vpp.on_boundary()) {
                printf("vertex --- vpp\n");
                add_entry(v_index, vpp_index, (-1./360.)*R);
            }

            if (!edge_110.on_boundary()) {
                printf("vertex --- edge_110\n");
                add_entry(v_index, edge_110_index, (-1./90.)*R);
            }
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }
    // For each psi^us at a midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = velocity_node_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // if (tri.null()) continue;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            int v_index = velocity_node_indices[v];
            int vp_index = velocity_node_indices[vp];
            int vpp_index = velocity_node_indices[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            int edge_110_index = velocity_node_indices[edge_110];
            int edge_011_index = velocity_node_indices[edge_011];
            int edge_101_index = velocity_node_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

	    if (!edge_110.on_boundary()) {
                printf("midpoint --- edge_110\n");
                add_entry(edge_110_index,
		      edge_110_index, (4./45.)*R);
            }
            if (!edge_011.on_boundary()) {
                printf("midpoint --- edge_011\n");
                add_entry(edge_110_index,
                      edge_011_index, (2./45.)*R);
            }
            
            if (!edge_101.on_boundary()) {
                printf("midpoint --- edge_101\n");
                add_entry(edge_110_index,
                      edge_101_index, (2./45.)*R);
            }
            
            if (!v.on_boundary()) {
                printf("midpoint --- v\n");
                add_entry(edge_110_index,
                      v_index, (-1./90.)*R);
            }
            
        }
    }

    auto matrix = SparseMatrix(N, N);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();

    // make_sparsity_image(matrix, DATA "P20_P20_gramian_navierstokes.ppm");

    return matrix;
}
