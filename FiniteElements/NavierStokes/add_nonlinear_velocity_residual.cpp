#include "NavierStokes/NavierStokesSolver.h"
#include "core.h"

#define RESIDUAL_ADVECTION false
#define WEIGHTS_MODE 1

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
#if RESIDUAL_ADVECTION
            vec2 grad_weights[6*6] = {
#if WEIGHTS_MODE == 0
-K1/140 - K2/140, 11*K1/2520 + 11*K2/2520, -K1/280 - K2/280, -2*K1/315 - 2*K2/315, -K1/126 - K2/126, -4*K1/315 - 4*K2/315,
11*K1/2520 + 11*K2/2520, -K1/140 - K2/140, -K1/280 - K2/280, -2*K1/315 - 2*K2/315, -4*K1/315 - 4*K2/315, -K1/126 - K2/126,
-K1/280 - K2/280, -K1/280 - K2/280, 13*K1/420 + 13*K2/420, K1/210 + K2/210, 2*K1/105 + 2*K2/105, 2*K1/105 + 2*K2/105,
-2*K1/315 - 2*K2/315, -2*K1/315 - 2*K2/315, K1/210 + K2/210, -4*K1/105 - 4*K2/105, 2*K1/315 + 2*K2/315, 2*K1/315 + 2*K2/315,
-K1/126 - K2/126, -4*K1/315 - 4*K2/315, 2*K1/105 + 2*K2/105, 2*K1/315 + 2*K2/315, 4*K1/63 + 4*K2/63, 2*K1/63 + 2*K2/63,
-4*K1/315 - 4*K2/315, -K1/126 - K2/126, 2*K1/105 + 2*K2/105, 2*K1/315 + 2*K2/315, 2*K1/63 + 2*K2/63, 4*K1/63 + 4*K2/63,
#else
K1/280, -11*K1/2520, K1/140, K1/126, 2*K1/315, 4*K1/315,
-11*K2/2520, K2/280, K2/140, K2/126, 4*K2/315, 2*K2/315,
-K1/280 - K2/280, -K1/280 - K2/280, 13*K1/420 + 13*K2/420, K1/210 + K2/210, 2*K1/105 + 2*K2/105, 2*K1/105 + 2*K2/105,
-K1/630 + 2*K2/315, 2*K1/315 - K2/630, -K1/105 - K2/105, 2*K1/105 + 2*K2/105, 4*K1/315 + 2*K2/315, 2*K1/315 + 4*K2/315,
K1/630 + K2/126, -2*K1/315, K1/105 - 4*K2/105, -2*K1/105 - 4*K2/315, -4*K1/315 - 2*K2/63, -2*K1/315 - 8*K2/315,
-2*K2/315, K1/126 + K2/630, -4*K1/105 + K2/105, -4*K1/315 - 2*K2/105, -8*K1/315 - 2*K2/315, -2*K1/63 - 4*K2/315,
#endif
            };
            for (int i = 0; i < 6; i++) {
                for (int K = 0; K < 6; K++) {
                    vec2 u = velocity_prev[elements[i]];
                    // vec2 u_transformed = u.x()*side1 + u.y()*side2;
                    vec2 u_transformed = u;
                    integral += u_transformed * vec2::dot(velocity_prev[elements[K]], grad_weights[6*i + K]);
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

#if RESIDUAL_ADVECTION
            vec2 grad_weights[6*6] = {
#if WEIGHTS_MODE == 0
-K1/105 - K2/21, 2*K1/315 + 2*K2/315, -K1/630 + 2*K2/315, 4*K1/315 - 2*K2/105, 2*K1/105 + 2*K2/315, 2*K1/315 - 2*K2/105,
2*K1/315 + 2*K2/315, -K1/21 - K2/105, 2*K1/315 - K2/630, -2*K1/105 + 4*K2/315, -2*K1/105 + 2*K2/315, 2*K1/315 + 2*K2/105,
-K1/630 + 2*K2/315, 2*K1/315 - K2/630, -K1/105 - K2/105, 2*K1/105 + 2*K2/105, 4*K1/315 + 2*K2/315, 2*K1/315 + 4*K2/315,
4*K1/315 - 2*K2/105, -2*K1/105 + 4*K2/315, 2*K1/105 + 2*K2/105, -16*K1/105 - 16*K2/105, -8*K1/105 - 16*K2/315, -16*K1/315 - 8*K2/105,
2*K1/105 + 2*K2/315, -2*K1/105 + 2*K2/315, 4*K1/315 + 2*K2/315, -8*K1/105 - 16*K2/315, -16*K1/105 - 16*K2/315, -16*K1/315 - 16*K2/315,
2*K1/315 - 2*K2/105, 2*K1/315 + 2*K2/105, 2*K1/315 + 4*K2/315, -16*K1/315 - 8*K2/105, -16*K1/315 - 16*K2/315, -16*K1/315 - 16*K2/105,
#else
-2*K1/105, 4*K1/315, K1/126, -4*K1/63, -2*K1/315, -2*K1/63,
4*K2/315, -2*K2/105, K2/126, -4*K2/63, -2*K2/63, -2*K2/315,
-2*K1/315 - 2*K2/315, -2*K1/315 - 2*K2/315, K1/210 + K2/210, -4*K1/105 - 4*K2/105, 2*K1/315 + 2*K2/315, 2*K1/315 + 2*K2/315,
4*K1/315 - 2*K2/105, -2*K1/105 + 4*K2/315, 2*K1/105 + 2*K2/105, -16*K1/105 - 16*K2/105, -8*K1/105 - 16*K2/315, -16*K1/315 - 8*K2/105,
-4*K1/315 - 2*K2/315, 2*K1/105 + 8*K2/315, -2*K1/105 - 4*K2/315, 16*K1/105 + 32*K2/315, 8*K1/105 + 8*K2/315, 16*K1/315,
8*K1/315 + 2*K2/105, -2*K1/315 - 4*K2/315, -4*K1/315 - 2*K2/105, 32*K1/315 + 16*K2/105, 16*K2/315, 8*K1/315 + 8*K2/105,
#endif
            };
            for (int i = 0; i < 6; i++) {
                for (int K = 0; K < 6; K++) {
                    // integral += velocity_prev[elements[i]] * vec2::dot(velocity_prev[elements[K]], grad_weights[6*i + K]);
                    vec2 u = velocity_prev[elements[i]];
                    // vec2 u_transformed = u.x()*side1 + u.y()*side2;
                    vec2 u_transformed = u;
                    integral += u_transformed * vec2::dot(velocity_prev[elements[K]], grad_weights[6*i + K]);
                }
            }
#endif



} // endif m_use_advection
        }
        velocity_residual[edge] += integral;
    }
}


// Gramian projection matrix for P2_0 (zero on the boundary).
SparseMatrix NavierStokesSolver::gramian_matrix_P2_0()
{
    std::vector<EigenTriplet> coefficients;
    int N = num_velocity_variation_nodes();
    auto add_entry = [&](int i, int j, double value) {
        // printf("%d %d in %dx%d\n", i, j, N, N);
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

	    // printf("vertex --- v\n");
            add_entry(v_index, v_index, (1./60.)*R);
            if (!vp.on_boundary()) {
                // printf("vertex --- vp\n");
                add_entry(v_index, vp_index, (-1./360.)*R);
            }
            if (!vpp.on_boundary()) {
                // printf("vertex --- vpp\n");
                add_entry(v_index, vpp_index, (-1./360.)*R);
            }

            if (!edge_110.on_boundary()) {
                // printf("vertex --- edge_110\n");
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
                // printf("midpoint --- edge_110\n");
                add_entry(edge_110_index,
		      edge_110_index, (4./45.)*R);
            }
            if (!edge_011.on_boundary()) {
                // printf("midpoint --- edge_011\n");
                add_entry(edge_110_index,
                      edge_011_index, (2./45.)*R);
            }
            
            if (!edge_101.on_boundary()) {
                // printf("midpoint --- edge_101\n");
                add_entry(edge_110_index,
                      edge_101_index, (2./45.)*R);
            }
            
            if (!v.on_boundary()) {
                // printf("midpoint --- v\n");
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


void NavierStokesSolver::explicit_advection_lagrangian()
{
    auto func = [&](P2Element element) {
        if (element.on_boundary()) return;
        vec2 position;
        if (element.is_vertex()) {
            position = vec2(geom.position[element.vertex].x(), geom.position[element.vertex].z());
        } else {
            auto midpoint = 0.5*geom.position[element.edge.a().vertex()] + 0.5*geom.position[element.edge.b().vertex()];
            position = vec2(midpoint.x(), midpoint.z());
        }
        vec2 vel = velocity_prev[element];
        vec2 prev_position = position - m_current_time_step_dt*vel;
        _test_point_1[element] = vec3(prev_position.x(), 0, prev_position.y());

        int prev_i = floor((0.5*prev_position.x() + 0.5) * velocity_grid_N);
        int prev_j = floor((0.5*prev_position.y() + 0.5) * velocity_grid_N);
    
        if (prev_i < 0) prev_i = 0;
        if (prev_j < 0) prev_j = 0;
        if (prev_i >= velocity_grid_N) prev_i = velocity_grid_N-1;
        if (prev_j >= velocity_grid_N) prev_j = velocity_grid_N-1;

        vec2 prev_vel = velocity_grid_samples[velocity_grid_N*prev_j + prev_i];
        velocity[element] = prev_vel;
    };
    for (auto v : geom.mesh.vertices()) {
        func(P2Element(v));
    }
    for (auto e : geom.mesh.edges()) {
        func(P2Element(e));
    }
}


std::tuple<Face, vec3> NavierStokesSolver::traverse(Face tri, vec3 origin, vec3 shift, int depth)
{
    const int MAX_DEPTH = 50;
    if (depth == MAX_DEPTH) return {tri, origin};

    Halfedge hes[3] = {
        tri.halfedge(),
        tri.halfedge().next(),
        tri.halfedge().next().next()
    };
    Vertex vs[3] = {
        tri.halfedge().vertex(),
        tri.halfedge().next().vertex(),
        tri.halfedge().next().next().vertex()
    };
    vec3 ps[3];
    for (int i = 0; i < 3; i++) {
        ps[i] = eigen_to_vec3(geom.position[vs[i]]);
    }
    // Gram-Schmidt for a basis on this triangle.
    vec3 e1 = (ps[1] - ps[0]).normalized();
    vec3 e2 = ps[2] - ps[0];
    e2 -= e1*vec3::dot(e2, e1);
    e2 = e2.normalized();

    // Express the triangle in this basis.
    vec2 ps_t[3];
    for (int i = 0; i < 3; i++) {
        ps_t[i] = vec2(vec3::dot(ps[i]-ps[0], e1), vec3::dot(ps[i]-ps[0], e2));
    }
    // Create a ray.
    vec2 o_t = vec2(vec3::dot(origin - ps[0], e1), vec3::dot(origin-ps[0], e2));
    vec2 d_t = vec2(vec3::dot(shift, e1), vec3::dot(shift, e2));

    int line_hit_index = -1;
    double min_t = std::numeric_limits<double>::max();
    // Intersect each line determined by the triangle sides.
    for (int i = 0; i < 3; i++) {
        vec2 line_n = (ps_t[(i+1)%3] - ps_t[i]).perp();
        vec2 line_p = ps_t[i];
        double t = vec2::dot(line_p - o_t, line_n)/vec2::dot(d_t, line_n);
        if (t >= 0 && t < min_t) {
            min_t = t;
            line_hit_index = i;
        }
    }

    if (line_hit_index == -1 || min_t > 1) {
        // Travel stops on this triangle.
        // printf("Travel stops.\n");
        return {tri, origin+shift};
    }
    Halfedge hit_he = hes[line_hit_index];
    if (hit_he.twin().face().null()) {
        // Hit the boundary. Stop at the boundary intersection.
        // printf("Hit the boundary.\n");
        return {tri, origin + min_t*shift};
    }
    // Travel proceeds on another triangle.
    // Create an orthonormal basis for each incident face to the edge being travelled over.
    // This basis shares the E1 vector.
    vec3 E1 = eigen_to_vec3(geom.vector(hit_he)).normalized();
    vec3 from_E2 = -eigen_to_vec3(geom.vector(hit_he.next()));
    from_E2 -= E1*vec3::dot(from_E2, E1);
    from_E2 = from_E2.normalized();
    vec3 to_E2 = -eigen_to_vec3(geom.vector(hit_he.twin().next().next()));
    to_E2 -= E1*vec3::dot(to_E2, E1);
    to_E2 = to_E2.normalized();
    
    vec3 new_shift = E1*vec3::dot((1-min_t)*shift, E1) + to_E2*vec3::dot((1-min_t)*shift, from_E2);

    double fix = 0.001; // (To prevent intersections with the edge just passed through.)
    // printf("Continuing travels...\n");
    return traverse(hit_he.twin().face(), origin + min_t*shift + new_shift*fix, (1-fix)*new_shift, depth+1);
}



void NavierStokesSolver::explicit_advection_traversal()
{
    std::function<double(double,double,double)> u_basis[6] = {
        [](double x, double y, double z)->double {
            return x - 2*x*y - 2*x*z;
        },
        [](double x, double y, double z)->double {
            return y - 2*y*z - 2*y*x;
        },
        [](double x, double y, double z)->double {
            return z - 2*z*x - 2*z*y;
        },
        [](double x, double y, double z)->double {
            return 4*x*y;
        },
        [](double x, double y, double z)->double {
            return 4*y*z;
        },
        [](double x, double y, double z)->double {
            return 4*z*x;
        },
    };


    auto new_velocity = P2Attachment<vec2>(geom.mesh);
    auto to_vec3 = [](vec2 vec)->vec3 { return vec3(vec.x(), 0, vec.y()); };
    auto to_vec2 = [](vec3 vec)->vec2 { return vec2(vec.x(), vec.z()); };

    auto barycentric_coeff = [&](vec3 c, vec3 pa, vec3 pb, vec3 nor)->double {
        vec3 cross = vec3::cross(pb - c, pa - c);
        if (vec3::dot(cross, nor) < 0) return cross.length();
        return -cross.length();
    };

    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            new_velocity[v] = vec2(0,0);
            continue;
        }
        vec3 u = to_vec3(velocity_prev[v]);
        vec3 n = normal[v];


        auto start = v.halfedge();
        auto he = start;
        vec3 p1 = eigen_to_vec3(geom.position[v]);
        do {
            auto tri = he.face();
            Vertex v1 = v;
            Vertex v2 = he.next().vertex();
            Vertex v3 = he.next().next().vertex();
            vec3 p2 = eigen_to_vec3(geom.position[v2]);
            vec3 p3 = eigen_to_vec3(geom.position[v3]);

            // Gram-Schmidt to get a basis on this triangle.
            vec3 tri_basis_1 = (p2 - p1).normalized();
            vec3 tri_basis_2 = p3 - p1;
            tri_basis_2 -= tri_basis_1 * vec3::dot(tri_basis_2, tri_basis_1);
            tri_basis_2 = tri_basis_2.normalized();
            // Rotate a vector 90 degrees anticlockwise on the triangle plane.
            auto perp = [&](vec3 vec)->vec3 {
                float a = vec3::dot(vec, tri_basis_1);
                float b = vec3::dot(vec, tri_basis_2);
                return -tri_basis_2*a + tri_basis_1*b;
            };

            vec3 shift = triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);
            // std::cout << triangle_projection_matrix[tri] << "\n";
            // std::cout << u << "\n";
            // printf("%.2g\n", m_current_time_step_dt);
            // std::cout << -m_current_time_step_dt*u << "\n";
            // std::cout << shift << "\n";getchar();
    
            if (vec3::dot(shift, perp(p2 - p1)) <= 0 &&
                    vec3::dot(shift, perp(p3 - p1)) >= 0) {
                double fix = 0.001;
                Face out_face;
                vec3 out_pos;
                std::tie(out_face, out_pos) = traverse(tri, p1+(1+fix)*shift, (1-fix)*shift);
                assert(!out_face.null());

                Edge edges[3] = {
                    out_face.halfedge().edge(),
                    out_face.halfedge().next().edge(),
                    out_face.halfedge().next().next().edge()
                };
                Vertex vs[3] = {
                    out_face.halfedge().vertex(),
                    out_face.halfedge().next().vertex(),
                    out_face.halfedge().next().next().vertex()
                };
                vec3 ps[3];
                for (int i = 0; i < 3; i++) ps[i] = eigen_to_vec3(geom.position[vs[i]]);
                
                // Compute barycentric coordinates.
                double x = barycentric_coeff(out_pos, ps[1], ps[2], n);
                double y = barycentric_coeff(out_pos, ps[2], ps[0], n);
                double z = barycentric_coeff(out_pos, ps[0], ps[1], n);
                double w = x+y+z;
                x /= w;
                y /= w;
                z /= w;
                if (_test_point_1_mode == 0) {
	            _test_point_1[v] = ps[0]*x + ps[1]*y + ps[2]*z;
                } else {
	            _test_point_1[v] = p1 + shift;
                }
	        _test_point_2[v] = eigen_to_vec3(geom.barycenter(out_face));
                vec3 val = vec3(0,0,0);

                P2Element elements[6] = {vs[0], vs[1], vs[2], edges[0], edges[1], edges[2]};
                for (int i = 0; i < 6; i++) {
                    val += u_basis[i](x,y,z) * to_vec3(velocity[elements[i]]);
                }
                new_velocity[v] = to_vec2(val);
                break;
            }
            he = he.twin().next();
        } while (he != start);
    }

    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            new_velocity[edge] = vec2(0,0);
            continue;
        }
        vec3 u = to_vec3(velocity_prev[edge]);
        vec3 n = normal[edge];

        vec3 edge_vec = eigen_to_vec3(geom.position[edge.b().vertex()] - geom.position[edge.a().vertex()]);
        vec3 midpoint = eigen_to_vec3(geom.midpoint(edge));
        Face tri;
        if (vec3::dot(n, vec3::cross(-u, edge_vec)) > 0) {
            tri = edge.b().face();
        } else {
            tri = edge.a().face();
        }

        auto v1 = tri.halfedge().vertex();
        auto v2 = tri.halfedge().next().vertex();
        auto v3 = tri.halfedge().next().next().vertex();
        vec3 p1 = eigen_to_vec3(geom.position[v1]);
        vec3 p2 = eigen_to_vec3(geom.position[v2]);
        vec3 p3 = eigen_to_vec3(geom.position[v3]);
        auto edge_110 = tri.halfedge().edge();
        auto edge_011 = tri.halfedge().next().edge();
        auto edge_101 = tri.halfedge().next().next().edge();
        vec3 tri_n = triangle_normal[tri];

        vec3 shift = triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);

	double fix = 0.001;
	Face out_face;
	vec3 out_pos;
	std::tie(out_face, out_pos) = traverse(tri, midpoint+(1+fix)*shift, (1-fix)*shift);
	assert(!out_face.null());

        Edge edges[3] = {
            out_face.halfedge().edge(),
            out_face.halfedge().next().edge(),
            out_face.halfedge().next().next().edge()
        };
        Vertex vs[3] = {
            out_face.halfedge().vertex(),
            out_face.halfedge().next().vertex(),
            out_face.halfedge().next().next().vertex()
        };
        vec3 ps[3];
        for (int i = 0; i < 3; i++) ps[i] = eigen_to_vec3(geom.position[vs[i]]);
        
        // Compute barycentric coordinates.
        double x = barycentric_coeff(out_pos, ps[1], ps[2], n);
        double y = barycentric_coeff(out_pos, ps[2], ps[0], n);
        double z = barycentric_coeff(out_pos, ps[0], ps[1], n);
        double w = x+y+z;
        x /= w;
        y /= w;
        z /= w;
        if (_test_point_1_mode == 0) {
	    _test_point_1[edge] = ps[0]*x + ps[1]*y + ps[2]*z;
        } else {
	    _test_point_1[edge] = midpoint + shift;
        }
        _test_point_2[edge] = eigen_to_vec3(geom.barycenter(out_face));
        vec3 val = vec3(0,0,0);

	P2Element elements[6] = {vs[0], vs[1], vs[2], edges[0], edges[1], edges[2]};
	for (int i = 0; i < 6; i++) {
	    val += u_basis[i](x,y,z) * to_vec3(velocity[elements[i]]);
	}
        new_velocity[edge] = to_vec2(val);
    }

    for (auto v : geom.mesh.vertices()) {
        velocity[v] = new_velocity[v];
    }
    for (auto edge : geom.mesh.edges()) {
        velocity[edge] = new_velocity[edge];
    }
}


void NavierStokesSolver::explicit_advection()
{
    if (m_advection_traversal) {
        explicit_advection_traversal();
    } else {
        explicit_advection_lagrangian();
    }
}
