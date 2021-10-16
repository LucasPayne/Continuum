
void Solver::scalar_poisson_system(SparseMatrix &mass_matrix, Eigen::VectorXd &rhs, PlaneFunction source, PlaneFunction dirichlet_boundary_function)
{
    // The indices of interior vertices and interior midpoints are concatenated.
    //     note: When accessing a midpoint, the index must be shifted to start at num_interior_vertices.
    int N = num_interior_vertices + num_interior_edges;
    rhs = Eigen::VectorXd(N);
    for (int i = 0; i < N; i++) rhs[i] = 0.;
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    /*--------------------------------------------------------------------------------
        Compute the integrals for interior vertices.
    --------------------------------------------------------------------------------*/
    for (Vertex v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = interior_vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom.position[v];

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Edge vp_edge = he.edge(); // contains midpoint_vp
            Vertex vpp = he.next().tip();
            Edge vpp_edge = he.next().next().edge(); // contains midpoint_vpp
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            double C = 1.0/(4.0*geom.triangle_area(tri));
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            double val = 0.;
            // Diagonal term.
            val = C * 0.5 * K3.dot(K3);
            add_entry(v_index, v_index, val);

            // vp contribution.
            val = -(1./6.) * C * K1.dot(K3);
            if (vp.on_boundary()) {
                double bv = dirichlet_boundary_function(vp_pos.x(), vp_pos.z());
                rhs[v_index] -= bv * val;
            } else {
                int vp_index = interior_vertex_indices[vp];
                add_entry(v_index, vp_index, val);
            }
            
            // vpp contribution.
            val = -(1./6.) * C * K2.dot(K3);
            if (vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(vpp_pos.x(), vpp_pos.z());
                rhs[v_index] -= bv * val;
            } else {
                int vpp_index = interior_vertex_indices[vpp];
                add_entry(v_index, vpp_index, val);
            }

            // midpoint_vp contribution.
            val = (2./3.)*C*K1.dot(K3);
            int midpoint_vp_index = interior_midpoint_indices[vp_edge];
            add_entry(v_index, num_interior_vertices + midpoint_vp_index, val);
            
            // midpoint_vpp contribution.
            val = (2./3.)*C*K2.dot(K3);
            int midpoint_vpp_index = interior_midpoint_indices[vpp_edge];
            add_entry(v_index, num_interior_vertices + midpoint_vpp_index, val);

            he = he.twin().next();
        } while (he != start);
    }

    /*--------------------------------------------------------------------------------
        Compute the integrals for interior midpoints.
    --------------------------------------------------------------------------------*/
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = interior_midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};
        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            double C = 1.0/(4.0*geom.triangle_area(tri));
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            int v_index = interior_vertex_indices[v];
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();
            auto midpoint_vp_index = interior_midpoint_indices[midpoint_vp];
            auto midpoint_vpp_index = interior_midpoint_indices[midpoint_vpp];
            auto midpoint_vp_pos = midpoints[midpoint_vp];
            auto midpoint_vpp_pos = midpoints[midpoint_vpp];

            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            // printf("self midpoint\n");
            double val = 0.;
            // 110, 110
            val = 4.*C/3. * (K3.dot(K3) - K1.dot(K2));
            // val = 4.*C/3. * (K1.dot(K1) + K1.dot(K2) + K2.dot(K2));
            add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_index, val);

            // printf("boundary midpoints\n");
            // 110, 011
            val = 4.*C/3. * (K1.dot(K3));
            if (midpoint_vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(midpoint_vpp_pos.x(), midpoint_vpp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vpp_index, val);
            }
            
            // 110, 101
            val = 4.*C/3. * (K2.dot(K3));
            if (midpoint_vp.on_boundary()) {
                double bv = dirichlet_boundary_function(midpoint_vp_pos.x(), midpoint_vp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vp_index, val);
            }
            
            // printf("vertices\n");
            // 110, 200
            val = 2.*C/3. * (K1.dot(K2));
            if (vp.on_boundary()) {
                double bv = dirichlet_boundary_function(vp_pos.x(), vp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, vp_index, val);
            }
            // 110, 020
            if (vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(vpp_pos.x(), vpp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, vpp_index, val);
            }
        }
    }
    /*--------------------------------------------------------------------------------
        Source term integration.
    --------------------------------------------------------------------------------*/
    // Sample the source function at each vertex and midpoint.
    // This determines a piecewise quadratic interpolation.
    VertexAttachment<double> g_vertex(geom.mesh);
    EdgeAttachment<double> g_edge(geom.mesh);
    for (auto v : geom.mesh.vertices()) {
        auto v_pos = geom.position[v];
        g_vertex[v] = source(v_pos.x(), v_pos.z());
    }
    for (auto edge : geom.mesh.edges()) {
        auto midpoint_pos = midpoints[edge];
        g_edge[edge] = source(midpoint_pos.x(), midpoint_pos.z());
    }

    // Vertex trial integrals.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = interior_vertex_indices[v];
        auto v_pos = geom.position[v];

        double integral = 0.;

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            double tri_area = geom.triangle_area(tri);
            double J = tri_area;
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];

            auto opposite_edge = he.next().edge();
            auto opposite_edge_pos = midpoints[he.next().edge()];

            // // v
            // // integral += (1./60.) * J * g_vertex[v];
            // integral += (1./60.) * J * source(v_pos.x(), v_pos.z());
            // // vp
            // integral += (-1./360.) * J * source(vp_pos.x(), vp_pos.z());
            // // vpp
            // integral += (-1./360.) * J * source(vpp_pos.x(), vpp_pos.z());
            // // opposite_edge midpoint
            // integral += (-1./90.) * J * source(opposite_edge_pos.x(), opposite_edge_pos.z());
            // v
            integral += (1./60.) * J * g_vertex[v];
            // vp
            integral += (-1./360.) * J * g_vertex[vp];
            // vpp
            integral += (-1./360.) * J * g_vertex[vpp];
            // opposite_edge midpoint
            integral += (-1./90.) * J * g_edge[opposite_edge];

            he = he.next();
        } while (he != start);
        rhs[v_index] += integral;
    }

    // Midpoint trial integrals.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = interior_midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};

        double integral = 0.;

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            double tri_area = geom.triangle_area(tri);
            double J = tri_area;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();

            auto edge_pos = midpoints[edge];
            auto midpoint_vp_pos = midpoints[midpoint_vp];
            auto midpoint_vpp_pos = midpoints[midpoint_vpp];
            
            // // edge
            // integral += (4./45.) * J * source(edge_pos.x(), 0, edge_pos.z());
            // // midpoint_vp
            // integral += (2./45.) * J * source(midpoint_vp_pos.x(), 0, midpoint_vp_pos.z());
            // // midpoint_vpp
            // integral += (2./45.) * J * source(midpoint_vpp_pos.x(), 0, midpoint_vpp_pos.z());
            // // v
            // integral += (-1./90.) * J * source(v_pos.x(), v_pos.z());
            // edge
            integral += (4./45.) * J * g_edge[edge];
            // midpoint_vp
            integral += (2./45.) * J * g_edge[midpoint_vp];
            // midpoint_vpp
            integral += (2./45.) * J * g_edge[midpoint_vpp];
            // v
            integral += (-1./90.) * J * g_vertex[v];
        }
        rhs[num_interior_vertices + midpoint_index] += integral;
    }


    // Finalize the mass matrix.
    mass_matrix = SparseMatrix(N, N);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();
}


SparseMatrix Solver::laplacian_matrix_P1()
{
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    /*--------------------------------------------------------------------------------
        Construct the system.
    --------------------------------------------------------------------------------*/
    for (Vertex v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = vertex_indices[v]; // Global interior vertex index.
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
            double C = 1.0/(4.0*geom.triangle_area(tri));

            // Diagonal term.
            double diagonal_term = C*(vpp_pos - vp_pos).dot(vpp_pos - vp_pos);
            add_entry(v_index, v_index, diagonal_term);

            // vp contribution.
            if (vp.on_boundary()) {
            } else {
                double vp_term = C*(v_pos - vpp_pos).dot(vpp_pos - vp_pos);
                int vp_index = vertex_indices[vp];
                add_entry(v_index, vp_index, vp_term);
            }
            
            // vpp contribution.
            if (vpp.on_boundary()) {
            } else {
                double vpp_term = C*(vp_pos - v_pos).dot(vpp_pos - vp_pos);
                int vpp_index = vertex_indices[vpp];
                add_entry(v_index, vpp_index, vpp_term);
            }

            he = he.twin().next();
        } while (he != start);
    }
    auto matrix = SparseMatrix(num_interior_vertices, num_interior_vertices);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    return matrix;
}
