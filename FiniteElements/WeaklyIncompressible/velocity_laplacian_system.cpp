
void Solver::velocity_laplacian_system(SparseMatrix &mass_matrix, Eigen::VectorXd &rhs)
{
    // u is approximated by P2 vector elements, which are trivial products of scalar P2 elements.

    // Linear system ordering:
    // - 2*N_u basis functions for Phi^u, alternating between the x component and y component of
    //     each aggregate basis function phi^u.

    // Initialize data used to construct the linear system.
    rhs = Eigen::VectorXd(system_N);
    for (int i = 0; i < system_N; i++) rhs[i] = 0.;
    std::vector<EigenTriplet> coefficients;

    auto add_entry = [&](int i, int j, double value) {
        // Add the value to entry (i,j) in the resulting matrix.
        // printf("%d %d %.6f\n", i,j,value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };
    
    // Construct the top-left block, consisting of 2x2 multiples of the identity.
    //================================================================================
    auto insert_top_left_block = [&](int psi_u_index, int phi_u_index, double value) {
        add_entry(2*psi_u_index, 2*phi_u_index, value);
        add_entry(2*psi_u_index+1, 2*phi_u_index+1, value);
    };
    // For each basis trial function psi^u ...
    //------------------------------------------------------------
    // For each psi^u based on a vertex.
    for (auto v : geom.mesh.vertices()) {
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
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            // Note: Multiplication by mu is here.
            double C = mu * 1.0/(4.0*geom.triangle_area(tri));

            double val = 0.;
            // Diagonal term.
            val = C * 0.5 * K3.dot(K3);
            insert_top_left_block(v_index, v_index, val);

            // vp contribution.
            val = -(1./6.) * C * K1.dot(K3);
            if (vp.on_boundary()) {
                vec2 bv = u_boundary[vp];
                rhs[2*v_index+0] -= bv.x() * val;
                rhs[2*v_index+1] -= bv.y() * val;
            } else {
                int vp_index = interior_vertex_indices[vp];
                insert_top_left_block(v_index, vp_index, val);
            }
            
            // vpp contribution.
            val = -(1./6.) * C * K2.dot(K3);
            if (vpp.on_boundary()) {
                vec2 bv = u_boundary[vpp];
                rhs[2*v_index+0] -= bv.x() * val;
                rhs[2*v_index+1] -= bv.y() * val;
            } else {
                int vpp_index = interior_vertex_indices[vpp];
                insert_top_left_block(v_index, vpp_index, val);
            }

            // midpoint_vp contribution.
            val = (2./3.)*C*K1.dot(K3);
            int midpoint_vp_index = interior_midpoint_indices[vp_edge];
            insert_top_left_block(v_index, num_interior_vertices + midpoint_vp_index, val);
            
            // midpoint_vpp contribution.
            val = (2./3.)*C*K2.dot(K3);
            int midpoint_vpp_index = interior_midpoint_indices[vpp_edge];
            insert_top_left_block(v_index, num_interior_vertices + midpoint_vpp_index, val);

            he = he.twin().next();
        } while (he != start);
    }
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = interior_midpoint_indices[edge];
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

            // Note: Multiplication by mu is here.
            double C = mu * 1.0/(4.0*geom.triangle_area(tri));

            double val = 0.;
            // 110, 110
            val = 4.*C/3. * (K3.dot(K3) - K1.dot(K2));
            // val = 4.*C/3. * (K1.dot(K1) + K1.dot(K2) + K2.dot(K2));
            insert_top_left_block(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_index, val);

            // 110, 011
            val = 4.*C/3. * (K1.dot(K3));
            if (midpoint_vpp.on_boundary()) {
                vec2 bv = u_boundary[midpoint_vpp];
                rhs[2*(num_interior_vertices + midpoint_index)+0] -= bv.x() * val;
                rhs[2*(num_interior_vertices + midpoint_index)+1] -= bv.y() * val;
            } else {
                insert_top_left_block(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vpp_index, val);
            }
            
            // 110, 101
            val = 4.*C/3. * (K2.dot(K3));
            if (midpoint_vp.on_boundary()) {
                vec2 bv = u_boundary[midpoint_vp];
                rhs[2*(num_interior_vertices + midpoint_index)+0] -= bv.x() * val;
                rhs[2*(num_interior_vertices + midpoint_index)+1] -= bv.y() * val;
            } else {
                insert_top_left_block(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vp_index, val);
            }
            
            // 110, 200
            val = 2.*C/3. * (K1.dot(K2));
            if (vp.on_boundary()) {
                vec2 bv = u_boundary[vp];
                rhs[2*(num_interior_vertices + midpoint_index)+0] -= bv.x() * val;
                rhs[2*(num_interior_vertices + midpoint_index)+1] -= bv.y() * val;
            } else {
                insert_top_left_block(num_interior_vertices + midpoint_index, vp_index, val);
            }
            // 110, 020
            if (vpp.on_boundary()) {
                vec2 bv = u_boundary[vpp];
                rhs[2*(num_interior_vertices + midpoint_index)+0] -= bv.x() * val;
                rhs[2*(num_interior_vertices + midpoint_index)+1] -= bv.y() * val;
            } else {
                insert_top_left_block(num_interior_vertices + midpoint_index, vpp_index, val);
            }
        }
    }

    // Finalize the mass matrix.
    mass_matrix = SparseMatrix(system_N, system_N);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();

    // (figure creation)
    // Write the sparsity pattern to a PPM file.
    if (write_sparsity_pattern) {
        int num_nonzeros = 0;
        for (int i = 0; i < system_N; i++) {
            for (int j = 0; j < system_N; j++) {
                if (fabs(mass_matrix.coeff(i, j)) >= 1e-4) {
                    num_nonzeros += 1;
                }
            }
        }

        FILE *ppm_file = fopen(DATA "sparsity_pattern.ppm", "w+");
        const int rhs_pixel_width = std::max(3, system_N / 7);
        fprintf(ppm_file, "P3\n");
        fprintf(ppm_file, "# %zu vertices, %zu triangles, %dx%d system with %d entries, %d non-zeros, fill %.4f\n",
            geom.mesh.num_vertices(),
            geom.mesh.num_faces(),
            system_N, system_N,
            system_N * system_N,
            num_nonzeros,
            num_nonzeros * 1.f/(system_N * system_N)
        );
        fprintf(ppm_file, "%d %d\n", system_N + rhs_pixel_width, system_N);
        fprintf(ppm_file, "255\n");
        for (int i = 0; i < system_N; i++) {
            for (int j = 0; j < system_N; j++) {
                if (fabs(mass_matrix.coeff(i, j)) >= 1e-4) {
                    if (fabs(mass_matrix.coeff(i, j) - mass_matrix.coeff(j, i)) <= 1e-4) {
                        // signify when this entry is symmetric (equal to its corresponding transpose entry).
                        fprintf(ppm_file, "255 0 0 ");
                    } else {
                        fprintf(ppm_file, "0 0 0 ");
                    }
                } else {
                    if (i >= 2*N_u && j >= 2*N_u) {
                        // signify the lower-right zero block.
                        fprintf(ppm_file, "0 255 0 ");
                    } else if (i >= 2*N_u || j >= 2*N_u) {
                        // signify the top-right and bottom-left blocks.
                        fprintf(ppm_file, "120 0 230 ");
                    } else if (i >= 2*num_interior_vertices || j >= 2*num_interior_vertices) {
                        // signify the midpoints division in the top-left block.
                        fprintf(ppm_file, "245 245 245 ");
                    } else {
                        fprintf(ppm_file, "255 255 255 ");
                    }
                }
            }
            // Draw RHS.
            for (int j = system_N; j < system_N + rhs_pixel_width; j++) {
                if (fabs(rhs[i]) >= 1e-4) {
		    fprintf(ppm_file, "0 0 0 ");
                } else {
		    fprintf(ppm_file, "255 255 255 ");
                }
            }
            fprintf(ppm_file, "\n");
        }
        fclose(ppm_file);
        exit(EXIT_SUCCESS);
    }
}
