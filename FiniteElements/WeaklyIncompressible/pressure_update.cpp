
// Form the Gramian matrix of integrals of phi_p_i*psi_p_j.
SparseMatrix Solver::compute_pressure_gramian_matrix()
{
    std::vector<EigenTriplet> coefficients;
    int N_p = geom.mesh.num_vertices();

    auto add_entry = [&](int i, int j, double value) {
        // printf("%d %d %.6f\n", i, j, value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^p.
    for (auto v : geom.mesh.vertices()) {
        auto start = v.halfedge();
        auto he = start;
        int v_index = vertex_indices[v];
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];

            double R = 2*geom.triangle_area(tri);

            add_entry(v_index, v_index, (1./12.)*R);
            add_entry(v_index, vp_index, (1./24.)*R);
            add_entry(v_index, vpp_index, (1./24.)*R);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    auto matrix = SparseMatrix(N_p, N_p);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    make_sparsity_image(matrix, DATA "pressure_gramian.ppm");

    return matrix;
}


void Solver::pressure_update()
{
    int N_p = geom.mesh.num_vertices();

    // While updating the pressure, the L2 projection of div(u) onto pressure space
    // is computed. Store it in this vector (so div(u) can be reconstructed for debugging and visualization later).
    auto u_div_l2_proj = Eigen::VectorXd(N_p);
    for (int i = 0; i < N_p; i++) u_div_l2_proj[i] = 0.;

    // Collect the current (to be previous) pressure into a vector.
    auto p_prev_vector = Eigen::VectorXd(N_p);
    int index = 0;
    for (auto v : geom.mesh.vertices()) {
        p_prev_vector[index] = p[v];
        index += 1;
    }
    Eigen::VectorXd rhs = pressure_gramian_matrix * p_prev_vector;
    // Compute the divergence term.
    // For each psi^p.
    for (auto v : geom.mesh.vertices()) {
        int v_index = vertex_indices[v];
        auto v_pos = geom.position[v];

        double integral = 0.;

        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            
            vec2 u110 = u[edge_110];
            vec2 u011 = u[edge_011];
            vec2 u101 = u[edge_101];
            integral += (1./6.) * vec2::dot(K3.perp(), u110 + u011 + u101);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);

        u_div_l2_proj[v_index] = integral;
        // u_div_l2_proj[v_index] = 0.5;

        rhs[v_index] += C*integral;
    }

    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(pressure_gramian_matrix);
    solver.factorize(pressure_gramian_matrix);
    Eigen::VectorXd p_vector = solver.solve(rhs);

    // Reassociate each p value to the mesh.
    index = 0;
    for (auto v : geom.mesh.vertices()) {
        p[v] = p_vector[index];
        index += 1;
    }

    // Reconstruct the P1 divergence of u, then associate it to the mesh.
    std::cout << Eigen::MatrixXd(pressure_gramian_matrix) << "\n";
    // Eigen::VectorXd div_u_vector = solver.solve(u_div_l2_proj);
    Eigen::VectorXd div_u_vector = u_div_l2_proj;
    index = 0;
    for (auto v : geom.mesh.vertices()) {
        div_u[v] = div_u_vector[index];
        index += 1;
    }
    for (auto v : geom.mesh.vertices()) {
        printf("%.6g\n", div_u[v]);
        std::cout << u[v] << "\n";
    }
    getchar();
}
