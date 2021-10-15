

void Solver::project_divergence()
{
    // Compute alpha, the average divergence.
    double w = 1. / geom.mesh.num_vertices();
    double alpha = 0.;
    for (auto v : geom.mesh.vertices()) {
        alpha += div_u[v] * w;
    }

    // Solve -laplacian(gamma) = alpha, zero on the boundary.
    // The result is gamma in P2 space.
    SparseMatrix matrix;
    Eigen::VectorXd rhs;
    scalar_poisson_system(matrix, rhs,
                          [&](double, double) { return alpha; },
                          [&](double, double) { return 0.; });

    // printf("%.6g\n", alpha);
    // getchar();

    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(matrix);
    solver.factorize(matrix);
    Eigen::VectorXd gamma_vector = solver.solve(rhs);
    // Reassociate gamma to the mesh.
    P2Attachment<double> gamma(geom.mesh);
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            gamma[v] = 0;
        } else {
            gamma[v] = gamma_vector[interior_vertex_index];
            interior_vertex_index += 1;
        }
    }
    int interior_midpoint_index = 0;
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            gamma[edge] = 0;
        } else {
            gamma[edge] = gamma_vector[num_interior_vertices + interior_midpoint_index];
            interior_midpoint_index += 1;
        }
    }
    
    
    // Project grad(gamma) to P2.
    int N = num_interior_vertices + num_interior_edges;
    auto grad_gamma_proj = std::vector<vec2>(N);
    for (int i = 0; i < N; i++) grad_gamma_proj[i] = vec2(0,0);

    // For each psi_u scalar at a vertex.
    for (Vertex v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = interior_vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom.position[v];

        vec2 integral = vec2(0.,0.);

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
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            double gamma_002 = gamma[v];
            double gamma_200 = gamma[vp];
            double gamma_020 = gamma[vpp];
            double gamma_101 = gamma[edge_101];
            double gamma_011 = gamma[edge_011];
            double gamma_110 = gamma[edge_110];
            integral += (1./15.)*K3.perp() * gamma_002;
            integral += (-1./30.)*K3.perp() * gamma_020;
            integral += (-1./30.)*K3.perp() * gamma_200;

            integral += (-1./10.)*K3.perp() * gamma_011;
            integral += (1./30.)*K3.perp() * gamma_110;
            integral += (-1./10.)*K3.perp() * gamma_101;

            he = he.twin().next();
        } while (he != start);

        grad_gamma_proj[v_index] += integral;
        
    }
    // For each psi_u scalar at a midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = interior_midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};

        vec2 integral = vec2(0,0);
        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            double gamma_002 = gamma[v];
            double gamma_200 = gamma[vp];
            double gamma_020 = gamma[vpp];
            double gamma_101 = gamma[edge_101];
            double gamma_011 = gamma[edge_011];
            double gamma_110 = gamma[edge_110];

            integral += (1./30.)*K3.perp() * gamma_002;
            integral += ((1./15.)*K1.perp() + (-1./30.)*K2.perp()) * gamma_020;
            integral += ((-1./30.)*K1.perp() + (1./15.)*K2.perp()) * gamma_200;

            integral += ((4./15.)*K1.perp() + (2./15.)*K2.perp()) * gamma_011;
            integral += ((4./15.)*K1.perp() + (4./15.)*K2.perp()) * gamma_110;
            integral += ((2./15.)*K1.perp() + (4./15.)*K2.perp()) * gamma_101;
        }
        grad_gamma_proj[num_interior_vertices + midpoint_index] += integral;
    }

}

