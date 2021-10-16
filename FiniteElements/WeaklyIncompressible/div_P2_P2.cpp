

// Gramian projection matrix for P2 (no boundary condition).
SparseMatrix Solver::gramian_matrix_P2()
{
    std::vector<EigenTriplet> coefficients;
    int N = geom.mesh.num_vertices() + geom.mesh.num_edges();
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^us on a vertex.
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
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            int edge_110_index = midpoint_indices[edge_110];
            int edge_011_index = midpoint_indices[edge_011];
            int edge_101_index = midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

            add_entry(v_index, v_index, (1./60.)*R);
            add_entry(v_index, vp_index, (-1./360.)*R);
            add_entry(v_index, vpp_index, (-1./360.)*R);

            add_entry(v_index, geom.mesh.num_vertices() + edge_110_index, (-1./90.)*R);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }
    // For each psi^us at a midpoint.
    for (auto edge : geom.mesh.edges()) {
        int midpoint_index = midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            if (tri.null()) continue;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            int v_index = vertex_indices[v];
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            int edge_110_index = midpoint_indices[edge_110];
            int edge_011_index = midpoint_indices[edge_011];
            int edge_101_index = midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

	    add_entry(geom.mesh.num_vertices() + edge_110_index,
		      geom.mesh.num_vertices() + edge_110_index, (4./45.)*R);
            add_entry(geom.mesh.num_vertices() + edge_110_index,
                      geom.mesh.num_vertices() + edge_011_index, (2./45.)*R);
            
            add_entry(geom.mesh.num_vertices() + edge_110_index,
                      geom.mesh.num_vertices() + edge_101_index, (2./45.)*R);
            
            add_entry(geom.mesh.num_vertices() + edge_110_index,
                      v_index, (-1./90.)*R);
            
        }
    }
    

    auto matrix = SparseMatrix(N, N);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    make_sparsity_image(matrix, DATA "P2_P2_gramian.ppm");

    return matrix;
}

// Gramian projection matrix for P2_0 (zero on the boundary).
SparseMatrix Solver::gramian_matrix_P2_0()
{
    std::vector<EigenTriplet> coefficients;
    int N = num_interior_vertices + num_interior_edges;
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^us on a vertex.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto start = v.halfedge();
        auto he = start;
        int v_index = interior_vertex_indices[v];
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            int edge_110_index = interior_midpoint_indices[edge_110];
            int edge_011_index = interior_midpoint_indices[edge_011];
            int edge_101_index = interior_midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

            add_entry(v_index, v_index, (1./60.)*R);
            if (!vp.on_boundary()) add_entry(v_index, vp_index, (-1./360.)*R);
            if (!vpp.on_boundary()) add_entry(v_index, vpp_index, (-1./360.)*R);

            if (!edge_110.on_boundary()) add_entry(v_index, num_interior_vertices + edge_110_index, (-1./90.)*R);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }
    // For each psi^us at a midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = interior_midpoint_indices[edge];
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
            int v_index = interior_vertex_indices[v];
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            int edge_110_index = interior_midpoint_indices[edge_110];
            int edge_011_index = interior_midpoint_indices[edge_011];
            int edge_101_index = interior_midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

	    if (!edge_110.on_boundary()) add_entry(num_interior_vertices + edge_110_index,
		      num_interior_vertices + edge_110_index, (4./45.)*R);
            if (!edge_011.on_boundary()) add_entry(num_interior_vertices + edge_110_index,
                      num_interior_vertices + edge_011_index, (2./45.)*R);
            
            if (!edge_101.on_boundary()) add_entry(num_interior_vertices + edge_110_index,
                      num_interior_vertices + edge_101_index, (2./45.)*R);
            
            if (!v.on_boundary()) add_entry(num_interior_vertices + edge_110_index,
                      v_index, (-1./90.)*R);
            
        }
    }
    

    auto matrix = SparseMatrix(N, N);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    make_sparsity_image(matrix, DATA "P20_P20_gramian.ppm");

    return matrix;
}


// Compute the divergence of vf in P2_2 projected into P2.
void Solver::div_P2_P2(P2Attachment<vec2> &vf, P2Attachment<double> &div)
{
    int N = num_interior_vertices + num_interior_edges;
    auto div_vector_proj = Eigen::VectorXd(N);
    for (int i = 0; i < N; i++) div_vector_proj[i] = 0.;
    
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

            vec2 vf_002 = vf[v];
            vec2 vf_200 = vf[vp];
            vec2 vf_020 = vf[vpp];
            vec2 vf_101 = vf[edge_101];
            vec2 vf_011 = vf[edge_011];
            vec2 vf_110 = vf[edge_110];

            // integral += (1./15.)*vec2::dot(K3.perp(), vf_002);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_020);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_200);

            // integral += (1./10.)*vec2::dot(K3.perp(), vf_011);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_110);
            // integral += (1./10.)*vec2::dot(K3.perp(), vf_101);
            
            integral += (1./15.)*vec2::dot(K3.perp(), vf_002);
            integral += (-1./30.)*vec2::dot(K2.perp(), vf_020);
            integral += (-1./30.)*vec2::dot(K1.perp(), vf_200);

            integral += vec2::dot((1./30.)*K1.perp() + (1./10.)*K2.perp(), vf_011);
            integral += vec2::dot((-1./30.)*K1.perp() + (-1./30.)*K2.perp(), vf_110);
            integral += vec2::dot((1./10.)*K1.perp() + (1./30.)*K2.perp(), vf_101);
            

            he = he.twin().next();
        } while (!he.face().null() && he != start);

        div_vector_proj[v_index] -= integral;
        
    }

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
            if (tri.null()) continue;
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

            vec2 vf_002 = vf[v];
            vec2 vf_200 = vf[vp];
            vec2 vf_020 = vf[vpp];
            vec2 vf_101 = vf[edge_101];
            vec2 vf_011 = vf[edge_011];
            vec2 vf_110 = vf[edge_110];
            // integral += (1./30.)*vec2::dot(K3.perp(), vf_002);
            // integral += vec2::dot((1./15.)*K1.perp() + (-1./30.)*K2.perp(), vf_020);
            // integral += vec2::dot((-1./30.)*K1.perp() + (1./15.)*K2.perp(), vf_200);

            // integral += vec2::dot((4./15.)*K1.perp() + (2./15.)*K2.perp(), vf_011);
            // integral += vec2::dot((4./15.)*K1.perp() + (4./15.)*K2.perp(), vf_110);
            // integral += vec2::dot((2./15.)*K1.perp() + (4./15.)*K2.perp(), vf_101);
            integral += vec2::dot((-1./30.)*K3.perp(), vf_002);
            integral += vec2::dot((1./10.)*K2.perp(), vf_020);
            integral += vec2::dot((1./10.)*K1.perp(), vf_200);

            integral += vec2::dot((-4./15.)*K1.perp() + (-2./15.)*K2.perp(), vf_011);
            integral += vec2::dot((4./15.)*K1.perp() + (4./15.)*K2.perp(), vf_110);
            integral += vec2::dot((-2./15.)*K1.perp() + (-4./15.)*K2.perp(), vf_101);
        }
        div_vector_proj[num_interior_vertices + midpoint_index] -= integral;
    }
    
    SparseMatrix gramian = gramian_matrix_P2_0();
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(gramian);
    solver.factorize(gramian);
    Eigen::VectorXd div_vector = solver.solve(div_vector_proj);

    // std::cout << "div_vector\n";
    // std::cout << div_vector << "\n";
    // std::cout << "div_vector_proj\n";
    // std::cout << div_vector_proj << "\n";
    // getchar();

    // Associate this divergence to the P2 mesh.
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (!v.on_boundary()) {
            div[v] = div_vector[interior_vertex_index];
            interior_vertex_index += 1;
        } else {
            div[v] = 0.;
        }
    }
    int interior_edge_index = 0;
    for (auto e : geom.mesh.edges()) {
        if (!e.on_boundary()) {
            div[e] = div_vector[num_interior_vertices + interior_edge_index];
            interior_edge_index += 1;
        } else {
            div[e] = 0.;
        }
    }
}

