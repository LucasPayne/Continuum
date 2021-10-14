
// A P2 triangle mesh discretization has nodal points at vertices and edge midpoints.
// Functions on these nodes are aggregated into a P2Attachment.
template <typename T>
class P2Attachment {
public:
    P2Attachment(SurfaceMesh &_mesh) :
        vertex_attachment(_mesh),
        edge_attachment(_mesh),
        mesh{_mesh}
    {}
    inline T operator[](Vertex v) const {
        return vertex_attachment[v];
    }
    inline T operator[](Edge e) const {
        return edge_attachment[e];
    }
    inline T &operator[](Vertex v) {
        return vertex_attachment[v];
    }
    inline T &operator[](Edge e) {
        return edge_attachment[e];
    }
    VertexAttachment<T> vertex_attachment;
    EdgeAttachment<T> edge_attachment;
private:
    SurfaceMesh &mesh;
};


struct Solver {
    Solver(SurfaceGeometry &_geom, double _mu);
    // Equation parameters.
    double mu;

    void iterate(); // Make one WI algorithm iteration.
    void velocity_laplacian_system(SparseMatrix &mass_matrix, Eigen::VectorXd &rhs);
    Eigen::VectorXd pressure_gradient_source();
    void pressure_update();
    
    int wi_iteration_number; // Start at n=0.
    // N_u: The number of vector coefficients of u.
    int N_u;
    // system_N: The size of the linear system (system_N x system_N matrix).
    int system_N;
    // Solution p_n, u_n (initializes to p0 = 0, u0 = 0).
    //------------------------------------------------------------
    // P2 coefficients for velocity u.
    P2Attachment<vec2> u;
    // P1 coefficients for pressure p.
    VertexAttachment<double> p;

    // Dirichlet boundary condition.
    P2Attachment<vec2> u_boundary;
    // Set boundary condition from a function.
    void set_u_boundary(PlaneVectorField _u_boundary);
    // Set pressure from a function. NOTE: This is only for testing.
    void set_pressure(PlaneFunction _pressure);

    // Additional mesh data.
    //------------------------------------------------------------
    // Flat index ordering of vertices and midpoints.
    VertexAttachment<int> vertex_indices;
    VertexAttachment<int> interior_vertex_indices;
    EdgeAttachment<int> interior_midpoint_indices;
    // Store precomputed midpoints for convenience.
    EdgeAttachment<Eigen::Vector3f> midpoints;

    // Mesh properties.
    int num_boundary_vertices;
    int num_interior_vertices;
    int num_boundary_edges;
    int num_interior_edges;

    SurfaceGeometry &geom;

    // Misc. additions for debugging.
    bool write_sparsity_pattern;
};

Solver::Solver(SurfaceGeometry &_geom, double _mu) :
    mu{_mu},
    u(_geom.mesh),
    p(_geom.mesh),
    u_boundary(_geom.mesh),
    vertex_indices(_geom.mesh),
    interior_vertex_indices(_geom.mesh),
    interior_midpoint_indices(_geom.mesh),
    midpoints(_geom.mesh),
    geom{_geom}
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    num_boundary_edges = 0;
    num_interior_edges = 0;
    int num_vertices = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            interior_vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            interior_vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
        vertex_indices[v] = num_vertices;
        num_vertices += 1;
    }
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            interior_midpoint_indices[edge] = -1;
            num_boundary_edges += 1;
        } else {
            interior_midpoint_indices[edge] = num_interior_edges;
            num_interior_edges += 1;
        }
    }
    // Compute midpoints.
    for (auto edge : geom.mesh.edges()) {
        auto p = geom.position[edge.a().vertex()];
        auto pp = geom.position[edge.b().vertex()];
        midpoints[edge] = 0.5*p + 0.5*pp;
    }

    // Zero-initialize the boundary function.
    set_u_boundary([](double,double) { return vec2(0.,0.); });
    // Initialize WI algorithm.
    wi_iteration_number = 0; // wi stands for "weakly incompressible".
    // Initialize pressure p0 to 0.
    set_pressure([](double,double) { return 0.; });
    // u0 doesn't matter, just set it to zero.
    for (auto v : geom.mesh.vertices()) {
        u[v] = vec2(0,0);
    }
    for (auto e : geom.mesh.edges()) {
        u[e] = vec2(0,0);
    }
    // N_u: The number of vector coefficients of u.
    N_u = num_interior_vertices + num_interior_edges;
    // system_N: The size of the linear system (system_N x system_N matrix).
    system_N = 2*N_u;
    
    // Misc. additions for debugging.
    write_sparsity_pattern = false;



    set_pressure([](double x,double y) { return 0.5*(x+1); });
}


void Solver::set_u_boundary(PlaneVectorField vf)
{
    for (auto v : geom.mesh.vertices()) {
        auto pos = geom.position[v];
        u_boundary[v] = vf(pos.x(), pos.z());
    }
    for (auto e : geom.mesh.edges()) {
        auto pos = midpoints[e];
        u_boundary[e] = vf(pos.x(), pos.z());
    }
}

// Testing function
void Solver::set_pressure(PlaneFunction _pressure)
{
    for (auto v : geom.mesh.vertices()) {
        auto pos = geom.position[v];
        p[v] = _pressure(pos.x(), pos.z());
    }
}




Eigen::VectorXd Solver::pressure_gradient_source()
{
    auto source = Eigen::VectorXd(2*N_u);
    for (int i = 0; i < 2*N_u; i++) source[i] = 0.;

    // Trial psi^u at vertices.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

        double integral_x = 0.;
        double integral_y = 0.;
        // For each triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double R = -0.5/geom.triangle_area(tri);

            vec2 val_v = vec2(0.,0.);

            // Integrate psi_u at v against p[v]*phi_p at v.
            // zero
            
            // val_v = p[v] * R * ((-1./6.)*K1 + (-1./6.)*K2); ---incorrect
            // integral_x += vec2::dot(val_v, K1);
            // integral_y += vec2::dot(val_v, K2);
            
            // Integrate psi_u at v against phi_p at vp.
            // zero
            
            // Integrate psi_u at v against phi_p at vpp.
            // zero
        
            he = he.twin().next();
        } while (he != start);

        int v_index = interior_vertex_indices[v];
        source[2*v_index + 0] = integral_x;
        source[2*v_index + 1] = integral_y;
    }
    // Trial psi^u at midpoints.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        double integral_x = 0.;
        double integral_y = 0.;
        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            printf("t = %d\n------------------------------------------------------------\n", t);
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.tip();
            auto vpp = he.vertex();

            int v_index = interior_vertex_indices[v];
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double R = -0.5/geom.triangle_area(tri);

            printf("pressure: %.6f %.6f %.6f\n", p[v], p[vp], p[vpp]);

            vec2 val_v = vec2(0.,0.);
            // Integrate psi_u at edge_110 against phi_p at 002 (v).
            val_v = p[v] * R * ((-1./6.)*K1 + (-1./6.)*K2);
            std::cout << "phi_p v " << val_v << "\n";
            printf("x, y %.6f %.6f\n", vec2::dot(val_v, K1), vec2::dot(val_v, K2));
            integral_x += vec2::dot(val_v, K1);
            integral_y += vec2::dot(val_v, K2);
            printf("integral: %.6f %.6f\n", integral_x, integral_y);
            getchar();

            // Integrate psi_u at edge_110 against phi_p at 200 (vp).
            val_v = p[vp] * R * ((1./6.)*K1);
            std::cout << "phi_p vp " << val_v << "\n";
            printf("x, y %.6f %.6f\n", vec2::dot(val_v, K1), vec2::dot(val_v, K2));
            integral_x += vec2::dot(val_v, K1);
            integral_y += vec2::dot(val_v, K2);
            printf("integral: %.6f %.6f\n", integral_x, integral_y);
            getchar();
            
            // Integrate psi_u at edge_110 against phi_p at 002 (vpp).
            val_v = p[vpp] * R * ((1./6.)*K2);
            std::cout << "phi_p_vpp" << val_v << "\n";
            printf("x, y %.6f %.6f\n", vec2::dot(val_v, K1), vec2::dot(val_v, K2));
            integral_x += vec2::dot(val_v, K1);
            integral_y += vec2::dot(val_v, K2);
            printf("integral: %.6f %.6f\n", integral_x, integral_y);
            getchar();
        }
        int index = num_interior_vertices + interior_midpoint_indices[edge];
        printf("integral x,y %.6f %.6f\n", integral_x, integral_y);
        source[2*index + 0] = integral_x;
        source[2*index + 1] = integral_y;
    }
    std::cout << source << "\n"; getchar();

    return source;
}

void Solver::pressure_update()
{

}



void Solver::iterate()
{
    // Compute the Laplacian matrix and boundary terms.
    SparseMatrix mass_matrix;
    Eigen::VectorXd rhs;
    velocity_laplacian_system(mass_matrix, rhs);
    // Compute the pressure gradient source term.
    rhs += pressure_gradient_source();

    /*--------------------------------------------------------------------------------
        Solve the system for u_{n+1}.
    --------------------------------------------------------------------------------*/
    //--------------------------------------------------------------------------------
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    //--------------------------------------------------------------------------------
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(mass_matrix);
    // Compute the numerical factorization 
    solver.factorize(mass_matrix);
    // Use the factors to solve the linear system 
    Eigen::VectorXd u_vector = solver.solve(rhs);
    /*--------------------------------------------------------------------------------
    // Reassociate each velocity coefficient (or boundary value) with the corresponding vertex or edge of the mesh.
    --------------------------------------------------------------------------------*/
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            u[v] = u_boundary[v];
        } else {
            u[v] = vec2(u_vector[2*interior_vertex_index+0],
		        u_vector[2*interior_vertex_index+1]);
            interior_vertex_index += 1;
        }
    }
    int interior_midpoint_index = 0;
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            u[edge] = u_boundary[edge];
        } else {
            u[edge] = vec2(u_vector[2*(num_interior_vertices + interior_midpoint_index) + 0],
                           u_vector[2*(num_interior_vertices + interior_midpoint_index) + 1]);
            interior_midpoint_index += 1;
        }
    }

    // Update the pressure.
    pressure_update();
}


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
        fprintf(ppm_file, "# %d vertices, %d triangles, %dx%d system with %d entries, %d non-zeros, fill %.4f\n",
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
