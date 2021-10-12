
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

    void solve();
    
    // Solution (initializes to 0).
    //------------------------------------------------------------
    // P2 coefficients for velocity u.
    P2Attachment<vec2> u;
    // P1 coefficients for pressure p.
    VertexAttachment<double> p;

    // Dirichlet boundary condition.
    P2Attachment<vec2> u_boundary;
    // Set boundary condition from a function.
    void set_u_boundary(PlaneVectorField _u_boundary);

    // Additional mesh data.
    //------------------------------------------------------------
    // Flat index ordering of vertices and midpoints.
    VertexAttachment<int> vertex_indices;
    EdgeAttachment<int> midpoint_indices;
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
    midpoint_indices(_geom.mesh),
    midpoints(_geom.mesh),
    geom{_geom}
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    num_boundary_edges = 0;
    num_interior_edges = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
    }
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            midpoint_indices[edge] = -1;
            num_boundary_edges += 1;
        } else {
            midpoint_indices[edge] = num_interior_edges;
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
    
    // Misc. additions for debugging.
    write_sparsity_pattern = false;
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


void Solver::solve()
{
    // DEBUGGING FLAGS
    const bool BUILD_TOP_LEFT = true;
    const bool BUILD_TOP_RIGHT = true;
    const bool BUILD_BOTTOM_LEFT = false;


    // u is approximated by P2 vector elements, which are trivial products of scalar P2 elements.
    // p is approximated by P1 elements.

    // N_u: The number of vector coefficients of u.
    int N_u = num_interior_vertices + num_interior_edges;
    // N_p: The number of coefficients of p.
    int N_p = num_interior_vertices;

    // system_N: The size of the linear system (system_N x system_N matrix).
    int system_N = 2*N_u + N_p;

    // Linear system ordering:
    // - Firstly, 2*N_u basis functions for Phi^u, alternating between the x component and y component of
    //     each aggregate basis function phi^u.
    // - Secondly, N_p basis functions for Phi^p.

    // Initialize data used to construct the linear system.
    auto rhs = Eigen::VectorXd(system_N);
    for (int i = 0; i < system_N; i++) rhs[i] = 0.;
    std::vector<EigenTriplet> coefficients;

    auto add_entry = [&](int i, int j, double value) {
        // Add the value to entry (i,j) in the resulting matrix.
        printf("%d %d %.6f\n", i,j,value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };

if (BUILD_TOP_LEFT) {
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
        int v_index = vertex_indices[v]; // Global interior vertex index.
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
                int vp_index = vertex_indices[vp];
                insert_top_left_block(v_index, vp_index, val);
            }
            
            // vpp contribution.
            val = -(1./6.) * C * K2.dot(K3);
            if (vpp.on_boundary()) {
                vec2 bv = u_boundary[vpp];
                rhs[2*v_index+0] -= bv.x() * val;
                rhs[2*v_index+1] -= bv.y() * val;
            } else {
                int vpp_index = vertex_indices[vpp];
                insert_top_left_block(v_index, vpp_index, val);
            }

            // midpoint_vp contribution.
            val = (2./3.)*C*K1.dot(K3);
            int midpoint_vp_index = midpoint_indices[vp_edge];
            insert_top_left_block(v_index, num_interior_vertices + midpoint_vp_index, val);
            
            // midpoint_vpp contribution.
            val = (2./3.)*C*K2.dot(K3);
            int midpoint_vpp_index = midpoint_indices[vpp_edge];
            insert_top_left_block(v_index, num_interior_vertices + midpoint_vpp_index, val);

            he = he.twin().next();
        } while (he != start);
    }
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = midpoint_indices[edge];
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
            int v_index = vertex_indices[v];
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();
            auto midpoint_vp_index = midpoint_indices[midpoint_vp];
            auto midpoint_vpp_index = midpoint_indices[midpoint_vpp];
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
} // end if (BUILD_TOP_LEFT)
if (BUILD_TOP_RIGHT) {
    // Construct the top-right block, consisting of 2x1 vectors.
    //================================================================================
    auto insert_top_right_block = [&](int psi_u_index, int phi_p_index, double val_x, double val_y) {
        add_entry(2*psi_u_index+0, 2*N_u + phi_p_index, val_x);
        add_entry(2*psi_u_index+1, 2*N_u + phi_p_index, val_y);

        //------------------------------------------------------------testing
        add_entry(2*N_u + phi_p_index, 2*psi_u_index+0, val_x);
        add_entry(2*N_u + phi_p_index, 2*psi_u_index+1, val_y);
    };
    // For each basis trial function psi^u ...
    //------------------------------------------------------------
    // For each psi^u based on a vertex.
    //     All integrals happen to be zero.
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = midpoint_indices[edge];
        int rhs_edge_index = num_interior_vertices + midpoint_index; // index into the RHS.
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
            int v_index = vertex_indices[v];
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();
            auto midpoint_vp_index = midpoint_indices[midpoint_vp];
            auto midpoint_vpp_index = midpoint_indices[midpoint_vpp];
            auto midpoint_vp_pos = midpoints[midpoint_vp];
            auto midpoint_vpp_pos = midpoints[midpoint_vpp];
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            const double coeff = 1./6.; // Integral of quadratic phi_110 on the reference triangle.
            double R = coeff * 0.5/geom.triangle_area(tri);

            double val_x = 0.;
            double val_y = 0.;
            // Integrate psi^u at edge against phi^p at v.
            val_x = R * K1.dot(K3);
            val_y = R * K2.dot(K3);
            if (v.on_boundary()) {
                vec2 bv = u_boundary[v];
                rhs[2*rhs_edge_index+0] -= bv.x() * val_x;
                rhs[2*rhs_edge_index+1] -= bv.y() * val_y;
            } else {
                insert_top_right_block(rhs_edge_index, v_index, val_x, val_y);
            }
            
            // Integrate psi^u at edge against phi^p at vp.
            val_x = R * K1.dot(K1);
            val_y = R * K2.dot(K1);
            if (vp.on_boundary()) {
                vec2 bv = u_boundary[vp];
                rhs[2*rhs_edge_index+0] -= bv.x() * val_x;
                rhs[2*rhs_edge_index+1] -= bv.y() * val_y;
            } else {
                insert_top_right_block(rhs_edge_index, vp_index, val_x, val_y);
            }
            
            // Integrate psi^u at edge against phi^p at vpp.
            val_x = R * K1.dot(K2);
            val_y = R * K2.dot(K2);
            if (vpp.on_boundary()) {
                vec2 bv = u_boundary[vpp];
                rhs[2*rhs_edge_index+0] -= bv.x() * val_x;
                rhs[2*rhs_edge_index+1] -= bv.y() * val_y;
            } else {
                insert_top_right_block(rhs_edge_index, vpp_index, val_x, val_y);
            }
        }
    }
} // end if (BUILD_TOP_RIGHT)
if (BUILD_BOTTOM_LEFT) {
    // Construct the bottom-left block, consisting of 1x2 vectors.
    //================================================================================
    auto insert_bottom_left_block = [&](int psi_p_index, int phi_u_index, double val_x, double val_y) {
        add_entry(2*N_u + psi_p_index, 2*phi_u_index + 0, val_x);
        add_entry(2*N_u + psi_p_index, 2*phi_u_index + 1, val_x);
    };
    // For each basis trial function psi^p ...
    //------------------------------------------------------------
    // The psi^p are only based on vertices.
    for (auto v : geom.mesh.vertices()) {
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
            Edge vp_edge = he.edge(); // contains midpoint_vp
            Edge opposite_edge = he.next().edge(); // contains opposite_midpoint
            Edge vpp_edge = he.next().next().edge(); // contains midpoint_vpp
            int vp_edge_index = midpoint_indices[vp_edge];
            int opposite_edge_index = midpoint_indices[opposite_edge];
            int vpp_edge_index = midpoint_indices[vpp_edge];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            const double coeff = 1./6.; // Integral of quadratic phi_110, phi_011, phi_101 on the reference triangle.
                                        // (integral of phi_002, phi_020, phi_200 is zero).
            double R = coeff * 0.5/geom.triangle_area(tri);

            double val_x = 0.;
            double val_y = 0.;
            // Integrate psi^p at v against phi^u at vp_edge.
            val_x = R * K1.dot(K3);
            val_y = R * K2.dot(K3);
            if (vp_edge.on_boundary()) {
                vec2 bv = u_boundary[vp_edge];
                rhs[2*N_u + v_index] -= bv.x() * val_x;
                rhs[2*N_u + v_index] -= bv.y() * val_y;
            } else {
                insert_bottom_left_block(v_index, num_interior_vertices + vp_edge_index, val_x, val_y);
            }
            // Integrate psi^p at v against phi^u at opposite_edge.
            if (opposite_edge.on_boundary()) {
                vec2 bv = u_boundary[opposite_edge];
                rhs[2*N_u + v_index] -= bv.x() * val_x;
                rhs[2*N_u + v_index] -= bv.y() * val_y;
            } else {
                insert_bottom_left_block(v_index, num_interior_vertices + opposite_edge_index, val_x, val_y);
            }
            // Integrate psi^p at v against phi^u at vpp_edge..
            if (vpp_edge.on_boundary()) {
                vec2 bv = u_boundary[vpp_edge];
                rhs[2*N_u + v_index] -= bv.x() * val_x;
                rhs[2*N_u + v_index] -= bv.y() * val_y;
            } else {
                insert_bottom_left_block(v_index, num_interior_vertices + vpp_edge_index, val_x, val_y);
            }

            he = he.twin().next();
        } while (he != start);
    }
} // end if (BUILD_BOTTOM_LEFT)
    

    // Finalize the mass matrix.
    auto mass_matrix = SparseMatrix(system_N, system_N);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();

    // std::cout << "mass matrix:\n" << Eigen::MatrixXd(mass_matrix) << "\n";
    // printf("rhs: ");
    // for (int i = 0; i < system_N; i++) printf("%.5g, ", rhs[i]);
    // printf("\n");
    
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
        fprintf(ppm_file, "P3\n");
        fprintf(ppm_file, "# %d vertices, %d triangles, %dx%d system with %d entries, %d non-zeros, fill %.4f\n",
            geom.mesh.num_vertices(),
            geom.mesh.num_faces(),
            system_N, system_N,
            system_N * system_N,
            num_nonzeros,
            num_nonzeros * 1.f/(system_N * system_N)
        );
        fprintf(ppm_file, "%d %d\n", system_N, system_N);
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
            fprintf(ppm_file, "\n");
        }
        fclose(ppm_file);
        exit(EXIT_SUCCESS);
    }

    /*--------------------------------------------------------------------------------
        Solve the system.
    --------------------------------------------------------------------------------*/
    //--------------------------------------------------------------------------------
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    //--------------------------------------------------------------------------------
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(mass_matrix);
    // Compute the numerical factorization 
    solver.factorize(mass_matrix);
    // Use the factors to solve the linear system 
    Eigen::VectorXd up = solver.solve(rhs);
    
    /*--------------------------------------------------------------------------------
    // Reassociate each coefficient (or boundary value) with the corresponding vertex or edge of the mesh.
    --------------------------------------------------------------------------------*/
    // Velocity
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            u[v] = u_boundary[v];
            p[v] = 0.; // Pressure boundary condition.
        } else {
            u[v] = vec2(up[2*interior_vertex_index+0],
		        up[2*interior_vertex_index+1]);
            p[v] = up[2*N_u + interior_vertex_index];
            interior_vertex_index += 1;
        }
        // std::cout << u[v] << "\n";
    }
    int interior_midpoint_index = 0;
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            u[edge] = u_boundary[edge];
        } else {
            u[edge] = vec2(up[2*(num_interior_vertices + interior_midpoint_index) + 0],
                           up[2*(num_interior_vertices + interior_midpoint_index) + 1]);
            interior_midpoint_index += 1;
        }
        // std::cout << u[edge] << "\n";
    }
}
