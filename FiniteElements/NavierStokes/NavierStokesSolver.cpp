/*================================================================================
    P2-P1 Taylor-Hood finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/


/*--------------------------------------------------------------------------------
    Initialize the solver.
--------------------------------------------------------------------------------*/
NavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity) :
    geom{_geom},
    velocity(_geom),
    pressure(_geom),
    velocity_prev(_geom),
    pressure_prev(_geom),
    source_samples(_geom),
    velocity_node_indices(_geom),
    pressure_node_indices(_geom),
{
    m_solving = false;
    m_iterating = false;
    m_kinematic_viscosity = _kinematic_viscosity;
    m_time = 0.;

    m_num_velocity_variation_nodes = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    m_num_pressure_variation_nodes = geom.mesh.num_vertices();
    m_system_N = 2*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes;
    
    /*--------------------------------------------------------------------------------
        Set up the indexing for nodes in the mesh.
        This gives a canonical ordering, which determines the ordering of the residual.
    --------------------------------------------------------------------------------*/
    // Velocity node ordering. Interior vertex nodes, then interior edge nodes.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) velocity_node_indices[v] = -1;
        velocity_node_indices[v] = counter;
        counter += 1;
    }
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) velocity_node_indices[e] = -1;
        velocity_node_indices[e] = counter;
        counter += 1;
    }
    // Pressure node ordering. Vertex nodes (interior or on the boundary).
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        pressure_node_indices[v] = counter;
        counter += 1;
    }

    // Initialize the pressure and velocity fields to 0.
    for (auto v : geom.mesh.vertices()) {
        velocity[v] = vec2(0.,0.);
        velocity_prev[v] = vec2(0.,0.);
        pressure[v] = 0.;
        pressure_prev[v] = 0.;
    }
    for (auto e : geom.mesh.edges()) {
        velocity[e] = vec2(0.,0.);
        velocity_prev[e] = vec2(0.,0.);
    }
    m_velocity_pressure_vector = Eigen::VectorXd(m_system_N);
    for (int i = 0; i < m_system_N; i++) m_velocity_pressure_vector[i] = 0.;
}

/*--------------------------------------------------------------------------------
    Set up the solver and model parameters.
--------------------------------------------------------------------------------*/
void NavierStokesSolver::set_source(PlaneVectorField vf)
{
    assert(!solving());
    source_function = vf;
}


/*--------------------------------------------------------------------------------
    Begin solving and iterating.
--------------------------------------------------------------------------------*/
void NavierStokesSolver::start_time_step(double dt)
{
    assert(!iterating());
    if (!solving()) {
        m_solving = true;
    }

    m_current_time_step_dt = dt;
    m_iterating = true;
}

void NavierStokesSolver::newton_iteration()
{
    assert(solving());

    Eigen::VectorXd residual = compute_residual();
    SparseMatrix J = compute_gateaux_matrix();
    
    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double> > linear_solver;
    linear_solver.compute(J);
    Eigen::VectorXd velocity_pressure_variation_vector = solver.solve(residual);

    velocity_pressure_vector -= velocity_pressure_variation_vector;
}


void NavierStokesSolver::end_time_step()
{
    assert(solving());
    assert(iterating());

    m_time += m_current_time_step_dt;
    m_iterating = false;
}

/*--------------------------------------------------------------------------------
    Matrix assembly.
--------------------------------------------------------------------------------*/
Eigen::VectorXd NavierStokesSolver::compute_residual()
{
    auto velocity_residual = P2Attachment<vec2>(geom);
    auto pressure_residual = P1Attachment<double>(geom);

    /*--------------------------------------------------------------------------------
        Time-step update term.
    dot((u - u_prev)/dt, psi^u)
    --------------------------------------------------------------------------------*/
    add_velocity_residual_time_step(velocity_residual);
    
    /*--------------------------------------------------------------------------------
        Source term.
    dot(-source_function, psi^u)
    --------------------------------------------------------------------------------*/
    add_velocity_residual_source(velocity_residual);

    /*--------------------------------------------------------------------------------
        Pressure term.
    -p * div(psi^u)
    --------------------------------------------------------------------------------*/
    add_velocity_residual_pressure(velocity_residual);

    /*--------------------------------------------------------------------------------
        Viscosity term.
    kinematic_viscosity * grad(u):grad(psi^u)
    --------------------------------------------------------------------------------*/
    add_velocity_residual_viscosity(velocity_residual);

    /*--------------------------------------------------------------------------------
        Advection term.
    dot(dot(u, grad(u)), psi^u)
    --------------------------------------------------------------------------------*/
    add_velocity_residual_advection(velocity_residual);


    // Put these residual attachments into a length (2*N_u + N_p) vector.
    auto residual = Eigen::VectorXd(m_system_N);;
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        residual[2*counter+0] = velocity_residual[v].x();
        residual[2*counter+1] = velocity_residual[v].y();
        counter += 1;
    }
    counter = 0;
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) continue;
        residual[2*geom.mesh.num_interior_vertices() + 2*counter+0] = velocity_residual[e].x();
        residual[2*geom.mesh.num_interior_vertices() + 2*counter+1] = velocity_residual[e].y();
        counter += 1;
    }
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        residual[2*num_velocity_variation_nodes() + counter] = pressure_residual[v];
        counter += 1;
    }
    return residual;
}

SparseMatrix NavierStokesSolver::compute_gateaux_matrix()
{
    struct TopLeftEntry {
        Element velocity_trial_node;
        int trial_component; // x: 0, y: 1
        Element velocity_test_node;
        int test_component;  // x: 0, y: 1
        double value;
    };
    struct BottomLeftEntry {
        Vertex pressure_trial_node;
        Element velocity_test_node;
        int test_component;  // x: 0, y: 1
        double value;
    };
    auto top_left_coefficients = std::vector<TopLeftEntry>();
    auto bottom_left_coefficients = std::vector<BottomLeftEntry>();

    // Construct the sparse matrix by converting the coefficient lists (which are in terms of the mesh)
    // into a list of triplets indexing into a matrix.
    auto eigen_coefficients = std::vector<EigenTriplet>();
    for (auto coeff : top_left_coefficients) {
        eigen_coefficients.push_back(EigenTriplet(
            2*P2_indices[coeff.velocity_trial_node] + coeff.trial_component,
            2*P2_indices[coeff.velocity_test_node] + coeff.test_component,
            coeff.value
        ));
    }
    for (auto coeff : bottom_left_coefficients) { // The matrix should be symmetric, so also set the top-right block.
        eigen_coefficients.push_back(EigenTriplet(
            2*m_num_velocity_variation_nodes + P1_indices[coeff.pressure_trial_node],
            2*P2_indices[coeff.velocity_test_node] + coeff.test_component,
            coeff.value
        ));
    }
    // Convert the list of triplets into a sparse matrix.
    auto gateaux_matrix = SparseMatrix(m_system_N, m_system_N);
    gateaux_matrix.setFromTriplets(eigen_coefficients.begin(), eigen_coefficients.end());
    gateaux_matrix.makeCompressed();
    return gateaux_matrix;
}
