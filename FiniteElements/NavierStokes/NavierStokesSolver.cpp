/*================================================================================
    P2-P1 Taylor-Hood finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/
#include "NavierStokes/NavierStokesSolver.h"
#include "NavierStokes/core.h"

// utilities
vec3 eigen_to_vec3(Eigen::Vector3f v)
{
    return vec3(v.x(), v.y(), v.z());
}
Eigen::Vector3f vec3_to_eigen(vec3 v)
{
    return Eigen::Vector3f(v.x(), v.y(), v.z());
}


/*--------------------------------------------------------------------------------
    Initialize the solver.
--------------------------------------------------------------------------------*/
NavierStokesSolver::NavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity) :
    geom{_geom},
    velocity(_geom.mesh),
    pressure(_geom.mesh),
    velocity_prev(_geom.mesh),
    pressure_prev(_geom.mesh),
    velocity_node_indices(_geom.mesh),
    pressure_node_indices(_geom.mesh),
    source_samples_P2(_geom.mesh)
{
    m_solving = false;
    m_iterating = false;
    m_kinematic_viscosity = _kinematic_viscosity;
    m_time = 0.;

    m_num_velocity_variation_nodes = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    m_num_pressure_variation_nodes = geom.mesh.num_vertices();
    // The size of the system is 2*N_u + N_p - 1.
    // The -1 is due to one pressure node being fixed.
    // (As a convention, the fixed node is the last in the ordering.)
    m_system_N = 2*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes - 1;
    
    /*--------------------------------------------------------------------------------
        Set up the indexing for nodes in the mesh.
        This gives a canonical ordering, which determines the ordering of the residual.
    --------------------------------------------------------------------------------*/
    // Velocity node ordering. Interior vertex nodes, then interior edge nodes.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            velocity_node_indices[v] = -1;
        } else {
            velocity_node_indices[v] = counter;
            counter += 1;
        }
    }
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) {
            velocity_node_indices[e] = -1;
        } else {
            velocity_node_indices[e] = counter;
            counter += 1;
        }
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
}

/*--------------------------------------------------------------------------------
    Set up the solver and model parameters.
--------------------------------------------------------------------------------*/
void NavierStokesSolver::set_source(PlaneVectorField vf)
{
    assert(!solving());
    source_function = vf;
    for (auto v : geom.mesh.vertices()) {
        auto pos = geom.position[v];
        source_samples_P2[v] = vf(pos.x(), pos.z());
    }
    for (auto e : geom.mesh.edges()) {
        auto pos = 0.5*geom.position[e.a().vertex()] + 0.5*geom.position[e.b().vertex()];
        source_samples_P2[e] = vf(pos.x(), pos.z());
    }
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

    // Initialize the previous state.
    // (The start of the Newton iteration is at the previous state.)
    for (auto v : geom.mesh.vertices()) {
        velocity_prev[v] = velocity[v];
        pressure_prev[v] = pressure[v];
    }
    for (auto e : geom.mesh.edges()) {
        velocity_prev[e] = velocity[e];
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
    Eigen::VectorXd velocity_pressure_variation_vector = linear_solver.solve(residual);
    
    std::cout << velocity_pressure_variation_vector << "\n";getchar();

    // Update the current velocity and pressure.
    // The vector of variations has to be reassociated to the mesh.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        velocity[v].x() -= velocity_pressure_variation_vector[2*counter+0];
        velocity[v].y() -= velocity_pressure_variation_vector[2*counter+1];
        counter += 1;
    }
    counter = 0;
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) continue;
        velocity[e].x() -= velocity_pressure_variation_vector[2*geom.mesh.num_interior_vertices() + 2*counter+0];
        velocity[e].y() -= velocity_pressure_variation_vector[2*geom.mesh.num_interior_vertices() + 2*counter+1];
        counter += 1;
    }
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (counter == m_num_pressure_variation_nodes-1) break; // skip the last pressure node
        pressure[v] -= velocity_pressure_variation_vector[2*num_velocity_variation_nodes() + counter];
        counter += 1;
    }
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
    auto velocity_residual = P2Attachment<vec2>(geom.mesh);
    auto pressure_residual = P2Attachment<double>(geom.mesh);
    compute_velocity_residual(velocity_residual);
    compute_pressure_residual(pressure_residual);

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
        if (counter == m_num_pressure_variation_nodes-1) break; // skip the last pressure node
        residual[2*num_velocity_variation_nodes() + counter] = pressure_residual[v];
        counter += 1;
    }
    return residual;
}



void make_sparsity_image(SparseMatrix &matrix, std::string name)
{
    assert(matrix.rows() == matrix.cols());
    int system_N = matrix.rows();
    int num_nonzeros = 0;
    for (int i = 0; i < system_N; i++) {
        for (int j = 0; j < system_N; j++) {
            if (fabs(matrix.coeff(i, j)) >= 1e-4) {
                num_nonzeros += 1;
            }
        }
    }
    FILE *ppm_file = fopen(const_cast<const char *>(name.c_str()), "w+");
    fprintf(ppm_file, "P3\n");
    fprintf(ppm_file, "%d %d\n", system_N, system_N);
    fprintf(ppm_file, "255\n");
    for (int i = 0; i < system_N; i++) {
        for (int j = 0; j < system_N; j++) {
            if (fabs(matrix.coeff(i, j)) >= 1e-4) {
                if (fabs(matrix.coeff(i, j) - matrix.coeff(j, i)) <= 1e-4) {
                    // signify when this entry is symmetric (equal to its corresponding transpose entry).
                    fprintf(ppm_file, "255 0 0 ");
                } else {
                    fprintf(ppm_file, "0 0 0 ");
                }
            } else {
	        fprintf(ppm_file, "255 255 255 ");
            }
        }
    }
    fclose(ppm_file);
}

SparseMatrix NavierStokesSolver::compute_gateaux_matrix()
{
    // auto top_left_coefficients = std::vector<TopLeftEntry>();
    // auto bottom_left_coefficients = std::vector<BottomLeftEntry>();
    auto top_left_coefficients = compute_gateaux_matrix_top_left();
    auto bottom_left_coefficients = compute_gateaux_matrix_bottom_left();


    // Construct the sparse matrix by converting the coefficient lists (which are in terms of the mesh)
    // into a list of triplets indexing into a matrix.
    auto eigen_coefficients = std::vector<EigenTriplet>();
    auto add_eigen_coefficient = [&](int i, int j, double val) {
        printf("%d %d %.6g\n", i, j, val);
        eigen_coefficients.push_back(EigenTriplet(i, j, val));
    };
    for (auto coeff : top_left_coefficients) {
        // printf("trial: %d\n", velocity_node_indices[coeff.velocity_trial_node]);
        add_eigen_coefficient(
            2*velocity_node_indices[coeff.velocity_trial_node] + coeff.trial_component,
            2*velocity_node_indices[coeff.velocity_test_node] + coeff.test_component,
            coeff.value
        );
    }
    for (auto coeff : bottom_left_coefficients) {
        add_eigen_coefficient(
            2*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
            2*velocity_node_indices[coeff.velocity_test_node] + coeff.test_component,
            coeff.value
        );
        // The matrix should be symmetric, so also set the top-right block.
        add_eigen_coefficient(
            2*velocity_node_indices[coeff.velocity_test_node] + coeff.test_component,
            2*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
            coeff.value
        );
    }
    printf("N_u: %d\n", m_num_velocity_variation_nodes);
    printf("N_p: %d\n", m_num_pressure_variation_nodes);
    printf("N: %d\n", m_system_N);
    // getchar();

    // Convert the list of triplets into a sparse matrix.
    auto gateaux_matrix = SparseMatrix(m_system_N, m_system_N);
    gateaux_matrix.setFromTriplets(eigen_coefficients.begin(), eigen_coefficients.end());
    gateaux_matrix.makeCompressed();

    make_sparsity_image(gateaux_matrix, std::string(DATA) + "navier_stokes_gateaux_sparsity.png");

    return gateaux_matrix;
}
