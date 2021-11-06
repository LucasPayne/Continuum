/*================================================================================
    P2-P1 Taylor-Hood mixed finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/
#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"
#include "core.h"


/*--------------------------------------------------------------------------------
    Initialize the solver.
--------------------------------------------------------------------------------*/
SurfaceNavierStokesSolver::SurfaceNavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity) :
    geom{_geom},

    velocity(_geom.mesh),
    pressure(_geom.mesh),
    centripetal(_geom.mesh),

    triangle_normal(_geom.mesh),
    triangle_projection_matrix(_geom.mesh),
    source_samples_P2(_geom.mesh),

    velocity_prev(_geom.mesh),
    pressure_prev(_geom.mesh),
    centripetal_prev(_geom.mesh),

    velocity_node_indices(_geom.mesh),
    pressure_node_indices(_geom.mesh)
{
    m_solving = false;
    m_kinematic_viscosity = _kinematic_viscosity;
    m_time = 0.;

    m_num_velocity_variation_nodes = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    m_num_pressure_variation_nodes = geom.mesh.num_vertices();
    m_num_centripetal_variation_nodes = geom.mesh.num_vertices();
    // The size of the system is 2*N_u + N_p - 1 + N_r.
    // The -1 is due to one pressure node being fixed.
    // (As a convention, the fixed node is the last in the ordering.)
    m_system_N = 2*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes - 1 + m_num_centripetal_variation_nodes;
    
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
    // These are also the centripetal node indices.
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        pressure_node_indices[v] = counter;
        counter += 1;
    }

    // Initialize the fields to 0.
    for (auto v : geom.mesh.vertices()) {
        velocity[v] = vec3(0.,0.,0.);
        velocity_prev[v] = vec3(0.,0.,0.);
        pressure[v] = 0.;
        pressure_prev[v] = 0.;
    }
    for (auto e : geom.mesh.edges()) {
        velocity[e] = vec3(0.,0.,0.);
        velocity_prev[e] = vec3(0.,0.,0.);
    }
    m_solution_vector = Eigen::VectorXd(m_system_N);
    for (int i = 0; i < m_system_N; i++) m_solution_vector[i] = 0.;

    // Set up surface normal data.
    for (auto tri : geom.mesh.faces()) {
        vec3 n = eigen_to_vec3(geom.triangle_normal(tri));
        triangle_normal[tri] = n;
        triangle_projection_matrix[tri] = mat3x3::identity() - vec3::outer(n, n);
    }

}

/*--------------------------------------------------------------------------------
    Solving.
--------------------------------------------------------------------------------*/
void SurfaceNavierStokesSolver::time_step(double delta_time)
{
    m_current_time_step_dt = delta_time;

    if (!solving()) {
        m_solving = true;
    }
    // Explicit advection.
    // explicit_advection_lagrangian();

    // Initialize the previous state.
    for (auto v : geom.mesh.vertices()) {
        velocity_prev[v] = velocity[v];
        pressure_prev[v] = pressure[v];
        centripetal_prev[v] = centripetal[v];
    }
    for (auto e : geom.mesh.edges()) {
        velocity_prev[e] = velocity[e];
    }
    
    // Initialize the solution vector.
    // Velocity.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        m_solution_vector[2*counter+0] = velocity[v].x();
        m_solution_vector[2*counter+1] = velocity[v].y();
        counter += 1;
    }
    counter = 0;
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) continue;
        m_solution_vector[2*geom.mesh.num_interior_vertices() + 2*counter+0] = velocity[e].x();
        m_solution_vector[2*geom.mesh.num_interior_vertices() + 2*counter+1] = velocity[e].y();
        counter += 1;
    }
    // Pressure.
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (counter == m_num_pressure_variation_nodes-1) break; // skip the last pressure node
        m_solution_vector[2*num_velocity_variation_nodes() + counter] = pressure[v];
        counter += 1;
    }
    // Centripetal.
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        m_solution_vector[2*num_velocity_variation_nodes() + num_pressure_variation_nodes()-1 + counter] = centripetal[v];
        counter += 1;
    }
    SparseMatrix matrix = compute_matrix();
    make_sparsity_image(matrix, DATA "upr_matrix.ppm");

    m_time += m_current_time_step_dt;
}




SparseMatrix SurfaceNavierStokesSolver::compute_matrix()
{
    auto velocity_block_coefficients = compute_velocity_block_coefficients();
    auto pressure_block_coefficients = compute_pressure_block_coefficients();
    auto centripetal_block_coefficients = compute_centripetal_block_coefficients();

    // Construct the sparse matrix by converting the coefficient lists (which are in terms of the mesh)
    // into a list of triplets indexing into a matrix.
    auto eigen_coefficients = std::vector<EigenTriplet>();
    auto add_eigen_coefficient = [&](int i, int j, double val) {
        // printf("%d %d, %.6g\n", i,j, val);
        eigen_coefficients.push_back(EigenTriplet(i, j, val));
    };
    for (auto coeff : velocity_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                add_eigen_coefficient(
                    2*velocity_node_indices[coeff.velocity_trial_node] + i,
                    2*velocity_node_indices[coeff.velocity_test_node] + j,
                    coeff.value.entry(i,j)
                );
            }
        }
    }
    getchar();
    for (auto coeff : pressure_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            add_eigen_coefficient(
                2*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
                2*velocity_node_indices[coeff.velocity_test_node] + i,
                coeff.value[i]
            );
            add_eigen_coefficient(
                2*velocity_node_indices[coeff.velocity_test_node] + i,
                2*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
                coeff.value[i]
            );
        }
    }
    getchar();
    for (auto coeff : centripetal_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            add_eigen_coefficient(
                2*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1 + pressure_node_indices[coeff.centripetal_trial_node],
                2*velocity_node_indices[coeff.velocity_test_node] + i,
                coeff.value[i]
            );
            add_eigen_coefficient(
                2*velocity_node_indices[coeff.velocity_test_node] + i,
                2*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1 + pressure_node_indices[coeff.centripetal_trial_node],
                coeff.value[i]
            );
        }
    }
    getchar();


    auto matrix = SparseMatrix(m_system_N, m_system_N);
    matrix.setFromTriplets(eigen_coefficients.begin(), eigen_coefficients.end());
    matrix.makeCompressed();
    return matrix;
}
