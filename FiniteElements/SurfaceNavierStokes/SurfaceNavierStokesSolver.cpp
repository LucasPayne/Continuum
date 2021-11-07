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
    normal(_geom.mesh),
    source_samples_P2(_geom.mesh),

    velocity_node_indices(_geom.mesh),
    pressure_node_indices(_geom.mesh)
{
    m_solving = false;
    m_kinematic_viscosity = _kinematic_viscosity;
    m_time = 0.;

    m_num_velocity_variation_nodes = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    m_num_pressure_variation_nodes = geom.mesh.num_vertices();
    m_num_centripetal_variation_nodes = geom.mesh.num_vertices();
    // The size of the system is 3*N_u + N_p - 1 + N_r.
    // The -1 is due to one pressure node being fixed.
    // (As a convention, the fixed node is the last in the ordering.)
    m_system_N = 3*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes - 1 + m_num_centripetal_variation_nodes;
    
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
        pressure[v] = 0.;
        centripetal[v] = 0.;
    }
    for (auto e : geom.mesh.edges()) {
        velocity[e] = vec3(0.,0.,0.);
    }
    m_solution_vector = Eigen::VectorXd(m_system_N);
    for (int i = 0; i < m_system_N; i++) m_solution_vector[i] = 0.;

    // Set up surface normal data.
    P2Attachment<int> num_tris(geom.mesh);
    for (auto v : geom.mesh.vertices()) {
        num_tris[v] = 0;
        normal[v] = vec3(0,0,0);
    }
    for (auto e : geom.mesh.edges()) {
        num_tris[e] = 0;
        normal[e] = vec3(0,0,0);
    }
    for (auto tri : geom.mesh.faces()) {
        vec3 n = eigen_to_vec3(geom.triangle_normal(tri));
        triangle_normal[tri] = n;
        triangle_projection_matrix[tri] = mat3x3::identity() - vec3::outer(n, n);
        auto start = tri.halfedge();
        auto he = start;
        do {
            normal[he.vertex()] += n;
            num_tris[he.vertex()] += 1;
            normal[he.edge()] += n;
            num_tris[he.edge()] += 1;
            he = he.next();
        } while (he != start);
    }
    for (auto v : geom.mesh.vertices()) {
        normal[v] /= num_tris[v];
    }
    for (auto e : geom.mesh.edges()) {
        normal[e] /= num_tris[e];
    }
    

    // Set up source data.
    for (auto v : geom.mesh.vertices()) {
        source_samples_P2[v] = vec3(0.,0.,0.);
    }
    for (auto e : geom.mesh.edges()) {
        source_samples_P2[e] = vec3(0.,0.,0.);
    }

}

/*--------------------------------------------------------------------------------
    Solving.
--------------------------------------------------------------------------------*/
void SurfaceNavierStokesSolver::time_step(double delta_time)
{
    if (!solving()) {
        m_solving = true;
    }
    m_current_time_step_dt = delta_time;

    // Explicit advection.
    // explicit_advection_lagrangian();

    printf("Constructing matrix...\n");
    SparseMatrix matrix = compute_matrix();
    make_sparsity_image(matrix, DATA "upr_matrix.ppm");
    printf("Constructing RHS...\n");
    Eigen::VectorXd rhs = compute_rhs();
    std::cout << rhs << "\n";

    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double> > linear_solver;
    printf("Factoring...\n");
    linear_solver.compute(matrix);
    printf("Solving...\n");
    m_solution_vector = linear_solver.solve(rhs);
    printf("Solved.\n");

    // Associate the solution to the mesh.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
    }
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        velocity[v].x() = m_solution_vector[3*counter+0];
        velocity[v].y() = m_solution_vector[3*counter+1];
        velocity[v].z() = m_solution_vector[3*counter+2];
        counter += 1;
    }
    counter = 0;
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) continue;
        velocity[e].x() = m_solution_vector[3*(geom.mesh.num_interior_vertices() + counter)+0];
        velocity[e].y() = m_solution_vector[3*(geom.mesh.num_interior_vertices() + counter)+1];
        velocity[e].z() = m_solution_vector[3*(geom.mesh.num_interior_vertices() + counter)+2];
        counter += 1;
    }
    // Pressure.
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (counter == m_num_pressure_variation_nodes-1) break; // skip the last pressure node
        pressure[v] = m_solution_vector[3*num_velocity_variation_nodes() + counter];
        counter += 1;
    }
    // Centripetal.
    counter = 0;
    for (auto v : geom.mesh.vertices()) {
        centripetal[v] = m_solution_vector[3*num_velocity_variation_nodes() + num_pressure_variation_nodes()-1 + counter];
        counter += 1;
    }
    
    m_time += m_current_time_step_dt;
}



void SurfaceNavierStokesSolver::set_source(std::function<vec3(double,double,double)> vf)
{
    for (auto v : geom.mesh.vertices()) {
        auto p = geom.position[v];
        source_samples_P2[v] = vf(p.x(), p.y(), p.z());
    }
    for (auto e : geom.mesh.edges()) {
        auto p = 0.5*geom.position[e.a().vertex()] + 0.5*geom.position[e.b().vertex()];
        source_samples_P2[e] = vf(p.x(), p.y(), p.z());
    }
}
