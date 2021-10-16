

#include "WeaklyIncompressible/sparsity_pattern.cpp"




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
    
    // Divergence elimination parameter.
    double C;

    bool solving; // When one iteration is made, some options (like setting the boundary) are locked.
    void iterate(); // Make one WI algorithm iteration.
    void project_divergence(); // Once the fixed point iteration is stable, project out a velocity component to make the divergence zero.

    void velocity_laplacian_system(SparseMatrix &mass_matrix, Eigen::VectorXd &rhs);
    Eigen::VectorXd pressure_gradient_source();
    SparseMatrix compute_pressure_gramian_matrix();
    SparseMatrix compute_scalar_velocity_gramian_matrix();
    SparseMatrix pressure_gramian_matrix;
    SparseMatrix scalar_velocity_gramian_matrix;
    void pressure_update(bool dont_actually_update = false);
    SparseMatrix velocity_laplacian_matrix;
    Eigen::VectorXd velocity_laplacian_rhs;
    
    // Poisson subproblem (different boundary conditions)
    void scalar_poisson_system(SparseMatrix &mass_matrix, Eigen::VectorXd &rhs, PlaneFunction source, PlaneFunction dirichlet_boundary_function);
    
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
    // Extra solution data that might be useful.
    // P1 coefficients for divergence of velocity, div(u).
    VertexAttachment<double> div_u;

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
    EdgeAttachment<int> midpoint_indices;
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

    //================================================================================
    // Chorin projection
    // void scalar_poisson_chorin(SparseMatrix &matrix, Eigen::VectorXd &rhs, VertexAttachment<double> source);

    SparseMatrix gramian_matrix_P2();
    void div_P2_P2(P2Attachment<vec2> &vf, P2Attachment<double> &div); // Compute the divergence of vf in P2_2 projected into P2.
};

Solver::Solver(SurfaceGeometry &_geom, double _mu) :
    mu{_mu},
    u(_geom.mesh),
    p(_geom.mesh),
    div_u(_geom.mesh),
    u_boundary(_geom.mesh),
    vertex_indices(_geom.mesh),
    midpoint_indices(_geom.mesh),
    interior_vertex_indices(_geom.mesh),
    interior_midpoint_indices(_geom.mesh),
    midpoints(_geom.mesh),
    geom{_geom}
{
    solving = false;

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
    int num_midpoints = 0;
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            interior_midpoint_indices[edge] = -1;
            num_boundary_edges += 1;
        } else {
            interior_midpoint_indices[edge] = num_interior_edges;
            num_interior_edges += 1;
        }
        midpoint_indices[edge] = num_midpoints;
        num_midpoints += 1;
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

    // Default divergence elimination parameter
    C = 0.001;

    // set_pressure([](double x,double y) { return 100; });

    pressure_gramian_matrix = compute_pressure_gramian_matrix();
    scalar_velocity_gramian_matrix = compute_scalar_velocity_gramian_matrix();
}


void Solver::set_u_boundary(PlaneVectorField vf)
{
    assert(!solving);
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





#include "WeaklyIncompressible/pressure_gradient_source.cpp"
#include "WeaklyIncompressible/pressure_update.cpp"
#include "WeaklyIncompressible/velocity_laplacian_system.cpp"
#include "WeaklyIncompressible/scalar_poisson_system.cpp"
#include "WeaklyIncompressible/project_divergence.cpp"

#include "WeaklyIncompressible/div_P2_P2.cpp"


// void Solver::scalar_poisson_chorin(SparseMatrix &matrix, Eigen::VectorXd &rhs, P2Attachment)
// {
//     
// }


#if 0
void Solver::iterate()
{
    // Compute the Laplacian matrix and boundary terms.
    velocity_laplacian_system(velocity_laplacian_matrix, velocity_laplacian_rhs);

    for (auto v : geom.mesh.vertices()) {
        p[v] = 0.;
    }
    // Compute the pressure gradient source term.
    Eigen::VectorXd rhs = velocity_laplacian_rhs + pressure_gradient_source();

    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(velocity_laplacian_matrix);
    solver.factorize(velocity_laplacian_matrix);
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
    pressure_update(true); // compute div(u)

    // Solve -Laplacian(gamma) = -
    SparseMatrix gamma_matrix;
    Eigen::VectorXd gamma_rhs;
    scalar_poisson_chorin(gamma_matrix, gamma_rhs, div_u);

}
#else
void Solver::iterate()
{
    if (!solving) {
        // Initialize.
        
        // Compute the Laplacian matrix and boundary terms.
        velocity_laplacian_system(velocity_laplacian_matrix, velocity_laplacian_rhs);
        
        // // Initialize the velocity.
        for (auto v : geom.mesh.vertices()) {
            if (v.on_boundary()) u[v] = u_boundary[v];
            else u[v] = vec2(0,0);
        }
        for (auto e : geom.mesh.edges()) {
            if (e.on_boundary()) u[e] = u_boundary[e];
            else u[e] = vec2(0,0);
        }

        // Initialize the pressure.
        // pressure_update(false);
        for (auto v : geom.mesh.vertices()) {
            // p[v] = 0.;
            // p[v] = frand();
            // if (v.on_boundary()) {
            //     p[v] = -2;
            // } else {
            //     p[v] = 2;
            // }
        }
        pressure_update(true); // compute div(u)

        solving = true;
        // return;
    }
    // Compute the pressure gradient source term.
    Eigen::VectorXd rhs = velocity_laplacian_rhs + pressure_gradient_source();

    /*--------------------------------------------------------------------------------
        Solve the system for u_{n+1}.
    --------------------------------------------------------------------------------*/
    //--------------------------------------------------------------------------------
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    //--------------------------------------------------------------------------------
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(velocity_laplacian_matrix);
    // Compute the numerical factorization 
    solver.factorize(velocity_laplacian_matrix);
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
    pressure_update(!solving); // don't actually update the pressure for the first iteration.
    solving = true;
}
#endif

