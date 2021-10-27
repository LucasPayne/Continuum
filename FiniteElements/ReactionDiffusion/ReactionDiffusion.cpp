/*================================================================================
================================================================================*/
#include "ReactionDiffusion/ReactionDiffusion.h"
#include "ReactionDiffusion/core.h"

// utilities
vec3 eigen_to_vec3(Eigen::Vector3f v)
{
    return vec3(v.x(), v.y(), v.z());
}
Eigen::Vector3f vec3_to_eigen(vec3 v)
{
    return Eigen::Vector3f(v.x(), v.y(), v.z());
}

ReactorDiffuser::ReactorDiffuser(SurfaceGeometry &_geom, float _mu, std::function<double(vec3, double, double)> _reaction_function) :
    geom{&_geom},
    mu{_mu},
    reaction_function{_reaction_function},
    interior_vertex_indices{_geom.mesh},
    u_mesh{_geom.mesh}
{
    time = 0.;
    num_nodes = geom->mesh.num_interior_vertices();
    u_vector = Eigen::VectorXd(num_nodes);
    for (int i = 0; i < num_nodes; i++) u_vector[i] = 0.;

    // Set up the node indices.
    int counter = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            interior_vertex_indices[v] = -1;
        } else {
            interior_vertex_indices[v] = counter;
            counter += 1;
        }
    }
    
    // Compute the Gramian matrix.
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        printf("%d %d %.6g\n", i,j,value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom->position[v];
        int v_index = interior_vertex_indices[v];

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];
	    int vp_index = interior_vertex_indices[vp];
	    int vpp_index = interior_vertex_indices[vpp];

            double tri_area = geom->triangle_area(tri);

            double val = 0.;

            val = 2*tri_area * (1./12.);
            add_entry(v_index, v_index, val);
            
            val = 2*tri_area * (1./24.);
            if (!vp.on_boundary()) add_entry(v_index, vp_index, val);
            if (!vpp.on_boundary()) add_entry(v_index, vpp_index, val);

            he = he.twin().next();
        } while (he != start);
    }
    gramian_matrix = SparseMatrix(num_nodes, num_nodes);
    gramian_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    // gramian_matrix.makeCompressed();
}


void ReactorDiffuser::time_step(double delta_time)
{
    int N = num_nodes;
    auto rhs = Eigen::VectorXd(N);
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };
    double inv_delta_time = 1./delta_time;

    /*--------------------------------------------------------------------------------
        Construct the system.
    --------------------------------------------------------------------------------*/
    for (Vertex v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = interior_vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom->position[v];

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];
            double C = mu * 1.0/(4.0*geom->triangle_area(tri));

            // Diagonal term.
            double diagonal_term = C*(vpp_pos - vp_pos).dot(vpp_pos - vp_pos);
            add_entry(v_index, v_index, diagonal_term);

            // vp contribution.
            if (vp.on_boundary()) {
            } else {
                double vp_term = C*(v_pos - vpp_pos).dot(vpp_pos - vp_pos);
                int vp_index = interior_vertex_indices[vp];
                add_entry(v_index, vp_index, vp_term);
            }
            
            // vpp contribution.
            if (vpp.on_boundary()) {
            } else {
                double vpp_term = C*(vp_pos - v_pos).dot(vpp_pos - vp_pos);
                int vpp_index = interior_vertex_indices[vpp];
                add_entry(v_index, vpp_index, vpp_term);
            }

            he = he.twin().next();
        } while (he != start);
    }
    
    
    // Finalize the mass matrix.
    auto mass_matrix = SparseMatrix(num_nodes, num_nodes);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());

    // u/delta_time term.
    mass_matrix += inv_delta_time * gramian_matrix;

    // mass_matrix.makeCompressed();

    // u_prev/delta_time term.
    rhs += inv_delta_time * gramian_matrix * u_vector;

    // Reaction term
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = interior_vertex_indices[v];
        vec3 v_pos = eigen_to_vec3(geom->position[v]);

        double total_triangle_area = 0.0;
        // For each adjacent triangle.
        auto start = v.halfedge(); //todo: provide this circulator in mesh_processing.
        auto he = start;
        do {
            auto tri = he.face();
            total_triangle_area += geom->triangle_area(tri);
            he = he.twin().next();
        } while (he != start);

        double reaction_sample = reaction_function(v_pos, u_mesh[v], time);

        double val = total_triangle_area/3.0 * reaction_sample;
        rhs[v_index] += val;
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
    u_vector = solver.solve(rhs);
    
    /*--------------------------------------------------------------------------------
    // Reassociate each coefficient with the corresponding vertex of the mesh.
    --------------------------------------------------------------------------------*/
    int interior_vertex_index = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            getchar();
        } else {
            u_mesh[v] = u_vector[interior_vertex_index];
            interior_vertex_index += 1;
        } 
    }

    time += delta_time;
}


void ReactorDiffuser::set_u(std::function<double(vec3)> func)
{
    int interior_vertex_index = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        vec3 pos = eigen_to_vec3(geom->position[v]);
        double val = func(pos);
        u_mesh[v] = val;
        u_vector[interior_vertex_index] = val;
        interior_vertex_index += 1;
    }
}
