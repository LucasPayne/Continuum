#ifndef HEADER_DEFINED_REACTION_DIFFUSION
#define HEADER_DEFINED_REACTION_DIFFUSION
#include "ReactionDiffusion/core.h"


class ReactorDiffuser {
public:
    ReactorDiffuser(SurfaceGeometry &_geom, float _mu, std::function<double(vec3, double, double)> _reaction_function);

    double mu;
    std::function<double(vec3, double, double)> reaction_function;
    SurfaceGeometry *geom;

    VertexAttachment<int> interior_vertex_indices;

    void time_step(double delta_time);

    double time;
    int num_nodes; // Variation nodes (not on the boundary).
    SparseMatrix gramian_matrix;

    Eigen::VectorXd u_vector;
    VertexAttachment<double> u_mesh;

    void set_u(std::function<double(vec3)> func);
private:
};


#endif // HEADER_DEFINED_NAVIER_STOKES_SOLVER
