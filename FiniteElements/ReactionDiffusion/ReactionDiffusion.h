#ifndef HEADER_DEFINED_REACTION_DIFFUSION
#define HEADER_DEFINED_REACTION_DIFFUSION
#include "ReactionDiffusion/core.h"


class ReactorDiffuser {
public:
    ReactorDiffuser(SurfaceGeometry &_geom,
                    float _mu1,
                    float _mu2,
                    std::function<double(vec3, double, double, double)> _u_reaction_function,
                    std::function<double(vec3, double, double, double)> _v_reaction_function
    );

    double mu1;
    double mu2;
    std::function<double(vec3, double, double, double)> u_reaction_function;
    std::function<double(vec3, double, double, double)> v_reaction_function;
    SurfaceGeometry *geom;

    VertexAttachment<int> interior_vertex_indices;

    void time_step(double delta_time);

    double time;
    int num_nodes; // Variation nodes (not on the boundary).
    SparseMatrix gramian_matrix;

    // Chemical concentrations.
    Eigen::VectorXd u_vector;
    VertexAttachment<double> u_mesh;
    Eigen::VectorXd v_vector;
    VertexAttachment<double> v_mesh;

    void set_u(std::function<double(vec3)> func);
    void set_v(std::function<double(vec3)> func);
private:
};


#endif // HEADER_DEFINED_NAVIER_STOKES_SOLVER
