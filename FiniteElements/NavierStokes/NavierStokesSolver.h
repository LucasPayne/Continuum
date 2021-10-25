/*================================================================================
    P2-P1 Taylor-Hood finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/
#ifndef HEADER_DEFINED_NAVIER_STOKES_SOLVER
#define HEADER_DEFINED_NAVIER_STOKES_SOLVER
#include "P2_P1.h"


struct NavierStokesSolver
{
public:
    NavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity);
    
    void start_time_step(double dt);
    void newton_iteration();
    void end_time_step();

    void set_source(PlaneVectorField vf);

    inline bool solving() const { return m_solving; }
    inline int num_velocity_variation_nodes() const { return m_num_velocity_variation_nodes; }
    inline int num_pressure_variation_nodes() const { return m_num_pressure_variation_nodes; }
    inline double kinematic_viscosity() const { return m_kinematic_viscosity; }
    inline double time() const { return m_time; }

private:
    // The velocity and pressure are changed during Newton iteration.
    P2Attachment<vec2> velocity;
    P1Attachment<double> pressure;
    // The previous velocity and pressure are only changed after a time step.
    P2Attachment<vec2> velocity_prev;
    P1Attachment<double> pressure_prev;

    P1Attachment<int> P1_indices; // 
    P2Attachment<int> P2_indices; // 

    inline int system_N() const { return m_system_N; } // The size of the velocity-pressure vectors and Gateaux matrix.
    SparseMatrix compute_gateaux_matrix();
    Eigen::VectorXd compute_residual();

    PlaneFunction<vec2> source_function; // The source form is an exact function, not sampled.

    int m_num_velocity_variation_nodes;
    int m_num_pressure_variation_nodes;
    Eigen::VectorXd m_velocity_pressure_vector;

    double m_kinematic_viscosity;
    bool m_solving; // Has the simulation started?
    bool m_iterating; // Is the algorithm in the middle of a Newton iteration?
    double m_current_time_step_dt;
    double m_time;

    // Helper functions for building up the residual vector.
    void add_velocity_residual_advection(P2Attachment<vec2> &velocity_residual);
    void add_velocity_residual_viscosity(P2Attachment<vec2> &velocity_residual);
    void add_velocity_residual_pressure(P2Attachment<vec2> &velocity_residual);
    void add_velocity_residual_source(P2Attachment<vec2> &velocity_residual);
    void add_velocity_residual_time_step(P2Attachment<vec2> &velocity_residual);
};

#endif // HEADER_DEFINED_NAVIER_STOKES_SOLVER
