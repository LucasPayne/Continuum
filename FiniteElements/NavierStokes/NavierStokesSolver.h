/*================================================================================
    P2-P1 Taylor-Hood finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/
#ifndef HEADER_DEFINED_NAVIER_STOKES_SOLVER
#define HEADER_DEFINED_NAVIER_STOKES_SOLVER
#include "NavierStokes/P2_P1.h"
#include "NavierStokes/core.h"

// Helper structs for building the Gateaux matrix.
struct TopLeftEntry {
    P2Element velocity_trial_node;
    int trial_component; // x: 0, y: 1
    P2Element velocity_test_node;
    int test_component;  // x: 0, y: 1
    double value;
    TopLeftEntry(P2Element _vtrialn, int _trial_c, P2Element _vtestn, int _tc, double val) :
        velocity_trial_node{_vtrialn}, trial_component{_trial_c}, velocity_test_node{_vtestn}, test_component{_tc}, value{val}
    {}
};
struct BottomLeftEntry {
    Vertex pressure_trial_node;
    P2Element velocity_test_node;
    int test_component;  // x: 0, y: 1
    double value;
    BottomLeftEntry(Vertex _ptn, P2Element _vtn, int _tc, double val) :
        pressure_trial_node{_ptn}, velocity_test_node{_vtn}, test_component{_tc}, value{val}
    {}
};


struct NavierStokesSolver
{
public:
    NavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity);
    
    void start_time_step(double dt);
    void newton_iteration();
    void end_time_step();

    void time_step(double dt);

    void set_source(TimeDependentPlaneVectorField vf);
    void set_velocity(PlaneVectorField vf);

    inline bool solving() const { return m_solving; }
    inline bool iterating() const { return m_iterating; }
    inline int num_velocity_variation_nodes() const { return m_num_velocity_variation_nodes; }
    inline int num_pressure_variation_nodes() const { return m_num_pressure_variation_nodes; }
    inline double kinematic_viscosity() const { return m_kinematic_viscosity; }
    inline double time() const { return m_time; }

    SurfaceGeometry &geom;

    // The velocity and pressure are changed during Newton iteration.
    P2Attachment<vec2> velocity;
    P1Attachment<double> pressure;

    TimeDependentPlaneVectorField source_function; // Exact function.

    bool m_use_advection; // for debugging whether advection actually works
private:
    // The previous velocity and pressure are only changed after a time step.
    P2Attachment<vec2> velocity_prev;
    P1Attachment<double> pressure_prev;

    P2Attachment<int> velocity_node_indices;
    P1Attachment<int> pressure_node_indices;

    inline int system_N() const { return m_system_N; } // The size of the velocity-pressure vectors and Gateaux matrix.
    std::tuple<SparseMatrix, SparseMatrix> compute_gateaux_matrix(); // returns (linear_term_matrix, gateaux_matrix).
    Eigen::VectorXd compute_residual(SparseMatrix &linear_term_matrix);
    void add_nonlinear_velocity_residual(P2Attachment<vec2> &velocity_residual);
    std::vector<TopLeftEntry> compute_linear_term_matrix_top_left();
    std::vector<BottomLeftEntry> compute_linear_term_matrix_bottom_left();

    P2Attachment<vec2> source_samples_P2; // Samples for approximate integration.
    void update_source_samples();

    int m_num_velocity_variation_nodes;
    int m_num_pressure_variation_nodes;
    int m_system_N;
    Eigen::VectorXd m_velocity_pressure_vector;

    double m_kinematic_viscosity;
    bool m_solving; // Has the simulation started?
    bool m_iterating; // Is the algorithm in the middle of a Newton iteration?
    double m_current_time_step_dt;
    double m_time;
};

#endif // HEADER_DEFINED_NAVIER_STOKES_SOLVER
