/*================================================================================
    P2-P1 Taylor-Hood finite element solver for the incompressible Navier-Stokes
    equations.
    
    The boundary condition is no-slip for the entire boundary.
================================================================================*/
#ifndef HEADER_DEFINED_SURFACE_NAVIER_STOKES_SOLVER
#define HEADER_DEFINED_SURFACE_NAVIER_STOKES_SOLVER
#include "core.h"
#include "P2_P1.h"

// Helper structs for building the matrix.
struct VelocityBlockEntry {
    P2Element velocity_trial_node;
    P2Element velocity_test_node;
    mat3x3 value;
    VelocityBlockEntry(P2Element _vtrialn, P2Element _vtestn, mat3x3 val) :
        velocity_trial_node{_vtrialn}, velocity_test_node{_vtestn}, value{val}
    {}
};
struct PressureBlockEntry {
    Vertex pressure_trial_node;
    P2Element velocity_test_node;
    vec3 value;
    PressureBlockEntry(Vertex _ptn, P2Element _vtn, vec3 val) :
        pressure_trial_node{_ptn}, velocity_test_node{_vtn}, value{val}
    {}
};
struct CentripetalBlockEntry {
    Vertex centripetal_trial_node;
    P2Element velocity_test_node;
    vec3 value;
    CentripetalBlockEntry(Vertex _ptn, P2Element _vtn, vec3 val) :
        centripetal_trial_node{_ptn}, velocity_test_node{_vtn}, value{val}
    {}
};


struct SurfaceNavierStokesSolver
{
public:
    SurfaceNavierStokesSolver(SurfaceGeometry &_geom, double _kinematic_viscosity);
    
    void time_step(double dt);

    inline bool solving() const { return m_solving; }
    inline int num_velocity_variation_nodes() const { return m_num_velocity_variation_nodes; }
    inline int num_pressure_variation_nodes() const { return m_num_pressure_variation_nodes; }
    inline int num_centripetal_variation_nodes() const { return m_num_centripetal_variation_nodes; }
    inline double kinematic_viscosity() const { return m_kinematic_viscosity; }
    inline double time() const { return m_time; }

    SurfaceGeometry &geom;

    P2Attachment<vec3> velocity;
    P1Attachment<double> pressure;
    P1Attachment<double> centripetal;

    FaceAttachment<vec3> triangle_normal;
    FaceAttachment<mat3x3> triangle_projection_matrix;

    P2Attachment<vec3> source_samples_P2; // Samples for approximate integration.

    void make_sparsity_image(SparseMatrix &matrix, std::string name);
private:
    // The previous velocity and pressure are only changed after a time step.
    P2Attachment<vec3> velocity_prev;
    P1Attachment<double> pressure_prev;
    P1Attachment<double> centripetal_prev;

    P2Attachment<int> velocity_node_indices;
    P1Attachment<int> pressure_node_indices;

    inline int system_N() const { return m_system_N; } // The size of the velocity-pressure-centripetal vectors.
    void explicit_advection_lagrangian();

    std::vector<PressureBlockEntry> compute_pressure_block_coefficients();
    std::vector<VelocityBlockEntry> compute_velocity_block_coefficients();
    std::vector<CentripetalBlockEntry> compute_centripetal_block_coefficients();
    SparseMatrix compute_matrix();

    int m_num_velocity_variation_nodes;
    int m_num_pressure_variation_nodes;
    int m_num_centripetal_variation_nodes;
    int m_system_N;
    Eigen::VectorXd m_solution_vector;

    double m_kinematic_viscosity;
    bool m_solving; // Has the simulation started?
    double m_current_time_step_dt;
    double m_time;
};

#endif // HEADER_DEFINED_SURFACE_NAVIER_STOKES_SOLVER
