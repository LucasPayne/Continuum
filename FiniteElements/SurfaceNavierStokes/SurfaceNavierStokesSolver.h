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
    P2Element centripetal_trial_node;
    P2Element velocity_test_node;
    vec3 value;
    CentripetalBlockEntry(P2Element _ptn, P2Element _vtn, vec3 val) :
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


    // Data for debugging.
    //--------------------------------------------------------------------------------
    P2Attachment<vec3> _test_point_1;
    P2Attachment<vec3> _test_point_2;
    //--------------------------------------------------------------------------------

    P2Attachment<vec3> velocity;
    P1Attachment<double> pressure;
    P2Attachment<double> centripetal;

    FaceAttachment<vec3> triangle_normal;
    FaceAttachment<mat3x3> triangle_projection_matrix;
    P2Attachment<vec3> normal;

    P2Attachment<vec3> source_samples_P2; // Samples for approximate integration.
    void set_source(std::function<vec3(double,double,double)> vf);

    // For debugging.
    void set_velocity(std::function<vec3(double,double,double)> vf);

    void make_sparsity_image(SparseMatrix &matrix, std::string name);

    void explicit_advection();

    double m_current_time_step_dt;

    bool m_advect;

    std::tuple<Face, vec3> traverse(Face tri, vec3 origin, vec3 shift, int depth=0, int ignore_index=-1);
private:
    P2Attachment<int> velocity_node_indices;
    P1Attachment<int> pressure_node_indices;
    P2Attachment<int> centripetal_node_indices;

    inline int system_N() const { return m_system_N; } // The size of the velocity-pressure-centripetal vectors.

    std::vector<PressureBlockEntry> compute_pressure_block_coefficients();
    std::vector<VelocityBlockEntry> compute_velocity_block_coefficients();
    std::vector<CentripetalBlockEntry> compute_centripetal_block_coefficients();
    SparseMatrix compute_matrix();
    Eigen::VectorXd compute_rhs();


    int m_num_velocity_variation_nodes;
    int m_num_pressure_variation_nodes;
    int m_num_centripetal_variation_nodes;
    int m_system_N;
    Eigen::VectorXd m_solution_vector;

    double m_kinematic_viscosity;
    bool m_solving; // Has the simulation started?
    double m_time;
};

#endif // HEADER_DEFINED_SURFACE_NAVIER_STOKES_SOLVER
