#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"



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
        if (i < 0 || j < 0 || i >= m_system_N || j >= m_system_N) {
            fprintf(stderr, "%d,%d,  %d x %d, index out of bounds.\n", i,j, m_system_N,m_system_N);
            exit(EXIT_FAILURE);
        }
        eigen_coefficients.push_back(EigenTriplet(i, j, val));
    };
    for (auto coeff : velocity_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                add_eigen_coefficient(
                    3*velocity_node_indices[coeff.velocity_trial_node] + i,
                    3*velocity_node_indices[coeff.velocity_test_node] + j,
                    coeff.value.entry(i,j)
                );
            }
        }
    }
    for (auto coeff : pressure_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            add_eigen_coefficient(
                3*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
                3*velocity_node_indices[coeff.velocity_test_node] + i,
                coeff.value[i]
            );
            add_eigen_coefficient(
                3*velocity_node_indices[coeff.velocity_test_node] + i,
                3*m_num_velocity_variation_nodes + pressure_node_indices[coeff.pressure_trial_node],
                coeff.value[i]
            );
        }
    }
    for (auto coeff : centripetal_block_coefficients) {
        for (int i = 0; i < 3; i++) {
            add_eigen_coefficient(
                3*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1 + centripetal_node_indices[coeff.centripetal_trial_node],
                3*velocity_node_indices[coeff.velocity_test_node] + i,
                coeff.value[i]
            );
            add_eigen_coefficient(
                3*velocity_node_indices[coeff.velocity_test_node] + i,
                3*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1 + centripetal_node_indices[coeff.centripetal_trial_node],
                coeff.value[i]
            );
        }
    }


    auto matrix = SparseMatrix(m_system_N, m_system_N);
    printf("%d x %d\n", m_system_N, m_system_N);
    matrix.setFromTriplets(eigen_coefficients.begin(), eigen_coefficients.end());
    matrix.makeCompressed();
    return matrix;
}
