#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"

void SurfaceNavierStokesSolver::make_sparsity_image(SparseMatrix &matrix, std::string name)
{
    assert(matrix.rows() == matrix.cols());
    int system_N = matrix.rows();
    int num_nonzeros = 0;
    for (int i = 0; i < system_N; i++) {
        for (int j = 0; j < system_N; j++) {
            if (fabs(matrix.coeff(i, j)) >= 1e-4) {
                num_nonzeros += 1;
            }
        }
    }
    FILE *ppm_file = fopen(const_cast<const char *>(name.c_str()), "w+");
    fprintf(ppm_file, "P3\n");
    fprintf(ppm_file, "%d %d\n", system_N, system_N);
    fprintf(ppm_file, "255\n");
    for (int i = 0; i < system_N; i++) {
        for (int j = 0; j < system_N; j++) {
            if (fabs(matrix.coeff(i, j)) >= 1e-4) {
                if (fabs(matrix.coeff(i, j) - matrix.coeff(j, i)) <= 1e-4) {
                    // signify when this entry is symmetric (equal to its corresponding transpose entry).
                    fprintf(ppm_file, "255 0 0 ");
                } else {
                    fprintf(ppm_file, "0 0 0 ");
                }
            } else {
                if (i < 3*m_num_velocity_variation_nodes && j < 3*m_num_velocity_variation_nodes) {
                    if (((i/3)%2 == 0) != ((j/3)%2 == 0)) {
	                fprintf(ppm_file, "255 255 255 ");
                    } else {
	                fprintf(ppm_file, "0 255 255 ");
                    }
                } else if (i < 3*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1
                            && j < 3*m_num_velocity_variation_nodes + m_num_pressure_variation_nodes-1) {
                    if (((i/3)%2 == 0) != ((j/3)%2 == 0)) {
	                fprintf(ppm_file, "230 230 230 ");
                    } else {
	                fprintf(ppm_file, "0 230 230 ");
                    }
                } else {
                    if (((i/3)%2 == 0) != ((j/3)%2 == 0)) {
	                fprintf(ppm_file, "255 255 255 ");
                    } else {
	                fprintf(ppm_file, "0 255 255 ");
                    }
                }
            }
        }
    }
    fclose(ppm_file);
}
