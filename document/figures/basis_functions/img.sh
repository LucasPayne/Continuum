#!/bin/bash

convert +append quadratic_basis_1_100.ppm quadratic_basis_2_010.ppm quadratic_basis_3_001.ppm tmp_1.png
convert +append quadratic_basis_4_110.ppm quadratic_basis_5_011.ppm quadratic_basis_6_101.ppm tmp_2.png
convert -append tmp_1.png tmp_2.png quadratic_basis.png
