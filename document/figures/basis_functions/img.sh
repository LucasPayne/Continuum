#!/bin/bash

convert +append quadratic_basis_100.ppm quadratic_basis_010.ppm quadratic_basis_001.ppm tmp_1.png
convert +append quadratic_basis_110.ppm quadratic_basis_011.ppm quadratic_basis_101.ppm tmp_2.png
convert -append tmp_1.png tmp_2.png quadratic_basis.png
