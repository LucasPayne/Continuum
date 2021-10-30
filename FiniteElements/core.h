#ifndef HEADER_DEFINED_CORE
#define HEADER_DEFINED_CORE

#include "cg_sandbox.h" //maybe shouldn't be core to the solver
#include "triangle_wrapper.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

using PlaneFunction = std::function<double(double x, double y)>; // Function of the XY plane.
using PlaneVectorField = std::function<vec2(double x, double y)>;
using TimeDependentPlaneVectorField = std::function<vec2(double x, double y, double t)>;
using PlaneFunctionNL1 = std::function<double(double x, double y, double u)>; // First-order non-linear plane function.
using SparseMatrix = Eigen::SparseMatrix<double>;
using EigenTriplet = Eigen::Triplet<double>;


vec3 eigen_to_vec3(Eigen::Vector3f v);
Eigen::Vector3f vec3_to_eigen(vec3 v);


#endif // HEADER_DEFINED_CORE
