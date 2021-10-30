#include "core.h"

// utilities
vec3 eigen_to_vec3(Eigen::Vector3f v)
{
    return vec3(v.x(), v.y(), v.z());
}
Eigen::Vector3f vec3_to_eigen(vec3 v)
{
    return Eigen::Vector3f(v.x(), v.y(), v.z());
}

