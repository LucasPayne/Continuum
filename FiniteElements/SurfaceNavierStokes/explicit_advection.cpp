#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"
#include "core.h"

void SurfaceNavierStokesSolver::explicit_advection()
{
    std::function<double(double,double,double)> u_basis[6] = {
        [](double x, double y, double z)->double {
            return x - 2*x*y - 2*x*z;
        },
        [](double x, double y, double z)->double {
            return y - 2*y*z - 2*y*x;
        },
        [](double x, double y, double z)->double {
            return z - 2*z*x - 2*z*y;
        },
        [](double x, double y, double z)->double {
            return 4*x*y;
        },
        [](double x, double y, double z)->double {
            return 4*y*z;
        },
        [](double x, double y, double z)->double {
            return 4*z*x;
        },
    };

    auto new_velocity = P2Attachment<vec3>(geom.mesh);

    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            new_velocity[v] = vec3(0,0,0);
            continue;
        }
        vec3 u = velocity[v];
        vec3 n = normal[v];

        auto start = v.halfedge();
        auto he = start;
        vec3 p1 = eigen_to_vec3(geom.position[v]);
        do {
            auto tri = he.face();
            Vertex v2 = he.next().vertex();
            Vertex v3 = he.next().next().vertex();
            vec3 p2 = eigen_to_vec3(geom.position[v2]);
            vec3 p3 = eigen_to_vec3(geom.position[v3]);

            // Gram-Schmidt to get a basis on this triangle.
            vec3 tri_basis_1 = (p2 - p1).normalized();
            vec3 tri_basis_2 = p3 - p1;
            tri_basis_2 -= tri_basis_1 * vec3::dot(tri_basis_2, tri_basis_1);
            tri_basis_2 = tri_basis_2.normalized();
            // Rotate a vector 90 degrees anticlockwise on the triangle plane.
            auto perp = [&](vec3 vec)->vec3 {
                float a = vec3::dot(vec, tri_basis_1);
                float b = vec3::dot(vec, tri_basis_2);
                return -tri_basis_2*a + tri_basis_1*b;
            };

            vec3 shift = triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);
            // std::cout << triangle_projection_matrix[tri] << "\n";
            // std::cout << u << "\n";
            // printf("%.2g\n", m_current_time_step_dt);
            // std::cout << -m_current_time_step_dt*u << "\n";
            // std::cout << shift << "\n";getchar();
    
	    _test_point_1[v] = shift;
	    _test_point_2[v] = eigen_to_vec3(geom.barycenter(tri));
            if (vec3::dot(shift, perp(p2 - p1)) <= 0 &&
                    vec3::dot(shift, perp(p3 - p1)) >= 0) {
                // Compute barycentric coordinates.
                vec3 c = shift + p1;
                double x = vec3::cross(p2 - c, p3 - c).length();
                double y = vec3::cross(p3 - c, p1 - c).length();
                double z = vec3::cross(p1 - c, p2 - c).length();
                double w = x+y+z;
                x /= w;
                y /= w;
                z /= w;
                // printf("%.2g %.2g %.2g\n", x, y, z);getchar();
	        _test_point_1[v] = p1*x + p2*y + p3*z;
            
                auto edge_110 = he.next().edge();
                auto edge_011 = he.next().next().edge();
                auto edge_101 = he.edge();
                P2Element elements[6] = {v, v2, v3, edge_110, edge_011, edge_101};
                vec3 val = vec3(0,0,0);
                for (int i = 0; i < 6; i++) {
                    val += u_basis[i](x,y,z) * velocity[elements[i]];
                }
                new_velocity[v] = val;
                break;
            }
            he = he.twin().next();
        } while (he != start);
    }

    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            new_velocity[edge] = vec3(0,0,0);
            continue;
        }
        vec3 u = velocity[edge];
        vec3 n = normal[edge];

        vec3 edge_vec = eigen_to_vec3(geom.position[edge.b().vertex()] - geom.position[edge.a().vertex()]);
        vec3 midpoint = eigen_to_vec3(geom.midpoint(edge));
        Face tri;
        if (vec3::dot(n, vec3::cross(-u, edge_vec)) > 0) {
            tri = edge.b().face();
        } else {
            tri = edge.a().face();
        }
        _test_point_2[edge] = eigen_to_vec3(geom.barycenter(tri));

        auto v1 = tri.halfedge().vertex();
        auto v2 = tri.halfedge().next().vertex();
        auto v3 = tri.halfedge().next().next().vertex();
        vec3 p1 = eigen_to_vec3(geom.position[v1]);
        vec3 p2 = eigen_to_vec3(geom.position[v2]);
        vec3 p3 = eigen_to_vec3(geom.position[v3]);
        auto edge_110 = tri.halfedge().edge();
        auto edge_011 = tri.halfedge().next().edge();
        auto edge_101 = tri.halfedge().next().next().edge();

        vec3 shift = triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);
        vec3 c = shift + midpoint;
	double x = vec3::cross(p2 - c, p3 - c).length();
	double y = vec3::cross(p3 - c, p1 - c).length();
	double z = vec3::cross(p1 - c, p2 - c).length();
	double w = x+y+z;
        x /= w;
        y /= w;
        z /= w;
	_test_point_1[edge] = p1*x + p2*y + p3*z;

	P2Element elements[6] = {v1, v2, v3, edge_110, edge_011, edge_101};
        vec3 val = vec3(0,0,0);
        for (int i = 0; i < 6; i++) {
            val += u_basis[i](x,y,z) * velocity[elements[i]];
        }
        new_velocity[edge] = val;
    }

    for (auto v : geom.mesh.vertices()) {
        velocity[v] = new_velocity[v];
    }
    for (auto edge : geom.mesh.edges()) {
        velocity[edge] = new_velocity[edge];
    }
}
