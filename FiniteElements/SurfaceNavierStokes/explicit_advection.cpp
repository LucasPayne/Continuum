#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"
#include "core.h"

std::tuple<Face, vec3> SurfaceNavierStokesSolver::traverse(Face tri, vec3 origin, vec3 shift, int depth, int ignore_index)
{
    const int MAX_DEPTH = 50;
    if (depth == MAX_DEPTH) return {tri, origin};

    Halfedge hes[3] = {
        tri.halfedge(),
        tri.halfedge().next(),
        tri.halfedge().next().next()
    };
    Vertex vs[3] = {
        tri.halfedge().vertex(),
        tri.halfedge().next().vertex(),
        tri.halfedge().next().next().vertex()
    };
    vec3 ps[3];
    for (int i = 0; i < 3; i++) {
        ps[i] = eigen_to_vec3(geom.position[vs[i]]);
    }
    // Gram-Schmidt for a basis on this triangle.
    vec3 e1 = (ps[1] - ps[0]).normalized();
    vec3 e2 = ps[2] - ps[0];
    e2 -= e1*vec3::dot(e2, e1);
    e2 = e2.normalized();

    // Express the triangle in this basis.
    vec2 ps_t[3];
    for (int i = 0; i < 3; i++) {
        ps_t[i] = vec2(vec3::dot(ps[i]-ps[0], e1), vec3::dot(ps[i]-ps[0], e2));
    }
    // Create a ray.
    vec2 o_t = vec2(vec3::dot(origin - ps[0], e1), vec3::dot(origin-ps[0], e2));
    vec2 d_t = vec2(vec3::dot(shift, e1), vec3::dot(shift, e2));

    int line_hit_index = -1;
    double min_t = std::numeric_limits<double>::max();
    // Intersect each line determined by the triangle sides.
    for (int i = 0; i < 3; i++) {
        if (i == ignore_index) continue;
        vec2 line_n = (ps_t[(i+1)%3] - ps_t[i]).perp();
        vec2 line_p = ps_t[i];
        double t = vec2::dot(line_p - o_t, line_n)/vec2::dot(d_t, line_n);
        if (t >= 0 && t < min_t) {
            min_t = t;
            line_hit_index = i;
        }
    }

    if (line_hit_index == -1 || min_t > 1) {
        // Travel stops on this triangle.
        // printf("Travel stops.\n");
        return {tri, origin+shift};
    }
    Halfedge hit_he = hes[line_hit_index];
    if (hit_he.twin().face().null()) {
        // Hit the boundary. Stop at the boundary intersection.
        // printf("Hit the boundary.\n");
        return {tri, origin + min_t*shift};
    }
    // Travel proceeds on another triangle.
    // Create an orthonormal basis for each incident face to the edge being travelled over.
    // This basis shares the E1 vector.
    vec3 E1 = eigen_to_vec3(geom.vector(hit_he)).normalized();
    vec3 from_E2 = -eigen_to_vec3(geom.vector(hit_he.next()));
    from_E2 -= E1*vec3::dot(from_E2, E1);
    from_E2 = from_E2.normalized();
    vec3 to_E2 = -eigen_to_vec3(geom.vector(hit_he.twin().next().next()));
    to_E2 -= E1*vec3::dot(to_E2, E1);
    to_E2 = to_E2.normalized();
    
    vec3 new_shift = E1*vec3::dot((1-min_t)*shift, E1) + to_E2*vec3::dot((1-min_t)*shift, from_E2);

    Face to_face = hit_he.twin().face();
    // Ignore the edge that was traversed over, when intersecting on the next triangle.
    int to_ignore_index = 0;
    {
        auto start = to_face.halfedge();
        auto he = start;
        do {
            if (he == hit_he.twin()) break;
            to_ignore_index += 1;
            he = he.next();
        } while (he != start);
        assert(to_ignore_index != 3);
    }
    float fix = -0.001;
    return traverse(to_face, origin + min_t*shift + fix*new_shift, new_shift, depth+1, to_ignore_index);
}



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

    auto barycentric_coeff = [&](vec3 c, vec3 pa, vec3 pb, vec3 nor)->double {
        vec3 cross = vec3::cross(pb - c, pa - c);
        if (vec3::dot(cross, nor) < 0) return cross.length();
        return -cross.length();
    };

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
            Vertex v1 = v;
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

            vec3 shift = 0.5 * triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);
            // std::cout << triangle_projection_matrix[tri] << "\n";
            // std::cout << u << "\n";
            // printf("%.2g\n", m_current_time_step_dt);
            // std::cout << -m_current_time_step_dt*u << "\n";
            // std::cout << shift << "\n";getchar();
    
	    _test_point_1[v] = shift;
            if (vec3::dot(shift, perp(p2 - p1)) <= 0 &&
                    vec3::dot(shift, perp(p3 - p1)) >= 0) {
                double fix = 0.001;
                Face out_face;
                vec3 out_pos;
                std::tie(out_face, out_pos) = traverse(tri, p1+(1+fix)*shift, (1-fix)*shift);
                assert(!out_face.null());

                Edge edges[3] = {
                    out_face.halfedge().edge(),
                    out_face.halfedge().next().edge(),
                    out_face.halfedge().next().next().edge()
                };
                Vertex vs[3] = {
                    out_face.halfedge().vertex(),
                    out_face.halfedge().next().vertex(),
                    out_face.halfedge().next().next().vertex()
                };
                vec3 ps[3];
                for (int i = 0; i < 3; i++) ps[i] = eigen_to_vec3(geom.position[vs[i]]);
                
                // Compute barycentric coordinates.
                double x = barycentric_coeff(out_pos, ps[1], ps[2], n);
                double y = barycentric_coeff(out_pos, ps[2], ps[0], n);
                double z = barycentric_coeff(out_pos, ps[0], ps[1], n);
                double w = x+y+z;
                x /= w;
                y /= w;
                z /= w;
	        _test_point_1[v] = ps[0]*x + ps[1]*y + ps[2]*z;
	        _test_point_2[v] = eigen_to_vec3(geom.barycenter(out_face));
                vec3 val = vec3(0,0,0);

                P2Element elements[6] = {vs[0], vs[1], vs[2], edges[0], edges[1], edges[2]};
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

        auto v1 = tri.halfedge().vertex();
        auto v2 = tri.halfedge().next().vertex();
        auto v3 = tri.halfedge().next().next().vertex();
        vec3 p1 = eigen_to_vec3(geom.position[v1]);
        vec3 p2 = eigen_to_vec3(geom.position[v2]);
        vec3 p3 = eigen_to_vec3(geom.position[v3]);
        auto edge_110 = tri.halfedge().edge();
        auto edge_011 = tri.halfedge().next().edge();
        auto edge_101 = tri.halfedge().next().next().edge();
        vec3 tri_n = triangle_normal[tri];

        vec3 shift = 0.5 * triangle_projection_matrix[tri] * (-m_current_time_step_dt*u);

	double fix = 0.001;
	Face out_face;
	vec3 out_pos;
	std::tie(out_face, out_pos) = traverse(tri, midpoint+(1+fix)*shift, (1-fix)*shift);
	assert(!out_face.null());

        Edge edges[3] = {
            out_face.halfedge().edge(),
            out_face.halfedge().next().edge(),
            out_face.halfedge().next().next().edge()
        };
        Vertex vs[3] = {
            out_face.halfedge().vertex(),
            out_face.halfedge().next().vertex(),
            out_face.halfedge().next().next().vertex()
        };
        vec3 ps[3];
        for (int i = 0; i < 3; i++) ps[i] = eigen_to_vec3(geom.position[vs[i]]);
        
        // Compute barycentric coordinates.
        double x = barycentric_coeff(out_pos, ps[1], ps[2], n);
        double y = barycentric_coeff(out_pos, ps[2], ps[0], n);
        double z = barycentric_coeff(out_pos, ps[0], ps[1], n);
        double w = x+y+z;
        x /= w;
        y /= w;
        z /= w;
	_test_point_1[edge] = ps[0]*x + ps[1]*y + ps[2]*z;
        _test_point_2[edge] = eigen_to_vec3(geom.barycenter(out_face));
        vec3 val = vec3(0,0,0);

	P2Element elements[6] = {vs[0], vs[1], vs[2], edges[0], edges[1], edges[2]};
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
