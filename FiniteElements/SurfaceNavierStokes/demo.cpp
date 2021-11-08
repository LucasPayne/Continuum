#include "core.h"
#include "SurfaceNavierStokes/demo.h"
#include "mesh_generators.cpp"
#include "mesh_processing/extensions/assimp_convert.h"

bool filming = false;
vec3 shift = vec3(0, 0.001, 0);

Face test_face;
double test_shift_t = 0.;

Demo::Demo()
{
    solver = nullptr;
    geom = nullptr;
}

std::tuple<Face, vec3> Demo::traverse(Face tri, vec3 origin, vec3 shift)
{
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
        ps[i] = eigen_to_vec3(geom->position[vs[i]]);
    }
    // Gram-Schmidt for a basis on this triangle.
    vec3 e1 = (ps[1] - ps[0]).normalized();
    vec3 e2 = ps[2] - ps[0];
    e2 -= e1*vec3::dot(e2, e1);
    e2 = e2.normalized();

    world->graphics.paint.line(origin, origin+e1, 0.01, vec4(1,1,0,1));
    world->graphics.paint.line(origin, origin+e2, 0.01, vec4(0,1,0,1));

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
        return {tri, origin+shift};
    }
    Halfedge hit_he = hes[line_hit_index];
    if (hit_he.twin().face().null()) {
        // Hit the boundary. Stop at the boundary intersection.
        return {tri, origin + min_t*shift};
    }
    // Travel proceeds on another triangle.
    // Create an orthonormal basis for each incident face to the edge being travelled over.
    // This basis shares the E1 vector.
    vec3 E1 = eigen_to_vec3(geom->vector(hit_he)).normalized();
    vec3 from_E2 = -eigen_to_vec3(geom->vector(hit_he.next()));
    from_E2 -= E1*vec3::dot(from_E2, E1);
    from_E2 = from_E2.normalized();
    vec3 to_E2 = -eigen_to_vec3(geom->vector(hit_he.twin().next().next()));
    to_E2 -= E1*vec3::dot(to_E2, E1);
    to_E2 = to_E2.normalized();
    
    vec3 new_shift = E1*vec3::dot((1-min_t)*shift, E1) + to_E2*vec3::dot((1-min_t)*shift, from_E2);

    double fix = 0.001; // (To prevent intersections with the edge just passed through.)
    return traverse(hit_he.twin().face(), origin + min_t*shift + new_shift*fix, (1-fix)*new_shift);
}

void Demo::init()
{
    // Create a camera controller.
    auto cameraman = world->entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    camera->background_color = vec4(1,1,1,1);
    auto t = cameraman.add<Transform>(0,2,0);
    // main_camera = camera;
    controller = world->add<CameraController>(cameraman);
    controller->angle = -M_PI/2;
    controller->azimuth = M_PI;

    // geom = assimp_to_surface_geometry(std::string(MODELS) + "bunny_head.stl");
    // // geom = assimp_to_surface_geometry(std::string(MODELS) + "tangram.stl");
    // // geom = assimp_to_surface_geometry(std::string(MODELS) + "cylinder.stl");
    // geom->mesh.lock();
    
    // geom = square_mesh(15);
    // for (auto v : geom->mesh.vertices()) {
    //     vec3 p = eigen_to_vec3(geom->position[v]);
    //     float theta = 0.95*(p.x()+1)*M_PI;
    //     // geom->position[v] += Eigen::Vector3f(0, 0.1*sin(4*p.x()), 0);
    //     // geom->position[v] = Eigen::Vector3f(cos(theta), sin(theta), p.z());
    // }
    //     double theta0 = 0.13;
    //     vec2 obstruction_position = vec2(0.2,0.2);
    // geom = square_minus_circle(0.25, theta0, 1, 1, 60, false, obstruction_position, false);

        // double theta0 = 0.1257;
        // vec2 obstruction_position = vec2(0,0);
        // geom = square_minus_circle(0.18, theta0, 1, 1, 60, true, obstruction_position, false);

        geom = square_mesh(25);

    // for (auto v : geom->mesh.vertices()) {
    //     vec3 p = eigen_to_vec3(geom->position[v]);
    //     // float theta = 0.95*(p.x()+1)*M_PI;
    //     // geom->position[v] = Eigen::Vector3f(cos(theta), sin(theta), p.z());
    //     geom->position[v] += Eigen::Vector3f(0, 0.1*sin(4*p.x()), 0);
    // }
    
    double viscosity = 0.001;
    solver = new SurfaceNavierStokesSolver(*geom, viscosity);

    solver->set_source([&](double x, double y, double z)->vec3 {
        const double r = 0.125;
        if ((vec2(x,z) - vec2(-0.85,0)).length() <= r) {
            return vec3(3000, 0, 0);
        }
        return vec3(0,0,0);
        // #if 0
        // double r = 0.125;
        // if ((x+0.8)*(x+0.8) + z*z <= r*r) return vec3(30,0,0);
        // // if (y*y + z*z <= r*r) return vec3(1,0,1);
        // return vec3(0,0,0);
        // #else
        // double r = 0.3;
        // if ((x+0.2)*(x+0.2) + z*z <= r*r) return vec3(300,0,300);
        // return vec3(0,0,0);
        // #endif
    });

    int count = 0;
    for (auto f : geom->mesh.faces()) {
        if (count == 27) {
            test_face = f;
            break;
        }
        count += 1;
    }
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_R) {
            solver->time_step(1./300.);
        }
        if (e.key.code == KEY_P) {
            solver->m_current_time_step_dt = 1./300.;
            solver->explicit_advection();
        }
        if (e.key.code == KEY_B) {
            solver->m_advect = !solver->m_advect;
        }
        if (e.key.code == KEY_1) {
            solver->set_velocity([&](double x, double y, double z)->vec3 {
                double r = 0.50;
                if ((x)*(x) + z*z <= r*r) return 6.7745353*vec3(1,0,1);
                return vec3(0,0,0);
            });
        }
        if (e.key.code == KEY_9) {
            filming = true;
        }
    }
}

void Demo::update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
    double velocity_mul = 0.02;


    auto paint_velocity = [&](vec3 p, vec3 u) {
        const float rho = 1;
        float len = 1 - exp(-rho*sqrt(u.x()*u.x() + u.z()*u.z()));
	len *= 0.035;
        world->graphics.paint.line(p + shift, p + shift + len*u/u.length(), 0.001, vec4(0,0,1,1));
    };

    for (auto v : geom->mesh.vertices()) {
        vec3 p = eigen_to_vec3(geom->position[v]);
        vec3 u = solver->velocity[v];
        // world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        paint_velocity(p, u);
    }
    for (auto e : geom->mesh.edges()) {
        vec3 a = eigen_to_vec3(geom->position[e.a().vertex()]);
        vec3 b = eigen_to_vec3(geom->position[e.b().vertex()]);
        // world->graphics.paint.line(a,b,0.001,vec4(0,0,0,1));
        vec3 p = 0.5*a + 0.5*b;
        vec3 u = solver->velocity[e];
        // world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        paint_velocity(p, u);
    }
    #if 0
    for (auto tri : geom->mesh.faces()) {
        vec3 c = eigen_to_vec3(geom->barycenter(tri));
        vec3 n = solver->triangle_normal[tri];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
        // vec3 k = solver->triangle_projection_matrix[tri] * vec3(1,0,0);
        // world->graphics.paint.line(c,c+0.1*k,0.001,vec4(0.5,0,0.5,1));
    }
    for (auto e : geom->mesh.edges()) {
        vec3 c = eigen_to_vec3(geom->midpoint(e));
        vec3 n = solver->normal[e];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
    }
    for (auto v : geom->mesh.vertices()) {
        vec3 c = eigen_to_vec3(geom->position[v]);
        vec3 n = solver->normal[v];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
    }
    #endif

    #if 1
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        vec3 c = eigen_to_vec3(geom->position[v]);
        vec3 r =  solver->_test_point_1[v];
        //world->graphics.paint.line(c,c+100*(r-c),0.005,vec4(0,0,0,1));
        world->graphics.paint.line(c+shift,r+shift,0.005,vec4(1,0,0,1));
        world->graphics.paint.sphere(r, 0.006, vec4(0,0,0,1));
        world->graphics.paint.line(c+shift,solver->_test_point_2[v]+shift,0.001,vec4(0,1,0,1));
    }
    for (auto e : geom->mesh.edges()) {
        if (e.on_boundary()) continue;
        vec3 c = eigen_to_vec3(geom->midpoint(e));
        vec3 r =  solver->_test_point_1[e];
        world->graphics.paint.line(c+shift,r+shift,0.005,vec4(1,0,0,1));
        world->graphics.paint.sphere(r, 0.006, vec4(0,0,0,1));
        world->graphics.paint.line(c+shift,solver->_test_point_2[e]+shift,0.005,vec4(0.5,1,0,1));
    }
    #endif

    static int counter = 0;
    if (filming) {
        solver->time_step(1./300.);
        // solver->time_step(1./10.);
        save_solution(std::string(DATA) + "flows/sns_2/navier_stokes_" + std::to_string(counter) + ".txt");
        counter += 1;
    }

    #if 0
    const float test_move_speed = 1;
    if (world->input.keyboard.down(KEY_LEFT_ARROW)) {
        test_shift_t += test_move_speed * dt;
    }
    if (world->input.keyboard.down(KEY_RIGHT_ARROW)) {
        test_shift_t -= test_move_speed * dt;
    }
    vec3 c = eigen_to_vec3(geom->barycenter(test_face));
    vec3 s = solver->triangle_projection_matrix[test_face]*vec3(test_shift_t, 0, test_shift_t);
    world->graphics.paint.sphere(c, 0.01, vec4(0.5,0.5,0.5,1));
    world->graphics.paint.line(c + shift, c+s + shift, 0.005,  vec4(0.5,1,0.5,1));

    Face out_face;
    vec3 out_pos;
    std::tie(out_face, out_pos) = traverse(test_face, c, s);
    world->graphics.paint.sphere(out_pos, 0.01, vec4(1,0.2,0.2,1));
    world->graphics.paint.sphere(eigen_to_vec3(geom->barycenter(out_face)), 0.02, vec4(0.5,1,0.2,1));
    #endif
}


void Demo::post_render_update()
{
}


void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
    }
}

void Demo::save_solution(std::string filename)
{
    FILE *file = fopen(filename.c_str(), "w+");
    for (auto v : geom->mesh.vertices()) {
        // fprintf(file, "%.7f %.7f %.7f\n", solver->velocity[v].x(), solver->velocity[v].y(), solver->velocity[v].z());
        fprintf(file, "%.7f %.7f\n", solver->velocity[v].x(), solver->velocity[v].z());
    }
    for (auto e : geom->mesh.edges()) {
        // fprintf(file, "%.7f %.7f %.7f\n", solver->velocity[e].x(), solver->velocity[e].y(), solver->velocity[e].z());
        fprintf(file, "%.7f %.7f\n", solver->velocity[e].x(), solver->velocity[e].z());
    }
    for (auto v : geom->mesh.vertices()) {
        fprintf(file, "%.7g\n", solver->pressure[v]);
    }
    // for (auto v : geom->mesh.vertices()) {
    //     fprintf(file, "%.7g\n", solver->centripetal[v]);
    // }
    // for (auto e : geom->mesh.edges()) {
    //     fprintf(file, "%.7g\n", solver->centripetal[e]);
    // }

    fclose(file);
}
