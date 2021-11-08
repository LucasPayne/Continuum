#include "core.h"
#include "SurfaceNavierStokes/demo.h"
#include "mesh_generators.cpp"
#include "mesh_processing/extensions/assimp_convert.h"

bool filming = false;
vec3 shift = vec3(0, 0.001, 0);
int render_solution_mode = 0;

Face test_face;
double test_shift_t = 0.;

Demo::Demo()
{
    solver = nullptr;
    geom = nullptr;
}

std::tuple<Face, vec3, mat3x3> Demo::traverse(Face tri, vec3 origin, vec3 shift, mat3x3 destination_to_origin_matrix, int ignore_index)
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

    // world->graphics.paint.line(origin, origin + vec3(0,1,0), 0.005, vec4(0,0,0,1));
    // world->graphics.paint.line(origin+vec3(0,0.001,0), origin+e1+vec3(0,0.001,0), 0.01, vec4(1,1,0,1));
    // world->graphics.paint.line(origin+vec3(0,0.001,0), origin+e2+vec3(0,0.001,0), 0.01, vec4(0,1,0,1));
    // if (ignore_index != -1) world->graphics.paint.line(ps[ignore_index], ps[(ignore_index+1)%3], 0.01, vec4(0,1,1,1));

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

        vec3 to = origin+shift;

        vec3 tri_normal = solver->triangle_normal[tri];
        world->graphics.paint.line(origin+0.001*tri_normal, to+0.001*tri_normal, 0.01, vec4(0,1,0,1));

        return {tri, to, destination_to_origin_matrix};
    }
    Halfedge hit_he = hes[line_hit_index];
    if (hit_he.twin().face().null()) {
        // Hit the boundary. Stop at the boundary intersection.

        vec3 to = origin + min_t * shift;

        vec3 tri_normal = solver->triangle_normal[tri];
        world->graphics.paint.line(origin+0.001*tri_normal, to+0.001*tri_normal, 0.01, vec4(0,1,0,1));

        return {tri, to, destination_to_origin_matrix};
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
    // Rotate the rest of the shift to this new triangle plane.
    vec3 new_shift = E1*vec3::dot((1-min_t)*shift, E1) + to_E2*vec3::dot((1-min_t)*shift, from_E2);
    
    // Construct a rotation matrix which transforms vectors from the new face to the old face.
    float d = vec3::dot(from_E2, to_E2);
    if (d < 0) d = 0;
    if (d > 1) d = 1;
    float theta = -acos(d);
    mat3x3 to_matrix = Quaternion::from_axis_angle(E1, theta).matrix().top_left();

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

    vec3 to = origin + min_t*shift;

    vec3 tri_normal = solver->triangle_normal[tri];
    world->graphics.paint.line(origin+0.001*tri_normal, to+0.001*tri_normal, 0.01, vec4(0,1,0,1));

    mat3x3 matrix = destination_to_origin_matrix * to_matrix;

    return traverse(to_face, to, new_shift, matrix, to_ignore_index);
}

void Demo::init()
{
    // Create a camera controller.
    auto cameraman = world->entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    camera->background_color = vec4(1,1,1,1);
    auto t = cameraman.add<Transform>(0,2,0);
    main_camera = camera;
    controller = world->add<CameraController>(cameraman);
    controller->angle = -M_PI/2;
    controller->azimuth = M_PI;

    geom = assimp_to_surface_geometry(std::string(MODELS) + "sphere_fine.off");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "icosahedron.ply");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "simple_gear.ply");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "side_gear.stl");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "drive_wheel.stl");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "bunny_head.stl");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "tangram.stl");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "cylinder.stl");
    Eigen::Vector3f avg = Eigen::Vector3f(0,0,0);
    for (auto v : geom->mesh.vertices()) {
        // std::cout << geom->position[v] << "\n";
        // getchar();
        avg += geom->position[v];
    }
    avg /= geom->mesh.num_vertices();
    double max_norm = 0;
    for (auto v : geom->mesh.vertices()) {
        geom->position[v] -= avg;
        if (geom->position[v].norm() > max_norm) {
            max_norm = geom->position[v].norm();
        }
    }
    for (auto v : geom->mesh.vertices()) {
        geom->position[v] /= max_norm;
    }
    
    geom->mesh.lock();
    
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

        // geom = square_mesh(25);
        // geom = torus_mesh(25);

    // for (auto v : geom->mesh.vertices()) {
    //     vec3 p = eigen_to_vec3(geom->position[v]);
    //     float theta = 0.95*(p.x()+1)*M_PI;
    //     geom->position[v] = Eigen::Vector3f(cos(theta), sin(theta), p.z());
    //     // geom->position[v] += Eigen::Vector3f(0, 0.1*sin(4*p.x()), 0);
    // }
    
    double viscosity = 0.001;
    solver = new SurfaceNavierStokesSolver(*geom, viscosity);

    solver->set_source([&](double x, double y, double z)->vec3 {
        //const double r = 0.125;
        // if ((vec2(x,z) - vec2(-0.85,0)).length() <= r) {
        //     return vec3(3000, 0, 0);
        // }
        const double r = 0.25;
        vec3 s = vec3(0,0,0);
        if (y > 0) {
            if ((vec2(x,z)).length() <= r) {
                s = vec3(3000, 0, 0);
            }
            // if ((vec2(x,z) - vec2(-0.3, 0)).length() <= r) {
            //     s = vec3(-3000, 0, 0);
            // }
        }
        return s;
        // return s + vec3(0,-2000,0);
        
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

    // Shaders
    //------------------------------------------------------------
    solution_shader.add_shader(GLShader(VertexShader, SHADERS "surface_flow_visualization_solution/solution.vert"));
    solution_shader.add_shader(GLShader(TessControlShader, SHADERS "surface_flow_visualization_solution/solution.tcs"));
    solution_shader.add_shader(GLShader(TessEvaluationShader, SHADERS "surface_flow_visualization_solution/solution.tes"));
    solution_shader.add_shader(GLShader(FragmentShader, SHADERS "surface_flow_visualization_solution/solution.frag"));
    solution_shader.link();
}

void Demo::render_solution()
{
    const double pressure_mul = 1.;

    // Scale the pressure.
    VertexAttachment<double> scaled_pressure(geom->mesh);
    double min_pressure = std::numeric_limits<double>::infinity();
    double max_pressure = 0;
    for (auto v : geom->mesh.vertices()) {
        if (solver->pressure[v] < min_pressure) {
            min_pressure = solver->pressure[v];
        }
        if (solver->pressure[v] > max_pressure) {
            max_pressure = solver->pressure[v];
        }
    }
    if (max_pressure != min_pressure) {
        for (auto v : geom->mesh.vertices()) {
            scaled_pressure[v] = pressure_mul*(solver->pressure[v] - min_pressure) / (max_pressure - min_pressure);
        }
    } else {
        for (auto v : geom->mesh.vertices()) {
            scaled_pressure[v] = 0;
        }
    }

    // Render the solution textures.
    auto position_data = std::vector<vec3>();
    auto velocity_data = std::vector<vec3>();
    auto pressure_data = std::vector<float>();
    // Create a 6-vertex patch per triangle.
    for (auto tri : geom->mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        do {
            auto v = he.vertex();
            auto e = he.edge();
            Eigen::Vector3f midpoint = 0.5*geom->position[e.a().vertex()] + 0.5*geom->position[e.b().vertex()];
            for (auto pos : {geom->position[v], midpoint}) {
                position_data.push_back(eigen_to_vec3(pos));
            }
            velocity_data.push_back(solver->velocity[v]);
            velocity_data.push_back(solver->velocity[e]);
            pressure_data.push_back(scaled_pressure[v]);
            pressure_data.push_back(0.f); // Dummy data, as there is no midpoint pressure.

            he = he.next();
        } while (he != start);
    }
    // Check that the right number of vertices were created.
    size_t data_num_vertices = 6*geom->mesh.num_faces();
    for (size_t len : {position_data.size(), velocity_data.size(), pressure_data.size()}) {
        assert(len == data_num_vertices);
    }

    GLuint vao;
    glCreateVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint vbos[3]; // position, velocity, pressure
    glGenBuffers(3, vbos);
    struct {
        const void *data;
        size_t data_size;
        size_t gl_data_number;
        GLenum gl_data_type;
    } data_to_upload[3] = {
        {&position_data[0], sizeof(vec3), 3, GL_FLOAT}, // position
        {&velocity_data[0], sizeof(vec3), 3, GL_FLOAT}, // velocity (P2)
        {&pressure_data[0], sizeof(float), 1, GL_FLOAT}, // pressure (P1)
    };
    // Upload the data.
    for (int i = 0; i < 3; i++) {
        auto metadata = data_to_upload[i];
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        glBufferData(GL_ARRAY_BUFFER, data_num_vertices * metadata.data_size, metadata.data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(i, metadata.gl_data_number, metadata.gl_data_type, GL_FALSE, 0, (const void *) 0);
        glEnableVertexAttribArray(i);
    }

    // Render the mesh.
    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT); //------
    solution_shader.bind();
    glUniform1i(solution_shader.uniform_location("mode"), render_solution_mode);
    auto vp_matrix = main_camera->view_projection_matrix();
    glUniformMatrix4fv(solution_shader.uniform_location("mvp_matrix"), 1, GL_FALSE, (const GLfloat *) &vp_matrix);
    glPatchParameteri(GL_PATCH_VERTICES, 6);
    glDrawArrays(GL_PATCHES, 0, data_num_vertices);
    solution_shader.unbind();

    // Clean up.
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(3, vbos);
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
        if (e.key.code == KEY_6) {
            solver->set_velocity([&](double x, double y, double z)->vec3 {
                double r = 0.50;
                if ((x)*(x) + z*z <= r*r) return 6.7745353*vec3(1,0,1);
                return vec3(0,0,0);
            });
        }
        if (e.key.code == KEY_1) render_solution_mode = 0;
        if (e.key.code == KEY_2) render_solution_mode = 1;
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

    #if 0
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
        save_solution(std::string(DATA) + "flows/icosahedron/surface_navier_stokes_" + std::to_string(counter) + ".txt");
        counter += 1;
    }

    #if 1
    const float test_move_speed = 0.25;
    if (world->input.keyboard.down(KEY_LEFT_ARROW)) {
        test_shift_t += test_move_speed * dt;
    }
    if (world->input.keyboard.down(KEY_RIGHT_ARROW)) {
        test_shift_t -= test_move_speed * dt;
    }
    vec3 c = eigen_to_vec3(geom->barycenter(test_face));
    vec3 s = solver->triangle_projection_matrix[test_face]*vec3(-10*test_shift_t, 0, test_shift_t);
    world->graphics.paint.sphere(c, 0.01, vec4(0.5,0.5,0.5,1));
    world->graphics.paint.line(c + shift, c+s + shift, 0.005,  vec4(0.5,1,0.5,1));

    Face out_face;
    vec3 out_pos;
    mat3x3 out_matrix;
    std::tie(out_face, out_pos, out_matrix) = traverse(test_face, c, s, mat3x3::identity());
    world->graphics.paint.sphere(out_pos, 0.03, vec4(1,0.2,0.2,1));
    world->graphics.paint.sphere(eigen_to_vec3(geom->barycenter(out_face)), 0.02, vec4(0.5,1,0.2,1));
    vec3 transformed = out_matrix * eigen_to_vec3(geom->vector(out_face.halfedge()));
    world->graphics.paint.line(c+shift, c+shift + transformed, 0.025, vec4(1,0,1,1));
    #endif
}


void Demo::post_render_update()
{
    render_solution();
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
        fprintf(file, "%.7f %.7f %.7f\n", solver->velocity[v].x(), solver->velocity[v].y(), solver->velocity[v].z());
        // fprintf(file, "%.7f %.7f\n", solver->velocity[v].x(), solver->velocity[v].z());
    }
    for (auto e : geom->mesh.edges()) {
        fprintf(file, "%.7f %.7f %.7f\n", solver->velocity[e].x(), solver->velocity[e].y(), solver->velocity[e].z());
        // fprintf(file, "%.7f %.7f\n", solver->velocity[e].x(), solver->velocity[e].z());
    }
    for (auto v : geom->mesh.vertices()) {
        fprintf(file, "%.7g\n", solver->pressure[v]);
    }
    for (auto v : geom->mesh.vertices()) {
        fprintf(file, "%.7g\n", solver->centripetal[v]);
    }
    for (auto e : geom->mesh.edges()) {
        fprintf(file, "%.7g\n", solver->centripetal[e]);
    }

    fclose(file);
}
