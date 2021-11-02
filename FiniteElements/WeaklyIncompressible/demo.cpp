
// temporary
int FIGURE_MODE = 0;

enum MeshModes {
    MM_cylinder,
    MM_cylinder_2,
    MM_lid_driven_cavity,
    MM_lid_driven_cavity_2,
    MM_flow,
    MM_trivial_disk,
    MM_vector_poisson,
    NUM_MESH_MODES
};

struct Demo : public IBehaviour {
    Demo();
    void init();

    CameraController* controller;

    void update();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);

    void recreate_solver();

    SurfaceGeometry *geom;
    Solver *solver;

    int mesh_N;
    double mu;
    float theta0; // for cylinder flow demos.
    bool solving_while_moving; // for visualizing the flow around the object.
    bool many_sample_curve; // for getting a figure of the real curve

    // mesh generation options.
    int mesh_mode;
    bool random;

    bool flow_mode;
    vec2 obstruction_position;

    // Visualization.
    GLShaderProgram solution_shader; // Render the solution (velocity and pressure) to textures.
    GLuint solution_fbo; // For render-to-texture.
    GLuint solution_texture; // 3 components: r:velocity_x, g:velocity_y, b:pressure.
    GLuint solution_depth_texture; //--- Is this needed for completeness?
    // Plotting.
    GLShaderProgram sprite_shader;
    GLuint sprite_vao;
    bool wireframe;
    bool vector_field;
    bool show_div_P2; // plot the divergence in P2 space

    // Screenshot
    int screenshot_blx;
    int screenshot_bly;
    int screenshot_trx;
    int screenshot_try;
    float f_screenshot_blx;
    float f_screenshot_bly;
    float f_screenshot_trx;
    float f_screenshot_try;

    // Screenshots for specific figures.
    void take_high_res_screenshot(int n);
    int high_res_screenshot_n;

    float C; // Iteration parameter.
};

#include "WeaklyIncompressible/take_high_res_screenshots.cpp"


void Demo::recreate_solver()
{
    if (geom != nullptr) delete geom;


    if (mesh_mode == MM_cylinder) { // Flow around a cylinder.
        geom = square_minus_circle(0.28 / 3, theta0, 1, 1, mesh_N, false, vec2(0,0), many_sample_curve, 1,1/1.61803);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                // if (fabs(y) > 0.25) return vec2(0,0);
                if (x < -0.99) return vec2(1,0);
                if (x > 0.99) return vec2(1,0);
                return vec2(0,0);
            }
        );
    } else if (mesh_mode == MM_cylinder_2) {
        // geom = square_minus_circle(0.28 / 6, theta0, 1, 1, mesh_N, false, vec2(0,0), many_sample_curve, 1,0.25/1.61803);
        geom = square_minus_circle(0.28 / 6, theta0, 0.7, 2.2*1.5, mesh_N, true, vec2(0,0), many_sample_curve, 1,0.25/1.61803);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                if (x < -0.99) return vec2(0.25,0);
                if (x > 0.99) return vec2(0.25,0);
                return vec2(0,0);
            }
        );
    } else if (mesh_mode == MM_lid_driven_cavity) { // Lid-driven cavity.
        // geom = circle_mesh(mesh_N, false);
        auto sq_mesh = SquareMesh(mesh_N);
        geom = sq_mesh.geom;

        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);

        // Set the lid boundary condition explicitly, on the vertex and midpoint sample points.
        // (This is to avoid possible errors at corners if the boundary condition was specified with a function.)
        for (int i = 0; i < mesh_N+1; i++) {
            solver->u_boundary[sq_mesh.vertex(i,mesh_N)] = vec2(-1,0); // x is inverted...
        }
        for (int i = 0; i < mesh_N; i++) {
            auto v1 = sq_mesh.vertex(i,mesh_N);
            auto v2 = sq_mesh.vertex(i+1,mesh_N);
            solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(-1,0);
        }
        // for (int i = 0; i < mesh_N+1; i++) {
        //     solver->u_boundary[sq_mesh.vertex(i,0)] = vec2(-1,0); // x is inverted...
        // }
        // for (int i = 0; i < mesh_N; i++) {
        //     auto v1 = sq_mesh.vertex(i,0);
        //     auto v2 = sq_mesh.vertex(i+1,0);
        //     solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(-1,0);
        // }
    } else if (mesh_mode == MM_lid_driven_cavity_2) { // Lid-driven cavity with an obstruction.
        geom = square_minus_circle(0.25, theta0, 1, 1, mesh_N, false, obstruction_position, many_sample_curve);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                if (y > 0.99) return vec2(-3,0);
                // if (y < -0.99) return vec2(-3,0);
                return vec2(0,0);
            }
        );
    } else if (mesh_mode == MM_trivial_disk) { // Axis-aligned flow on the disk.
        geom = circle_mesh(mesh_N, random);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                return vec2(1, 0);
            }
        );
    } else if (mesh_mode == MM_vector_poisson) { // Testing the vector Poisson equation.
        geom = circle_mesh(mesh_N, random);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                return vec2(x*x - y*y, x*x*x - y*y*y);
            }
        );
    } else if (mesh_mode == MM_flow) {
        geom = circle_mesh(mesh_N, random);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        // solver->set_u_boundary(
        //     [](double x, double y)->vec2 {
        //         // if (fabs(y) < 0.2) return vec2(1,0);
        //         // return vec2(0,0);
        //         if (fabs(y) < 0.25) {
        //             return vec2(1,0);
        //         }
        //         if (fabs(x) < 0.25) {
        //             return vec2(0,1);
        //         }
        //         return vec2(0,0);
        //     }
        // );
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                float angles[6];
                float c = 2*M_PI * 30;
                for (int i = 0; i < 6; i++) {
                    angles[i] = c + i*2*M_PI/6;
                }
                float flows[6] = {
                    1,
                    -2,
                    1,
                    2,
                    -1,
                    -1,
                };
                float theta = atan2(y,x) + M_PI;
                for (int i = 0; i < 6; i++) {
                    if (fabs(c + theta - angles[i]) < (2*M_PI/6) / 4) return flows[i]*vec2(cos(angles[i]), sin(angles[i]));
                }
                return vec2(0,0);
            }
        );
        
        // geom = square_mesh_with_square_hole(mesh_N);
        // if (solver != nullptr) delete solver;
        // solver = new Solver(*geom, mu);
        // solver->set_u_boundary(
        //     [](double x, double y)->vec2 {
        //         if (x > 0.05 && x < 0.95 && y > 0.05 && y < 0.95) {
        //             return vec2(0,0);
        //         } else {
        //             return vec2(1,0);
        //         }
        //     }
        // );
    }
}


Demo::Demo()
{
    // Set solver parameters.
    mu = 1;
    flow_mode = false;
    C = 0.008;

    // Mesh generation options
    mesh_N = 4;
    mesh_mode = 0;
    random = false;
    theta0 = 1.2;
    solving_while_moving = false;
    obstruction_position = vec2(0,0);
    // Plotting options
    wireframe = false;
    vector_field = true;
    show_div_P2 = false;
    many_sample_curve = false;

    high_res_screenshot_n = 0;

    recreate_solver();

    // Visualization.
    solution_shader.add_shader(GLShader(VertexShader, SHADERS "solution/solution.vert"));
    solution_shader.add_shader(GLShader(TessControlShader, SHADERS "solution/solution.tcs"));
    solution_shader.add_shader(GLShader(TessEvaluationShader, SHADERS "solution/solution.tes"));
    solution_shader.add_shader(GLShader(FragmentShader, SHADERS "solution/solution.frag"));
    solution_shader.link();
    
    // Create framebuffer for render-to-texture.
    glGenFramebuffers(1, &solution_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, solution_fbo);
    glGenTextures(1, &solution_texture);
    glBindTexture(GL_TEXTURE_2D, solution_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 1024, 1024, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, solution_texture, 0);
    glGenTextures(1, &solution_depth_texture);
    glBindTexture(GL_TEXTURE_2D, solution_depth_texture);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, 1024, 1024, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, solution_depth_texture, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    
    // Plotting sprites.
    sprite_shader.add_shader(GLShader(VertexShader, SHADERS "plot/plot.vert"));
    sprite_shader.add_shader(GLShader(FragmentShader, SHADERS "plot/plot.frag"));
    sprite_shader.link();
    glGenVertexArrays(1, &sprite_vao);

    GLuint sprite_vbo;
    glGenBuffers(1, &sprite_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sprite_vbo);
    float sprite_uvs[2*4] = {
        0,0, 1,0, 1,1, 0,1
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4 * 2, sprite_uvs, GL_STATIC_DRAW);

    glGenVertexArrays(1, &sprite_vao);
    glBindVertexArray(sprite_vao);
    glBindBuffer(GL_ARRAY_BUFFER, sprite_vbo);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
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

}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_P) {
            mesh_N += 1;
            recreate_solver();
        }
        if (e.key.code == KEY_O) {
            mesh_N -= 1;
            if (mesh_N < 2) mesh_N = 2;
            recreate_solver();
        }
        if (e.key.code == KEY_R) {
            random = !random;
            recreate_solver();
        }
        if (e.key.code == KEY_1) {
            wireframe = !wireframe;
        }
        if (e.key.code == KEY_2) {
            vector_field = !vector_field;
        }
        if (e.key.code == KEY_5) {
            solver->div_P2_P2(solver->u, solver->div_u_P2);
            
            // test for zero divergence
            // P2Attachment<vec2> vf(solver->geom.mesh);
            // for (auto v : solver->geom.mesh.vertices()) vf[v] = vec2(1,0);
            // for (auto e : solver->geom.mesh.edges()) vf[e] = vec2(1,0);
            // solver->div_P2_P2(vf, solver->div_u_P2);
        }
        if (e.key.code == KEY_6) {
            show_div_P2 = !show_div_P2;
        }
        if (e.key.code == KEY_M) {
            mesh_mode = (mesh_mode + 1) % NUM_MESH_MODES;
            recreate_solver();
        }
        if (e.key.code == KEY_N) {
            mesh_mode -= 1;
            if (mesh_mode < 0) mesh_mode = NUM_MESH_MODES-1;
            recreate_solver();
        }
        if (e.key.code == KEY_9) {
            mesh_N = 50;
            recreate_solver();
        }
        // if (e.key.code == KEY_I) {
        //     solver->C = C;
        //     solver->iterate();
        // }
        if (e.key.code == KEY_7) {
            solver->solve_taylor_hood();
            solver->pressure_update(true);
        }
        if (e.key.code == KEY_Y) {
            solver->project_divergence();
        }
        if (e.key.code == KEY_U) {
            flow_mode = !flow_mode;
        }
        if (e.key.code == KEY_8) {
            solver->write_sparsity_pattern = true;
        }
        if (e.key.code == KEY_Z) {
            solving_while_moving = !solving_while_moving;
        }
        if (e.key.code == KEY_0) {
            many_sample_curve = !many_sample_curve;
            recreate_solver();
        }
        if (e.key.code == KEY_G) { // temporary, figure creation
            FIGURE_MODE = (FIGURE_MODE + 1) % 2;
        }
        if (e.key.code == KEY_I) {
            take_high_res_screenshot(high_res_screenshot_n);
            high_res_screenshot_n = (high_res_screenshot_n+1)%4;
        }
        static int counter = 0;
        std::string pre = DATA + ("stokes_" + std::to_string(mesh_mode) + "_" + (wireframe ? "wireframe" : "velocity") + "_" + std::to_string(counter));
        if (e.key.code == KEY_T) {
            world->graphics.screenshot(pre + ".ppm",
        			      screenshot_blx, screenshot_bly, screenshot_trx - screenshot_blx, screenshot_try - screenshot_bly);
            counter += 1;
        }
    }
}

void Demo::update()
{
    if (world->input.keyboard.down(KEY_G)) {
        // Draw the screenshot rectangle.
        std::vector<vec2> ps = {
            vec2(f_screenshot_blx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_bly)
        };
        world->graphics.paint.chain_2D(ps, 1, vec4(1,0,0,1));
    }
    const float sp = 3;
    if (world->input.keyboard.down(KEY_3)) C *= 1 - sp*dt;
    if (world->input.keyboard.down(KEY_4)) C *= 1 + sp*dt;

    if (world->input.keyboard.down(KEY_X)) {
        theta0 -= dt;
        recreate_solver();
        if (solving_while_moving) {
            solver->solve_taylor_hood();
            solver->pressure_update(true);
        }
    }
    if (world->input.keyboard.down(KEY_C)) {
        theta0 += dt;
        recreate_solver();
        if (solving_while_moving) {
            solver->solve_taylor_hood();
            solver->pressure_update(true);
        }
    }
    const float ob_speed = 1.f;
    if (world->input.keyboard.down(KEY_LEFT_ARROW)) obstruction_position.x() += ob_speed * dt; // swapped
    if (world->input.keyboard.down(KEY_RIGHT_ARROW)) obstruction_position.x() -= ob_speed * dt;
    if (world->input.keyboard.down(KEY_DOWN_ARROW)) obstruction_position.y() -= ob_speed * dt;
    if (world->input.keyboard.down(KEY_UP_ARROW)) obstruction_position.y() += ob_speed * dt;

    if (world->input.keyboard.down(KEY_LEFT_ARROW)
            || world->input.keyboard.down(KEY_UP_ARROW)
            || world->input.keyboard.down(KEY_DOWN_ARROW)
            || world->input.keyboard.down(KEY_RIGHT_ARROW)) {
        recreate_solver();
        if (solving_while_moving) {
            solver->solve_taylor_hood();
            solver->pressure_update(true);
        }
    }
}


void Demo::post_render_update()
{

    // // Draw boundary velocity.
    // for (auto v : geom->mesh.vertices()) {
    //     vec3 pos = eigen_to_vec3(geom->position[v]);
    //     vec2 vec = solver->u_boundary[v];
    // }
    // for (auto e : geom->mesh.edges()) {
    //     vec3 pos = eigen_to_vec3(solver->midpoints[e]);
    //     vec2 vec = solver->u_boundary[e];
    //     world->graphics.paint.line(pos, pos + 0.06*vec3(vec.x(), 0, vec.y()), 0.01, vec4(1,0,0,1));
    // }
    
    if (flow_mode) {
        solver->C = (0.05/0.008)*C*dt;
        solver->iterate();
    }

    // Scale the pressure.
    VertexAttachment<double> scaled_pressure(geom->mesh);
    double min_pressure = std::numeric_limits<double>::infinity();
    double max_pressure = 0;
    for (auto v : geom->mesh.vertices()) {
        if (solver->p[v] < min_pressure) {
            min_pressure = solver->p[v];
        }
        if (solver->p[v] > max_pressure) {
            max_pressure = solver->p[v];
        }
    }
    if (max_pressure != min_pressure) {
        for (auto v : geom->mesh.vertices()) {
            scaled_pressure[v] = (solver->p[v] - min_pressure) / (max_pressure - min_pressure);
        }
    } else {
        for (auto v : geom->mesh.vertices()) {
            scaled_pressure[v] = 0;
        }
    }

    // Render the solution textures.
    std::vector<vec2> position_data;
    std::vector<vec2> velocity_data;
    std::vector<float> pressure_data;
    std::vector<float> div_u_data;
    std::vector<float> div_u_P2_data;
    // Create a 6-vertex patch per triangle.
    for (auto tri : geom->mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        do {
            auto v = he.vertex();
            auto e = he.edge();
            for (auto pos : {geom->position[v], solver->midpoints[e]}) {
                position_data.push_back(vec2(pos.x(), pos.z()));
            }
            velocity_data.push_back(solver->u[v]);
            velocity_data.push_back(solver->u[e]);
            pressure_data.push_back(scaled_pressure[v]);
            pressure_data.push_back(0.f); // Dummy data, as there is no midpoint pressure.
            div_u_data.push_back(solver->div_u[v]);
            div_u_data.push_back(0.f); // Dummy data, as there is no midpoint div(u).
            div_u_P2_data.push_back(solver->div_u_P2[v]);
            div_u_P2_data.push_back(solver->div_u_P2[e]);

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
    GLuint vbos[5]; // position, velocity, pressure
    glGenBuffers(5, vbos);
    struct {
        const void *data;
        size_t data_size;
        size_t gl_data_number;
        GLenum gl_data_type;
    } data_to_upload[5] = {
        {&position_data[0], sizeof(vec2), 2, GL_FLOAT}, // position
        {&velocity_data[0], sizeof(vec2), 2, GL_FLOAT}, // velocity (P2)
        {&pressure_data[0], sizeof(float), 1, GL_FLOAT}, // pressure (P1)
        {&div_u_data[0], sizeof(float), 1, GL_FLOAT}, // divergence of velocity (P1)
        {&div_u_P2_data[0], sizeof(float), 1, GL_FLOAT} // divergence of velocity (P2)
    };
    // Upload the data.
    for (int i = 0; i < 5; i++) {
        auto metadata = data_to_upload[i];
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        glBufferData(GL_ARRAY_BUFFER, data_num_vertices * metadata.data_size, metadata.data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(i, metadata.gl_data_number, metadata.gl_data_type, GL_FALSE, 0, (const void *) 0);
        glEnableVertexAttribArray(i);
    }

    // Render to texture.
    solution_shader.bind();
    glUniform1i(solution_shader.uniform_location("show_div_P2"), show_div_P2 ? 1 : 0);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport); // save the viewport to restore after rendering sample geometry.
    glDisable(GL_SCISSOR_TEST);
    glViewport(0,0,1024,1024);
    glBindFramebuffer(GL_FRAMEBUFFER, solution_fbo);
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    glPatchParameteri(GL_PATCH_VERTICES, 6);
    glDrawArrays(GL_PATCHES, 0, data_num_vertices);
    glEnable(GL_DEPTH_TEST);

    glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);
    glEnable(GL_SCISSOR_TEST);
    glEnable(GL_DEPTH_TEST);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    solution_shader.unbind();

    // Clean up.
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(3, vbos);
    
    // Draw velocity field.
    const float velocity_mul = 0.12;
    if (vector_field) {
        glBindFramebuffer(GL_FRAMEBUFFER, solution_fbo);
        auto solution_pixels = std::vector<float>(4*1024*1024);
        glReadPixels(0,0,1024,1024, GL_RGBA, GL_FLOAT, &solution_pixels[0]);
        glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);
        std::vector<float> pressures;
        // const int skip = 20;
        // const int skip = 14;
        const int skip = 35;
	std::vector<vec2> positions_2D;
        for (int i = 0; i < 1024; i += skip) {
            float x = -1 + i*2.f/(1024-1.f);
            for (int j = 0; j < 1024; j += skip) {
                float y = -1 + j*2.f/(1024-1.f);
                float velocity_x = solution_pixels[4*(1024*j + i) + 0];
                float velocity_y = solution_pixels[4*(1024*j + i) + 1];
                float pressure = solution_pixels[4*(1024*j + i) + 2];
                const float thickness = 0.005;
                const vec4 color = vec4(0,0,0,1);
                const float epsilon = 1e-5;
                if (fabs(velocity_x) < epsilon && fabs(velocity_y) < epsilon) {
                } else {
                    // world->graphics.paint.sphere(vec3(x, 0.005, y), 0.003, vec4(0,0,0,1));
                    world->graphics.paint.line(vec3(x, 0.005, y), vec3(x + velocity_mul*velocity_x, 0, y + velocity_mul*velocity_y), thickness, color);
	            pressures.push_back(pressure);
                    positions_2D.push_back(vec2(x,y));
                }
            }
        }
        const float circle_wid = 0.006;
        float x_half_wid = 1.25;
        float extent = 1.05*x_half_wid/4;
        float w = extent;
        float h = 0.566*extent;
        vec2 screenshot_positions[4] = {
            vec2(5*x_half_wid/8, 0),
            vec2(x_half_wid/4, 0),
            vec2(-x_half_wid/4, 0),
            vec2(-5*x_half_wid/8, 0)
        };
        for (int i = 0; i < positions_2D.size(); i++) {
            positions_2D[i].x() = 0.5*(positions_2D[i].x() - screenshot_positions[(high_res_screenshot_n+3)%4].x())/w + 0.5;
            positions_2D[i].y() = 1 - (0.5*(positions_2D[i].y() - screenshot_positions[(high_res_screenshot_n+3)%4].y())/h + 0.5);
        }
        
        #if 0
        std::vector<vec2> pos(1);
        for (int i = 0; i < positions_2D.size(); i++) {
            pos[0] = positions_2D[i] - vec2(0,0.0038);
            world->graphics.paint.circles(main_camera, pos, circle_wid, vec4(pressures[i],pressures[i],pressures[i],1), 0.3, vec4(0,0,0,1));
        }
        #endif
    }

    // Visual helper lines.
    // world->graphics.paint.line(vec3(-1,0,1+0.05), vec3(1,0,1+0.05), 0.01, vec4(0,0,1,1));

    // Draw the solution texture.
    const float asp = 0.566;
    const float height = 1.f/4.f;
    const float width = height * asp;
    struct {
        float bl_x;
        float bl_y;
        float tr_x;
        float tr_y;
    } img_datas[4] = {
        {0,0, height*asp,height},
        {0,height, height*asp,2*height},
        {0,2*height, height*asp,3*height},
        {0,3*height, height*asp,4*height},

        // {0,0, width,height},
        // {width,0, 2*width,height},
        // {2*width,0, 3*width,height},
        // {3*width,0, 4*width,height},
    };
    
    // Show solution sprites.
#if 0
    sprite_shader.bind();
    // glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, solution_texture);
    glUniform1i(sprite_shader.uniform_location("tex"), 0);
    glBindVertexArray(sprite_vao);
    for (int i = 0; i < 4; i++) {
        auto img_data = img_datas[i];
        glUniform2f(sprite_shader.uniform_location("bottom_left"), img_data.bl_x,img_data.bl_y);
        glUniform2f(sprite_shader.uniform_location("top_right"), img_data.tr_x, img_data.tr_y);
        glUniform1i(sprite_shader.uniform_location("mode"), i);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    }
    // glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    sprite_shader.unbind();
#endif

    // Draw wireframe or base.
    if (wireframe) {
        float thickness = 0.00125;
        world->graphics.paint.wireframe(*geom, mat4x4::translation(0,0,0), thickness);
        for (auto v : geom->mesh.vertices()) {
            // world->graphics.paint.sphere(eigen_to_vec3(geom->position[v]), 0.0075, vec4(0.9,0.9,0.9,1));
        }
        for (auto e : geom->mesh.edges()) {
            // world->graphics.paint.sphere(eigen_to_vec3(solver->midpoints[e]), 0.0075, vec4(0.9,0.9,0.9,1));
        }
    }
        // Draw boundaries.
        for (auto start : geom->mesh.boundary_loops()) {
            auto ps = std::vector<vec3>();
            auto he = start;
            do {
                ps.push_back(vec3(0,0,0)+eigen_to_vec3(geom->position[he.vertex()]));
                he = he.next();
            } while (he != start);
            if (wireframe) ps.push_back(ps[0]);
            else ps.push_back(ps[0]);
            world->graphics.paint.chain(ps, 0.001, vec4(0,0,0,1));
        }
    // Draw boundary condition.
    if (many_sample_curve) { //drawing specific figure
        for (auto start : geom->mesh.boundary_loops()) {
            auto ps = std::vector<vec3>();
            auto he = start;
            do {
                auto v = he.vertex();
                auto e = he.edge();
                vec3 bv[2] = {
                    vec3(solver->u_boundary[v].x(), 0, solver->u_boundary[v].y()),
                    vec3(solver->u_boundary[e].x(), 0, solver->u_boundary[e].y()) };
                vec3 pos[2] = {eigen_to_vec3(geom->position[v]),
                            eigen_to_vec3(solver->midpoints[e])};
                vec3 shift = vec3(0,0,0);
                const float epsilon = 1e-5;
                for (int i = 0; i < 2; i++) {
                    if (bv[i].length() > epsilon) {
                        // world->graphics.paint.sphere(pos[i]+shift, 0.003, vec4(1,0.6,0.6,1));
                        const float line_wid = 0.005;
                        vec4 col = vec4(0.7,0.9,1,1);
                        world->graphics.paint.line(pos[i]+shift, shift+pos[i] + bv[i]*velocity_mul, line_wid, col);
                        // arrow head
                        const float arrow_wid = 0.015;
                        vec3 tip = 1.0001*shift+pos[i]+bv[i]*velocity_mul;
                        vec3 bv_perp = vec3(-bv[i].z(), 0, bv[i].x());
                        vec3 arrow_bit_1 = tip - arrow_wid*bv[i].normalized() + arrow_wid*bv_perp.normalized();
                        vec3 arrow_bit_2 = tip - arrow_wid*bv[i].normalized() - arrow_wid*bv_perp.normalized();
                        world->graphics.paint.line(tip, arrow_bit_1, line_wid, col);
                        world->graphics.paint.line(tip, arrow_bit_2, line_wid, col);
                    }
                }
                he = he.next();
            } while (he != start);
        }
    }

    SurfaceMesh tri_mesh;
    SurfaceGeometry tri_geom(tri_mesh);
#if 0
    auto v = tri_mesh.add_vertex();
    auto vp = tri_mesh.add_vertex();
    auto vpp = tri_mesh.add_vertex();
    tri_mesh.add_triangle(v, vp, vpp);
    vec3 c = vec3(-3,0,0);
    vec3 p1 = c + vec3(0,0,0);
    vec3 p2 = c + vec3(1,0,0);
    vec3 p3 = c + vec3(1.25,0,2);
    tri_geom.position[v] = vec3_to_eigen(p1);
    tri_geom.position[vp] = vec3_to_eigen(p2);
    tri_geom.position[vpp] = vec3_to_eigen(p3);
    world->graphics.paint.wireframe(tri_geom, mat4x4::translation(0,0,0), 0.001);
    tri_geom.position[vp] += Eigen::Vector3f(0,1,0);
    world->graphics.paint.wireframe(tri_geom, mat4x4::translation(0,0,0), 0.005);
    world->graphics.paint.line(
            p2,
            p2 + vec3(0,1,0),
            0.002,
            vec4(0,0,0,1));
    vec3 drop = p1 + (p3 - p1).normalized() * vec3::dot(p2 - p1, (p3 - p1).normalized());
    world->graphics.paint.line(
            drop,
            p2,
            0.01,
            vec4(1,0,0,1));
    vec3 corner1 = p1 + p2 - drop;
    vec3 corner2 = p3 + p2 - drop;
    world->graphics.paint.line(
            p1,
            corner1,
            0.01,
            vec4(1,0,0,1));
    world->graphics.paint.line(
            p3,
            corner2,
            0.01,
            vec4(1,0,0,1));
    world->graphics.paint.line(
            p1,
            p3,
            0.01,
            vec4(1,0,0,1));
    world->graphics.paint.line(
            corner1,
            corner2,
            0.01,
            vec4(1,0,0,1));
#endif
    

#if 0
    auto v = tri_mesh.add_vertex();
    vec3 c = vec3(-3,0,0);
    tri_geom.position[v] = vec3_to_eigen(c);
    #define AMOUNT 4 
    // float angles[AMOUNT] = {0, 0.1, 0.3, 0.4, 0.7};
    float angles[AMOUNT] = {0, 0.15, 0.35, 0.7};
    Vertex vs[AMOUNT];
    vec3 vs_p[AMOUNT];
    for (int i = 0; i < AMOUNT; i++) {
        float theta = 2*M_PI*angles[i];
        float co = cos(theta);
        float si = sin(theta);
        vs[i] = tri_mesh.add_vertex();
        vec3 p = c + vec3(co, 0, si);
        vs_p[i] = p;
        tri_geom.position[vs[i]] = vec3_to_eigen(p);
    }
    for (int i = 0; i < AMOUNT; i++) {
        tri_mesh.add_triangle(v, vs[i], vs[(i+1)%AMOUNT]);
    }
    world->graphics.paint.wireframe(tri_geom, mat4x4::translation(0,0,0), 0.01);

    for (int i = 0; i < AMOUNT; i++) {
        int ii = (i+1)%AMOUNT;
        vec3 m1 = 0.5*c + 0.5*vs_p[i];
        vec3 m2 = 0.5*c + 0.5*vs_p[ii];
        auto line = [&](vec3 a, vec3 b) {
            vec3 shift = vec3(0,0.01,0);
            world->graphics.paint.line(a+shift, b+shift, 0.008, vec4(0,0,0,1));
        };
        line(m1, m2);

        auto circumcenter = [&](vec3 A, vec3 B, vec3 C) {
            // https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
            Eigen::Vector3f a = vec3_to_eigen(A);
            Eigen::Vector3f b = vec3_to_eigen(B);
            Eigen::Vector3f c = vec3_to_eigen(C);
            Eigen::Vector3f ac = c - a ;
            Eigen::Vector3f ab = b - a ;
            Eigen::Vector3f abXac = ab.cross(ac);

            // this is the vector from a TO the circumsphere center
            Eigen::Vector3f toCircumsphereCenter = (abXac.cross( ab )*ac.dot(ac) + ac.cross( abXac )*ab.dot(ab)) / (2.f*abXac.dot(abXac)) ;
            float circumsphereRadius = toCircumsphereCenter.norm();
            // The 3 space coords of the circumsphere center then:
            Eigen::Vector3f ccs = a  +  toCircumsphereCenter ; // now this is the actual 3space location
            return eigen_to_vec3(ccs);
        };
        vec3 cc = circumcenter(c, vs_p[i], vs_p[ii]);
        world->graphics.paint.sphere(cc, 0.03, vec4(1,0,0,1));
        line(c, m1);
        line(c, m2);
        if (vec2::dot(vec2(cc.x() - vs_p[i].x(), cc.z() - vs_p[i].z()),
                vec2(vs_p[ii].x() - vs_p[i].x(), vs_p[ii].z() - vs_p[i].z()).perp()) < 0) {
            vec3 m = 0.5*vs_p[i] + 0.5*vs_p[ii];
            line(m, cc);
            line(m1, m);
            line(m2, m);
        } else {
            line(m1, cc);
            line(m2, cc);
        }
    }
#endif

#if 1
    for (int R = 0; R < 3; R++) {
        auto v = tri_mesh.add_vertex();
        auto vp = tri_mesh.add_vertex();
        auto vpp = tri_mesh.add_vertex();
        Vertex vs[3] = {v,vp,vpp};
        tri_mesh.add_triangle(vp, v, vpp);
        vec3 c = vec3(-3-1.6*R,0,0);

        vec3 p1;
        vec3 p2;
        vec3 p3;
        if (FIGURE_MODE == 0) {
            p1 = c + vec3(0,0,0);
            p2 = c + vec3(-1.25,0,-0.6);
            p3 = c + vec3(0.66,0,1.7);
        } else {
            p1 = c + vec3(0,0,0);
            p2 = c + vec3(-1,0,0);
            p3 = c + vec3(0,0,1);
        }
        vec3 ps[3] = {p1,p2,p3};
        tri_geom.position[v] = vec3_to_eigen(p1);
        tri_geom.position[vp] = vec3_to_eigen(p2);
        tri_geom.position[vpp] = vec3_to_eigen(p3);
        #if 0
        world->graphics.paint.wireframe(tri_geom, mat4x4::translation(0,0,0), 0.001);
        tri_geom.position[vs[R]] += Eigen::Vector3f(0,1,0);
        world->graphics.paint.wireframe(tri_geom, mat4x4::translation(0,0,0), 0.005);
        world->graphics.paint.line(
                ps[R],
                ps[R] + vec3(0,1,0),
                0.002,
                vec4(0,0,0,1));
        #else
        world->graphics.paint.line(ps[0], ps[1], 0.005, vec4(0,0,0,1));
        world->graphics.paint.line(ps[1], ps[2], 0.005, vec4(0,0,0,1));
        world->graphics.paint.line(ps[2], ps[0], 0.005, vec4(0,0,0,1));
        
        #endif
    }
#endif
}



void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
        if (e.button.code == MOUSE_LEFT) {
            screenshot_blx = int(world->graphics.window_viewport.w * e.cursor.x);
            screenshot_bly = int(world->graphics.window_viewport.h * e.cursor.y);
            f_screenshot_blx = e.cursor.x;
            f_screenshot_bly = e.cursor.y;
        }
        if (e.button.code == MOUSE_RIGHT) {
            screenshot_trx = int(world->graphics.window_viewport.w * e.cursor.x);
            screenshot_try = int(world->graphics.window_viewport.h * e.cursor.y);
            f_screenshot_trx = e.cursor.x;
            f_screenshot_try = e.cursor.y;
        }
    }
}
