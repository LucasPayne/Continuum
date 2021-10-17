

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

    float C; // Iteration parameter.
};


void Demo::recreate_solver()
{
    if (geom != nullptr) delete geom;


    if (mesh_mode == MM_cylinder) { // Flow around a cylinder.
        geom = square_minus_circle(0.28, theta0, 1, 1, mesh_N);
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
        geom = square_minus_circle(0.28, theta0, 0.3, 2.5, mesh_N, true);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                if (x < -0.99) return vec2(1,0);
                if (x > 0.99) return vec2(1,0);
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
            solver->u_boundary[sq_mesh.vertex(i,mesh_N)] = vec2(1,0); // x is inverted...
        }
        for (int i = 0; i < mesh_N; i++) {
            auto v1 = sq_mesh.vertex(i,mesh_N);
            auto v2 = sq_mesh.vertex(i+1,mesh_N);
            solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(1,0);
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
        geom = square_minus_circle(0.12, theta0, 1, 1, mesh_N, false, obstruction_position);
        if (solver != nullptr) delete solver;
        solver = new Solver(*geom, mu);
        solver->set_u_boundary(
            [](double x, double y)->vec2 {
                if (y > 0.99) return vec2(3,0);
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
            if (mesh_mode < 0) mesh_mode = 0;
            recreate_solver();
        }
        if (e.key.code == KEY_9) {
            mesh_N = 50;
            recreate_solver();
        }
        if (e.key.code == KEY_I) {
            solver->C = C;
            solver->iterate();
        }
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
        for (int i = 0; i < 1024; i += 28) {
            float x = -1 + i*2.f/(1024-1.f);
            for (int j = 0; j < 1024; j += 28) {
                float y = -1 + j*2.f/(1024-1.f);
                float velocity_x = solution_pixels[4*(1024*j + i) + 0];
                float velocity_y = solution_pixels[4*(1024*j + i) + 1];
                const float thickness = 0.005;
                const vec4 color = vec4(0,0,0,1);
                const float epsilon = 1e-5;
                if (fabs(velocity_x) < epsilon && fabs(velocity_y) < epsilon) {
                } else {
                    world->graphics.paint.sphere(vec3(x, 0, y), 0.0075, vec4(0,0,0,1));
                    world->graphics.paint.line(vec3(x, 0, y), vec3(x + velocity_mul*velocity_x, 0, y + velocity_mul*velocity_y), thickness, color);
                }
            }
        }
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


    // Draw wireframe or base.
    if (wireframe) {
        world->graphics.paint.wireframe(*geom, mat4x4::translation(0,-0.01,0), 0.001);
        for (auto v : geom->mesh.vertices()) {
            world->graphics.paint.sphere(eigen_to_vec3(geom->position[v]), 0.0075, vec4(0.9,0.9,0.9,1));
        }
        for (auto e : geom->mesh.edges()) {
            world->graphics.paint.sphere(eigen_to_vec3(solver->midpoints[e]), 0.0075, vec4(0.9,0.9,0.9,1));
        }
    } else {
        // Draw boundaries.
        for (auto start : geom->mesh.boundary_loops()) {
            auto ps = std::vector<vec3>();
            auto he = start;
            do {
                ps.push_back(eigen_to_vec3(geom->position[he.vertex()]));
                he = he.next();
            } while (he != start);
	    ps.push_back(ps[0]);
            world->graphics.paint.chain(ps, 0.005, vec4(0,0,0,1));
        }
    }
    // Draw boundary condition.
    for (auto start : geom->mesh.boundary_loops()) {
        auto ps = std::vector<vec3>();
        auto he = start;
        do {
            auto v = he.vertex();
            auto e = he.edge();
            auto v_bv = vec3(solver->u_boundary[v].x(), 0, solver->u_boundary[v].y());
            auto e_bv = vec3(solver->u_boundary[e].x(), 0, solver->u_boundary[e].y());
            vec3 v_pos = eigen_to_vec3(geom->position[v]);
            vec3 e_pos = eigen_to_vec3(solver->midpoints[e]);
            vec3 shift = vec3(0,0.02,0);
            world->graphics.paint.line(v_pos+shift, shift+v_pos + v_bv*velocity_mul, 0.008, vec4(1,0.6,0.6,1));
            world->graphics.paint.line(e_pos+shift, shift+e_pos + e_bv*velocity_mul, 0.008, vec4(1,0.6,0.6,1));
            he = he.next();
        } while (he != start);
    }
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
