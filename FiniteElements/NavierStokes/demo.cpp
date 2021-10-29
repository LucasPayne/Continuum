#include "NavierStokes/demo.h"
#include "NavierStokes/mesh_generators.h"

enum MeshModes {
    MM_square,
    MM_square_with_obstruction,
    NUM_MESH_MODES
};

Demo::Demo()
{
    solver = nullptr;
    geom = nullptr;
}

void Demo::recreate_solver()
{
    if (solver != nullptr) delete solver;
    if (geom != nullptr) delete geom;
    double kinematic_viscosity = 0.01;

    if (mesh_mode == MM_square) {
        geom = square_mesh(12);
        solver = new NavierStokesSolver(*geom, kinematic_viscosity);
        solver->set_source(
            [&](double x, double y, double t)->vec2 {
                // const double r = 0.175;
                // if (x*x + y*y <= r*r) return vec2(25,0);
                const double r = 0.125;
                // if ((x-0.5)*(x-0.5) + y*y <= r*r) return vec2(-25, 0);

                // if ((x-0.5)*(x-0.5) + y*y <= r*r) return vec2(0, 50);
                // if ((x+0.5)*(x+0.5) + y*y <= r*r) return vec2(25, 0);
                if (t < 0.1) {
                    if ((x-0.5)*(x-0.5) + y*y <= r*r) return vec2(0, 50);
                }

                return vec2(0,0);
            }
        );
    } else if (mesh_mode == MM_square_with_obstruction) {
        double theta0 = 0.13;
        vec2 obstruction_position = vec2(0,0);
        geom = square_minus_circle(0.25, theta0, 1, 1, 8, false, obstruction_position, false);
        solver = new NavierStokesSolver(*geom, kinematic_viscosity);
        solver->set_source(
            [&](double x, double y, double t)->vec2 {
                // const double r = 0.175;
                // if (x*x + y*y <= r*r) return vec2(25,0);
                const double r = 0.125;
                if ((x-0.5)*(x-0.5) + y*y <= r*r) return vec2(0, 50);
                if ((x+0.5)*(x+0.5) + y*y <= r*r) return vec2(25, 0);
                return vec2(0,0);
            }
        );
    }
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

    
    /*--------------------------------------------------------------------------------
        Mesh selection
    --------------------------------------------------------------------------------*/
    mesh_mode = 0;
    // Create the solver.
    recreate_solver();
    
    /*--------------------------------------------------------------------------------
    Visualization
    --------------------------------------------------------------------------------*/
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

    show_wireframe = true;
    show_vector_field = true;
    show_source = true;
    filming = false;
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        // Separate Newton iterations for debugging.
        if (e.key.code == KEY_I) solver->start_time_step(0.01);
        if (e.key.code == KEY_O) solver->newton_iteration();
        if (e.key.code == KEY_P) solver->end_time_step();

        // One iteration.
        if (e.key.code == KEY_R) solver->time_step(0.0025);

        // Take a screenshot.
        if (e.key.code == KEY_T) {
            take_screenshot();
        }
        // Take a sequence of screenshots, to be converted to a video.
        if (e.key.code == KEY_9) {
            film_seconds = 15;
            // film_dt = 1./25.;
            film_dt = 1./100.;
            film_num_frames = ceil(film_seconds / film_dt);
            film_frame = 0;
            filming = true;
        }

        // Select mesh.
        if (e.key.code == KEY_M) {
            mesh_mode = (mesh_mode + 1) % NUM_MESH_MODES;
            recreate_solver();
        }
        if (e.key.code == KEY_N) {
            mesh_mode -= 1;
            if (mesh_mode < 0) mesh_mode = NUM_MESH_MODES-1;
            recreate_solver();
        }

        // Solver toggles
        if (e.key.code == KEY_M) {
            solver->m_use_advection = !solver->m_use_advection;
        }
        
        // Rendering toggles.
        if (e.key.code == KEY_1) {
            show_wireframe = !show_wireframe;
        }
        if (e.key.code == KEY_2) {
            show_vector_field = !show_vector_field;
        }
        if (e.key.code == KEY_3) {
            show_source = !show_source;
        }
        if (e.key.code == KEY_4) {
            solver->set_velocity([&](double x, double y)->vec2 {
                // return 0.06*vec2(exp(x), exp(3*y+x));
                const double r = 0.125;
                if ((vec2(x,y) - source_position).length() <= r) {
                    return vec2(2, 0);
                }
                return vec2(0,0);
            });
        }
    }
}

void Demo::update()
{
    double source_move_speed = 1.;
    if (world->input.keyboard.down(KEY_LEFT_ARROW)) {
        source_position.x() += source_move_speed * dt;
    }
    if (world->input.keyboard.down(KEY_RIGHT_ARROW)) {
        source_position.x() -= source_move_speed * dt;
    }
    if (world->input.keyboard.down(KEY_DOWN_ARROW)) {
        source_position.y() -= source_move_speed * dt;
    }
    if (world->input.keyboard.down(KEY_UP_ARROW)) {
        source_position.y() += source_move_speed * dt;
    }
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

    // Filming
    if (filming) {
        take_screenshot();
        film_frame += 1;
        if (film_frame == film_num_frames) {
            filming = false;
            exit(EXIT_SUCCESS);
        }
    }

    solver->set_source(
        [&](double x, double y, double t)->vec2 {
            const double r = 0.125;
            if ((vec2(x,y) - source_position).length() <= r) {
                return vec2(0, 30);
            }
            return vec2(0,0);
        }
    );
}


void Demo::post_render_update()
{
    // Filming a video, update simulation
    if (filming) {
        solver->time_step(film_dt);
    }

    if (show_wireframe) {
        double thickness = 0.005;
        world->graphics.paint.wireframe(*geom, mat4x4::translation(0,-0.01,0), thickness);
    }
    // Draw boundaries.
    for (auto start : geom->mesh.boundary_loops()) {
        auto ps = std::vector<vec3>();
        auto he = start;
        do {
            ps.push_back(vec3(0,0.0001,0)+eigen_to_vec3(geom->position[he.vertex()]));
            he = he.next();
        } while (he != start);
        ps.push_back(ps[0]);
        world->graphics.paint.chain(ps, 0.005, vec4(0,0,0,1));
    }
    
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
            scaled_pressure[v] = (solver->pressure[v] - min_pressure) / (max_pressure - min_pressure);
        }
    } else {
        for (auto v : geom->mesh.vertices()) {
            scaled_pressure[v] = 0;
        }
    }

    // Render the solution textures.
    auto position_data = std::vector<vec2>();
    auto velocity_data = std::vector<vec2>();
    auto pressure_data = std::vector<float>();
    auto div_u_data = std::vector<float>();
    auto div_u_P2_data = std::vector<float>();
    // Create a 6-vertex patch per triangle.
    for (auto tri : geom->mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        do {
            auto v = he.vertex();
            auto e = he.edge();
            Eigen::Vector3f midpoint = 0.5*geom->position[e.a().vertex()] + 0.5*geom->position[e.b().vertex()];
            for (auto pos : {geom->position[v], midpoint}) {
                position_data.push_back(vec2(pos.x(), pos.z()));
            }
            velocity_data.push_back(solver->velocity[v]);
            velocity_data.push_back(solver->velocity[e]);
            pressure_data.push_back(scaled_pressure[v]);
            pressure_data.push_back(0.f); // Dummy data, as there is no midpoint pressure.
            // --- filler. Don't have this data.
            div_u_data.push_back(0.f);
            div_u_data.push_back(0.f);
            div_u_P2_data.push_back(0.f);
            div_u_P2_data.push_back(0.f);

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
    // glUniform1i(solution_shader.uniform_location("show_div_P2"), show_div_P2 ? 1 : 0);
    glUniform1i(solution_shader.uniform_location("show_div_P2"), 0); //---no divergence currently
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
    // Also draw source velocity field.
    if (show_vector_field) {
        glBindFramebuffer(GL_FRAMEBUFFER, solution_fbo);
        auto solution_pixels = std::vector<float>(4*1024*1024);
        glReadPixels(0,0,1024,1024, GL_RGBA, GL_FLOAT, &solution_pixels[0]);
        glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);
        std::vector<float> pressures;
        // const int skip = 20;
        // const int skip = 14;
        const int skip = 23;
	std::vector<vec2> positions_2D;
        for (int i = 0; i < 1024; i += skip) {
            float x = -1 + i*2.f/(1024-1.f);
            for (int j = 0; j < 1024; j += skip) {
                float y = -1 + j*2.f/(1024-1.f);
                float velocity_x = solution_pixels[4*(1024*j + i) + 0];
                float velocity_y = solution_pixels[4*(1024*j + i) + 1];
                float pressure = solution_pixels[4*(1024*j + i) + 2];
                const float thickness = 0.002;
                const vec4 color = vec4(0,0,0,1);
                const float epsilon = 1e-5;

                // Draw velocity vector.
                const float velocity_mul = 0.12;
                if (fabs(velocity_x) >= epsilon || fabs(velocity_y) >= epsilon) {
                    world->graphics.paint.line(vec3(x, 0.005, y), vec3(x + velocity_mul*velocity_x, 0, y + velocity_mul*velocity_y), thickness, color);
	            pressures.push_back(pressure);
                    positions_2D.push_back(vec2(x,y));
                }
                // Draw source vector.
                if (show_source) {
                    const float source_mul = velocity_mul / 25.f;
                    vec2 source_velocity = solver->source_function(x,y, solver->time());
                    const vec4 source_color = vec4(0,0,1,1);
                    if (fabs(source_velocity.x()) >= epsilon || fabs(source_velocity.y()) >= epsilon) {
                        world->graphics.paint.line(vec3(x, 0.005, y),
                                                   vec3(x + source_mul*source_velocity.x(), 0, y + source_mul*source_velocity.y()),
                                                    thickness*0.6, source_color);
                    }
                }
            }
        }
    } //endif show_vector_field


}


void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
        // Set the screenshot rectangle.
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

void Demo::take_screenshot()
{
    static int counter = 0;
    std::string pre = std::string(DATA) + "navier_stokes_" + std::to_string(counter);
    world->graphics.screenshot(pre + ".ppm",
                               screenshot_blx,
                               screenshot_bly,
                               screenshot_trx - screenshot_blx,
                               screenshot_try - screenshot_bly);
    counter += 1;
}
