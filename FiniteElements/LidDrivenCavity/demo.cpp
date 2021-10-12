


struct Demo : public IBehaviour {
    Demo();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);

    void recreate_solver();

    SurfaceGeometry *geom;
    Solver *solver;

    int mesh_N;
    double mu;

    // Visualization.
    GLShaderProgram solution_shader; // Render the solution (velocity and pressure) to textures.
    GLuint solution_fbo; // For render-to-texture.
    GLuint solution_texture; // 3 components: r:velocity_x, g:velocity_y, b:pressure.
    GLuint solution_depth_texture; //--- Is this needed for completeness?
};


void Demo::recreate_solver()
{
    if (geom != nullptr) delete geom;
    // geom = circle_mesh(mesh_N, false);
    auto sq_mesh = SquareMesh(mesh_N);
    geom = sq_mesh.geom;

    if (solver != nullptr) delete solver;
    solver = new Solver(*geom, mu);

    // Set the lid boundary condition explicitly, on the vertex and midpoint sample points.
    // (This is to avoid possible errors at corners if the boundary condition was specified with a function.)
    for (int i = 0; i < mesh_N+1; i++) {
        solver->u_boundary[sq_mesh.vertex(i,0)] = vec2(1,0);
        // solver->u_boundary[sq_mesh.vertex(i,0)] = vec2(0,-1);
    }
    for (int i = 0; i < mesh_N; i++) {
        auto v1 = sq_mesh.vertex(i,0);
        auto v2 = sq_mesh.vertex(i+1,0);
        solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(1,0);
        // solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(0,-1);
    }
}


Demo::Demo()
{
    mesh_N = 4;
    mu = 1.0;
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
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_P) {
            mesh_N += 1;
            recreate_solver();
            solver->solve();
        }
        if (e.key.code == KEY_O) {
            mesh_N -= 1;
            if (mesh_N < 2) mesh_N = 2;
            recreate_solver();
            solver->solve();
        }
        if (e.key.code == KEY_R) {
            solver->solve();
        }
        if (e.key.code == KEY_T) {
            solver->write_sparsity_pattern = true;
            solver->solve();
        }
    }
}


void Demo::post_render_update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);

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


    // Render the solution textures.
    std::vector<vec2> position_data;
    std::vector<vec2> velocity_data;
    std::vector<float> pressure_data;
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
            pressure_data.push_back(solver->p[v]);
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
        {&position_data[0], sizeof(vec2), 2, GL_FLOAT}, // position
        {&velocity_data[0], sizeof(vec2), 2, GL_FLOAT}, // velocity
        {&pressure_data[0], sizeof(float), 1, GL_FLOAT} // pressure
    };
    // Upload the data.
    for (int i = 0; i < 3; i++) {
        auto metadata = data_to_upload[i];
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        glBufferData(GL_ARRAY_BUFFER, data_num_vertices * metadata.data_size, metadata.data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(i, metadata.gl_data_number, metadata.gl_data_type, GL_FALSE, 0, (const void *) 0);
        glEnableVertexAttribArray(i);
    }

    // Render to texture.
    solution_shader.bind();
    glBindFramebuffer(GL_FRAMEBUFFER, solution_fbo);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    glPatchParameteri(GL_PATCH_VERTICES, 6);
    glDrawArrays(GL_PATCHES, 0, data_num_vertices);
    glEnable(GL_DEPTH_TEST);

    glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);
    solution_shader.unbind();

    // Clean up.
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(3, vbos);

    


    // Draw the solution.
    double u_multiplier = 0.1;
    auto draw_u_vec = [&](Eigen::Vector3f pos, vec2 u_val) {
        // world->graphics.paint.sphere(eigen_to_vec3(pos), 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(eigen_to_vec3(pos), eigen_to_vec3(pos) + u_multiplier*vec3(u_val.x(), 0.05, u_val.y()), 0.005, vec4(0,0,0,1));
    };
    // Velocity
    for (auto v : geom->mesh.vertices()) {
        draw_u_vec(geom->position[v], solver->u[v]);
    }
    for (auto e : geom->mesh.edges()) {
        draw_u_vec(solver->midpoints[e], solver->u[e]);
    }
    // Pressure
    for (auto v : geom->mesh.vertices()) {
        double p = solver->p[v];
        vec4 color = vec4(p, p, p, 1);
        world->graphics.paint.sphere(eigen_to_vec3(geom->position[v]), 0.02, color);
    }

}
