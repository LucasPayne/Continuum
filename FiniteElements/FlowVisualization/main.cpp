#include "core.h"
#include "CameraController.h"
#include "P2_P1.h"
#include "mesh_generators.cpp"
#include "Filmer.h"

Aspect<Camera> main_camera;

bool HOLD_MODE = false;


struct Solution {
    Solution(SurfaceGeometry &_geom) :
        geom{_geom},
        velocity{_geom.mesh},
        pressure{_geom.mesh},
        divergence{_geom.mesh},
        divergence_linear{_geom.mesh}
    {
    }
    SurfaceGeometry &geom;
    P2Attachment<vec2> velocity;
    P1Attachment<double> pressure;
    P2Attachment<double> divergence;
    P1Attachment<double> divergence_linear;
};


class App : public IGC::Callbacks {
public:
    World &world;
    App(World &world, char *_directory_path, SurfaceGeometry &_geom);

    void close();
    void loop();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void window_handler(WindowEvent e);

    FILE *divergence_residual_file;

    char *directory_path;

    SurfaceGeometry &geom;

    // Load one frame of the solution.
    bool load_solution(int number);

    Solution sol;
    int solution_index;

    Filmer *filmer;
    
    GLShaderProgram solution_shader; // Render the solution (velocity and pressure) to textures.
    GLuint solution_fbo; // For render-to-texture.
    GLuint solution_texture; // 3 components: r:velocity_x, g:velocity_y, b:pressure.
    GLuint solution_depth_texture; //--- Is this needed for completeness?
    GLShaderProgram render_solution_shader;
    GLuint render_solution_vao;

    int render_solution_mode;

    void render_solution_texture();
    void render_solution(int mode); //velocity_x, velocity_y, pressure
    
    // Divergence computation.
    SparseMatrix gramian_matrix_P2_0();
    void div_P2_P2(P2Attachment<vec2> &vf, P2Attachment<double> &div);

    SparseMatrix gramian_matrix_P1_0();
    void div_P2_P1(P2Attachment<vec2> &vf, P1Attachment<double> &div);
};



App::App(World &_world, char *_directory_path, SurfaceGeometry &_geom) :
    world{_world},
    directory_path{_directory_path},
    geom{_geom},
    sol{geom}
{
    std::string fname = std::string(DATA) + "divergence_residuals.txt";
    divergence_residual_file = fopen(fname.c_str(), "w+");

    // OpenGL ----------------------------------------------------
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //------------------------------------------------------------
    //
    printf("Loading from path: %s\n", directory_path);

    // Create a camera controller.
    auto cameraman = world.entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    camera->background_color = vec4(1,1,1,1);
    auto t = cameraman.add<Transform>(0,2,0);
    main_camera = camera;
    CameraController *controller = world.add<CameraController>(cameraman);
    controller->angle = -M_PI/2;
    controller->azimuth = M_PI;

    // Create the screenshotter.
    auto filmer_e = world.entities.add();
    filmer = world.add<Filmer>(filmer_e, std::string(DATA) + "navier_render", [&](int film_frame) {
        solution_index = film_frame;
        load_solution(solution_index);
        render_solution(render_solution_mode);
    }, 1);

    //---example solution frame
    // load_solution(std::string(directory_path) + "navier_stokes_1.txt");

    solution_index = 0;
    load_solution(solution_index);

    solution_shader.add_shader(GLShader(VertexShader, SHADERS "flow_visualization_solution/solution.vert"));
    solution_shader.add_shader(GLShader(TessControlShader, SHADERS "flow_visualization_solution/solution.tcs"));
    solution_shader.add_shader(GLShader(TessEvaluationShader, SHADERS "flow_visualization_solution/solution.tes"));
    solution_shader.add_shader(GLShader(FragmentShader, SHADERS "flow_visualization_solution/solution.frag"));
    solution_shader.link();

    render_solution_shader.add_shader(GLShader(VertexShader, SHADERS "flow_visualization_render_solution/render_solution.vert"));
    render_solution_shader.add_shader(GLShader(FragmentShader, SHADERS "flow_visualization_render_solution/render_solution.frag"));
    render_solution_shader.link();
    
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

void App::close()
{
}

void App::loop()
{
    if (!filmer->filming) {
        render_solution(render_solution_mode);
    }

    // world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.001);

    // for (auto v : geom.mesh.vertices()) {
    //     world.graphics.paint.sphere(eigen_to_vec3(geom.position[v]), 0.01, vec4(sol.velocity[v].length(), 0,0,1));
    // }

    std::cout << "Solution index: " << solution_index << "\n";
    
    if (HOLD_MODE) {
        if (world.input.keyboard.down(KEY_M)) {
            solution_index = (solution_index+1);
            load_solution(solution_index);
        }
        if (world.input.keyboard.down(KEY_N)) {
            solution_index = (solution_index-1);
            load_solution(solution_index);
        }
    }
        
    // Draw boundaries.
    #if 0
    for (auto start : geom.mesh.boundary_loops()) {
        auto ps = std::vector<vec3>();
        auto he = start;
        do {
            ps.push_back(vec3(0,0,0)+eigen_to_vec3(geom.position[he.vertex()]));
            he = he.next();
        } while (he != start);
        ps.push_back(ps[0]);
        world.graphics.paint.chain(ps, 0.005, vec4(0,0,0,1));
    }
    #endif

    render_solution_texture();
}

void App::window_handler(WindowEvent e)
{
}

void App::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) {
            fclose(divergence_residual_file);
            exit(EXIT_SUCCESS);
        }
        if (e.key.code == KEY_1) render_solution_mode = 0;
        if (e.key.code == KEY_2) render_solution_mode = 1;
        if (e.key.code == KEY_3) render_solution_mode = 2;
        if (e.key.code == KEY_4) render_solution_mode = 3;
        if (!HOLD_MODE) {
            if (e.key.code == KEY_M) {
                solution_index = (solution_index+1);
                load_solution(solution_index);
            }
            if (e.key.code == KEY_N) {
                solution_index = (solution_index-1);
                load_solution(solution_index);
            }
        }
    }
}

void App::mouse_handler(MouseEvent e)
{
}



bool App::load_solution(int number)
{
    std::string path = std::string(directory_path) + "navier_stokes_" + std::to_string(number) + ".txt";
    printf("Loading \"%s\"\n", path.c_str());
    FILE *file = fopen(path.c_str(), "r");
    assert(file != nullptr);
    // Load velocity
    for (auto v : geom.mesh.vertices()) {
        float x, y;
        assert(fscanf(file, "%f %f\n", &x, &y) == 2);
        sol.velocity[v] = vec2(x,y);
    }
    for (auto e : geom.mesh.edges()) {
        float x, y;
        assert(fscanf(file, "%f %f\n", &x, &y) == 2);
        sol.velocity[e] = vec2(x,y);
    }
    // Load pressure
    for (auto v : geom.mesh.vertices()) {
        float p;
        assert(fscanf(file, "%f\n", &p) == 1);
        sol.pressure[v] = p;
    }

    // Compute the divergence
    // div_P2_P2(sol.velocity, sol.divergence);
    // div_P2_P1(sol.velocity, sol.divergence_linear);

    fclose(file);
    return true;
}


void App::render_solution_texture()
{
    const double pressure_mul = 1.;

    // Scale the pressure.
    VertexAttachment<double> scaled_pressure(geom.mesh);
    double min_pressure = std::numeric_limits<double>::infinity();
    double max_pressure = 0;
    for (auto v : geom.mesh.vertices()) {
        if (sol.pressure[v] < min_pressure) {
            min_pressure = sol.pressure[v];
        }
        if (sol.pressure[v] > max_pressure) {
            max_pressure = sol.pressure[v];
        }
    }
    if (max_pressure != min_pressure) {
        for (auto v : geom.mesh.vertices()) {
            scaled_pressure[v] = pressure_mul*(sol.pressure[v] - min_pressure) / (max_pressure - min_pressure);
        }
    } else {
        for (auto v : geom.mesh.vertices()) {
            scaled_pressure[v] = 0;
        }
    }

    // Render the solution textures.
    auto position_data = std::vector<vec2>();
    auto velocity_data = std::vector<vec2>();
    auto pressure_data = std::vector<float>();
    auto divergence_data = std::vector<float>();
    auto divergence_linear_data = std::vector<float>();
    // Create a 6-vertex patch per triangle.
    for (auto tri : geom.mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        do {
            auto v = he.vertex();
            auto e = he.edge();
            Eigen::Vector3f midpoint = 0.5*geom.position[e.a().vertex()] + 0.5*geom.position[e.b().vertex()];
            for (auto pos : {geom.position[v], midpoint}) {
                position_data.push_back(vec2(pos.x(), pos.z()));
            }
            velocity_data.push_back(sol.velocity[v]);
            velocity_data.push_back(sol.velocity[e]);
            pressure_data.push_back(scaled_pressure[v]);
            pressure_data.push_back(0.f); // Dummy data, as there is no midpoint pressure.
            divergence_data.push_back(sol.divergence[v]);
            divergence_data.push_back(sol.divergence[e]);
            divergence_linear_data.push_back(sol.divergence_linear[v]);
            divergence_linear_data.push_back(0.f); // Dummy data.

            he = he.next();
        } while (he != start);
    }
    // Check that the right number of vertices were created.
    size_t data_num_vertices = 6*geom.mesh.num_faces();
    for (size_t len : {position_data.size(), velocity_data.size(), pressure_data.size()}) {
        assert(len == data_num_vertices);
    }

    GLuint vao;
    glCreateVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint vbos[5]; // position, velocity, pressure, divergence, divergence_linear
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
        {&divergence_data[0], sizeof(float), 1, GL_FLOAT}, // divergence (P2)
        {&divergence_linear_data[0], sizeof(float), 1, GL_FLOAT}, // divergence_linear (P1)
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

    glBindFramebuffer(GL_FRAMEBUFFER, world.graphics.screen_buffer.id);
    glEnable(GL_SCISSOR_TEST);
    glEnable(GL_DEPTH_TEST);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    solution_shader.unbind();

    // Clean up.
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(3, vbos);
}


void App::render_solution(int mode)
{
    GLuint vao;
    glCreateVertexArrays(1, &vao);
    glBindVertexArray(vao);

    vec3 positions[4] = {
        vec3(-1,0,-1),
        vec3(1,0,-1),
        vec3(1,0,1),
        vec3(-1,0,1)
    };
    vec2 uvs[4] = {
        vec2(0,0),
        vec2(1,0),
        vec2(1,1),
        vec2(0,1)
    };
    GLuint vbos[2];
    glGenBuffers(2, vbos);
    struct {
        const void *data;
        size_t data_size;
        size_t gl_data_number;
        GLenum gl_data_type;
    } data_to_upload[2] = {
        {&positions[0], sizeof(vec3), 3, GL_FLOAT},
        {&uvs[0], sizeof(vec2), 2, GL_FLOAT},
    };
    // Upload the data.
    for (int i = 0; i < 2; i++) {
        auto metadata = data_to_upload[i];
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        glBufferData(GL_ARRAY_BUFFER, 4 * metadata.data_size, metadata.data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(i, metadata.gl_data_number, metadata.gl_data_type, GL_FALSE, 0, (const void *) 0);
        glEnableVertexAttribArray(i);
    }

    render_solution_shader.bind();
    glUniform1i(render_solution_shader.uniform_location("mode"), mode);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, solution_texture);
    glUniform1i(render_solution_shader.uniform_location("solution"), 0);
    auto vp_matrix = main_camera->view_projection_matrix();
    glUniformMatrix4fv(render_solution_shader.uniform_location("mvp_matrix"), 1, GL_FALSE, (const GLfloat *) &vp_matrix);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    render_solution_shader.unbind();

    // Clean up
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(2, vbos);
}


// Gramian projection matrix for P2_0 (zero on the boundary).
SparseMatrix App::gramian_matrix_P2_0()
{
    VertexAttachment<int> interior_vertex_indices(geom.mesh);
    EdgeAttachment<int> interior_midpoint_indices(geom.mesh);
    int index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            interior_vertex_indices[v] = -1;
            continue;
        }
        interior_vertex_indices[v] = index;
        index += 1;
    }
    index = 0;
    for (auto e : geom.mesh.edges()) {
        if (e.on_boundary()) {
            interior_midpoint_indices[e] = -1;
            continue;
        }
        interior_midpoint_indices[e] = index;
        index += 1;
    }

    std::vector<EigenTriplet> coefficients;
    int N = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    auto add_entry = [&](int i, int j, double value) {
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^us on a vertex.
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto start = v.halfedge();
        auto he = start;
        int v_index = counter;
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            int edge_110_index = interior_midpoint_indices[edge_110];
            int edge_011_index = interior_midpoint_indices[edge_011];
            int edge_101_index = interior_midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

            add_entry(v_index, v_index, (1./60.)*R);
            if (!vp.on_boundary()) add_entry(v_index, vp_index, (-1./360.)*R);
            if (!vpp.on_boundary()) add_entry(v_index, vpp_index, (-1./360.)*R);

            if (!edge_110.on_boundary()) add_entry(v_index, geom.mesh.num_interior_vertices() + edge_110_index, (-1./90.)*R);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
        counter += 1;
    }

    counter = 0;
    // For each psi^us at a midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = counter;
        Halfedge hes[2] = {edge.a(), edge.b()};

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // if (tri.null()) continue;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            int v_index = interior_vertex_indices[v];
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            int edge_110_index = interior_midpoint_indices[edge_110];
            int edge_011_index = interior_midpoint_indices[edge_011];
            int edge_101_index = interior_midpoint_indices[edge_101];

            double R = 2*geom.triangle_area(tri);

	    if (!edge_110.on_boundary()) add_entry(geom.mesh.num_interior_vertices() + edge_110_index,
		      geom.mesh.num_interior_vertices() + edge_110_index, (4./45.)*R);
            if (!edge_011.on_boundary()) add_entry(geom.mesh.num_interior_vertices() + edge_110_index,
                      geom.mesh.num_interior_vertices() + edge_011_index, (2./45.)*R);
            
            if (!edge_101.on_boundary()) add_entry(geom.mesh.num_interior_vertices() + edge_110_index,
                      geom.mesh.num_interior_vertices() + edge_101_index, (2./45.)*R);
            
            if (!v.on_boundary()) add_entry(geom.mesh.num_interior_vertices() + edge_110_index,
                      v_index, (-1./90.)*R);
            
        }
        counter += 1;
    }
    
    auto matrix = SparseMatrix(N, N);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    // make_sparsity_image(matrix, DATA "P20_P20_gramian.ppm");

    return matrix;
}


// Compute the divergence of vf in P2_2 projected into P2.
void App::div_P2_P2(P2Attachment<vec2> &vf, P2Attachment<double> &div)
{
    int N = geom.mesh.num_interior_vertices() + geom.mesh.num_interior_edges();
    auto div_vector_proj = Eigen::VectorXd(N);
    for (int i = 0; i < N; i++) div_vector_proj[i] = 0.;
    
    int counter = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = counter;
        auto v_pos = geom.position[v];

        double integral = 0.;

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            vec2 vf_002 = vf[v];
            vec2 vf_200 = vf[vp];
            vec2 vf_020 = vf[vpp];
            vec2 vf_101 = vf[edge_101];
            vec2 vf_011 = vf[edge_011];
            vec2 vf_110 = vf[edge_110];

            // integral += (1./15.)*vec2::dot(K3.perp(), vf_002);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_020);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_200);

            // integral += (1./10.)*vec2::dot(K3.perp(), vf_011);
            // integral += (-1./30.)*vec2::dot(K3.perp(), vf_110);
            // integral += (1./10.)*vec2::dot(K3.perp(), vf_101);
            
            integral += (1./15.)*vec2::dot(K3.perp(), vf_002);
            integral += (-1./30.)*vec2::dot(K2.perp(), vf_020);
            integral += (-1./30.)*vec2::dot(K1.perp(), vf_200);

            integral += vec2::dot((1./30.)*K1.perp() + (1./10.)*K2.perp(), vf_011);
            integral += vec2::dot((-1./30.)*K1.perp() + (-1./30.)*K2.perp(), vf_110);
            integral += vec2::dot((1./10.)*K1.perp() + (1./30.)*K2.perp(), vf_101);
            

            he = he.twin().next();
        } while (!he.face().null() && he != start);

        div_vector_proj[v_index] -= integral;
        counter += 1;
    }

    counter = 0;
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = counter;
        Halfedge hes[2] = {edge.a(), edge.b()};

        double integral = 0.;

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            if (tri.null()) continue;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            auto edge_110 = he.edge(); // vp to vpp
            auto edge_011 = he.next().edge(); // vpp to v
            auto edge_101 = he.next().next().edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            vec2 vf_002 = vf[v];
            vec2 vf_200 = vf[vp];
            vec2 vf_020 = vf[vpp];
            vec2 vf_101 = vf[edge_101];
            vec2 vf_011 = vf[edge_011];
            vec2 vf_110 = vf[edge_110];
            // integral += (1./30.)*vec2::dot(K3.perp(), vf_002);
            // integral += vec2::dot((1./15.)*K1.perp() + (-1./30.)*K2.perp(), vf_020);
            // integral += vec2::dot((-1./30.)*K1.perp() + (1./15.)*K2.perp(), vf_200);

            // integral += vec2::dot((4./15.)*K1.perp() + (2./15.)*K2.perp(), vf_011);
            // integral += vec2::dot((4./15.)*K1.perp() + (4./15.)*K2.perp(), vf_110);
            // integral += vec2::dot((2./15.)*K1.perp() + (4./15.)*K2.perp(), vf_101);
            integral += vec2::dot((-1./30.)*K3.perp(), vf_002);
            integral += vec2::dot((1./10.)*K2.perp(), vf_020);
            integral += vec2::dot((1./10.)*K1.perp(), vf_200);

            integral += vec2::dot((-4./15.)*K1.perp() + (-2./15.)*K2.perp(), vf_011);
            integral += vec2::dot((4./15.)*K1.perp() + (4./15.)*K2.perp(), vf_110);
            integral += vec2::dot((-2./15.)*K1.perp() + (-4./15.)*K2.perp(), vf_101);
        }
        div_vector_proj[geom.mesh.num_interior_vertices() + midpoint_index] -= integral;
        counter += 1;
    }
    
    SparseMatrix gramian = gramian_matrix_P2_0();
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(gramian);
    solver.factorize(gramian);
    Eigen::VectorXd div_vector = solver.solve(div_vector_proj);

    // std::cout << "div_vector\n";
    // std::cout << div_vector << "\n";
    // std::cout << "div_vector_proj\n";
    // std::cout << div_vector_proj << "\n";
    // getchar();

    // Associate this divergence to the P2 mesh.
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (!v.on_boundary()) {
            div[v] = div_vector[interior_vertex_index];
            interior_vertex_index += 1;
        } else {
            div[v] = 0.;
        }
    }
    int interior_edge_index = 0;
    for (auto e : geom.mesh.edges()) {
        if (!e.on_boundary()) {
            div[e] = div_vector[geom.mesh.num_interior_vertices() + interior_edge_index];
            interior_edge_index += 1;
        } else {
            div[e] = 0.;
        }
    }
}

SparseMatrix App::gramian_matrix_P1_0()
{
    std::vector<EigenTriplet> coefficients;
    int N_p = geom.mesh.num_vertices();

    VertexAttachment<int> vertex_indices(geom.mesh);
    int index = 0;
    for (auto v : geom.mesh.vertices()) {
        vertex_indices[v] = index;
        index += 1;
    }

    auto add_entry = [&](int i, int j, double value) {
        // printf("%d %d %.6f\n", i, j, value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    // For each psi^p.
    for (auto v : geom.mesh.vertices()) {
        auto start = v.halfedge();
        auto he = start;
        int v_index = vertex_indices[v];
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];

            double R = 2*geom.triangle_area(tri);

            add_entry(v_index, v_index, (1./12.)*R);
            add_entry(v_index, vp_index, (1./24.)*R);
            add_entry(v_index, vpp_index, (1./24.)*R);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);
    }

    auto matrix = SparseMatrix(N_p, N_p);
    matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    matrix.makeCompressed();
    // make_sparsity_image(matrix, DATA "pressure_gramian.ppm");

    return matrix;
}

void App::div_P2_P1(P2Attachment<vec2> &vf, P1Attachment<double> &div)
{
    int N_p = geom.mesh.num_vertices();
    auto u_div_l2_proj = Eigen::VectorXd(N_p);
    for (int i = 0; i < N_p; i++) u_div_l2_proj[i] = 0.;


    VertexAttachment<int> vertex_indices(geom.mesh);
    int index = 0;
    for (auto v : geom.mesh.vertices()) {
        vertex_indices[v] = index;
        index += 1;
    }

    // Compute the divergence term.
    // For each psi^p.
    for (auto v : geom.mesh.vertices()) {
        int v_index = vertex_indices[v];
        auto v_pos = geom.position[v];

        double integral = 0.;

        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            auto edge_110 = he.next().edge(); // vp to vpp
            auto edge_011 = he.next().next().edge(); // vpp to v
            auto edge_101 = he.edge(); // v to vp
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            vec2 u110 = vf[edge_110];
            vec2 u011 = vf[edge_011];
            vec2 u101 = vf[edge_101];
            vec2 u002 = vf[v];
            integral += (-1./6.) * vec2::dot(K3.perp(), u002);
            integral += (1./6.) * vec2::dot(K3.perp(), u110);
            integral += (-1./6.) * vec2::dot((K2 - K1).perp(), u011);
            integral += (-1./6.) * vec2::dot((K1 - K2).perp(), u101);
            
            he = he.twin().next();
        } while (!he.face().null() && he != start);

        u_div_l2_proj[v_index] = integral;
    }

    // Reconstruct the P1 divergence of u, then associate it to the mesh.
    SparseMatrix gramian = gramian_matrix_P1_0();
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(gramian);
    solver.factorize(gramian);
    Eigen::VectorXd div_u_vector = solver.solve(u_div_l2_proj);
    index = 0;
    float abs_max_div = 0;
    for (auto v : geom.mesh.vertices()) {
        div[v] = div_u_vector[index];
        if (fabs(div[v]) > abs_max_div) abs_max_div = fabs(div[v]);
        index += 1;
    }
    // Write the divergence residual to a text file.
    fprintf(divergence_residual_file, "%.12f\n", abs_max_div);
    
}

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4) {
        fprintf(stderr, "give good args\n");
        exit(EXIT_FAILURE);
    }
    char *directory_path = argv[1];
    int GEOM = 0;
    sscanf(argv[2], "%d", &GEOM);
    int square_n = 0;
    if (argc == 4) {
        sscanf(argv[3], "%d", &square_n);
    }

    // navier_3
    SurfaceGeometry *geom = nullptr;
    if (GEOM == 3) {
        double theta0 = 0.13;
        vec2 obstruction_position = vec2(0.2,0.2);
        geom = square_minus_circle(0.25, theta0, 1, 1, 60, false, obstruction_position, false);
    } else if (GEOM == 4) {
        assert(square_n > 0);
        double theta0 = 0.1257;
        vec2 obstruction_position = vec2(0,0);
        geom = square_minus_circle(0.18, theta0, 1, 1, square_n, true, obstruction_position, false);
    } else if (GEOM == 5) {
        assert(square_n > 0);
        geom = square_mesh(square_n);
    } else if (GEOM == 6) {
        double theta0 = 0.13;
        vec2 obstruction_position = vec2(0.2,0.2);
        geom = square_minus_circle(0.25, theta0, 1, 1, 60, false, obstruction_position, false);
    }

    assert(geom != nullptr);

    printf("[main] Creating context...\n");
    IGC::Context context("A world");
    printf("[main] Creating world...\n");
    World world(context);

    printf("[main] Creating app...\n");
    App app(world, directory_path, *geom);
    printf("[main] Adding app callbacks...\n");
    context.add_callbacks(&app);

    printf("[main] Entering loop...\n");
    context.enter_loop();
    context.close();
}
