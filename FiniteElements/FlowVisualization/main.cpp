#include "core.h"
#include "CameraController.h"
#include "P2_P1.h"
#include "mesh_generators.cpp"

Aspect<Camera> main_camera;


struct Solution {
    Solution(SurfaceGeometry &_geom) :
        geom{_geom},
        velocity{_geom.mesh},
        pressure{_geom.mesh}
    {
    }
    SurfaceGeometry &geom;
    P2Attachment<vec2> velocity;
    P1Attachment<double> pressure;
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

    Entity demo_e;
    Aspect<Behaviour> demo_b;

    char *directory_path;

    SurfaceGeometry &geom;

    // Load one frame of the solution.
    bool load_solution(int number);

    Solution sol;
    int solution_index;

    
    GLShaderProgram solution_shader; // Render the solution (velocity and pressure) to textures.
    GLuint solution_fbo; // For render-to-texture.
    GLuint solution_texture; // 3 components: r:velocity_x, g:velocity_y, b:pressure.
    GLuint solution_depth_texture; //--- Is this needed for completeness?
    GLShaderProgram render_solution_shader;
    GLuint render_solution_vao;

    void render_solution_texture();
    void render_solution(int mode); //velocity_x, velocity_y, pressure
};



App::App(World &_world, char *_directory_path, SurfaceGeometry &_geom) :
    world{_world},
    directory_path{_directory_path},
    geom{_geom},
    sol{geom}
{
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

    //---example solution frame
    // load_solution(std::string(directory_path) + "navier_stokes_1.txt");

    solution_index = 1;
    load_solution(1);


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
    render_solution(0);

    // world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.001);

    // for (auto v : geom.mesh.vertices()) {
    //     world.graphics.paint.sphere(eigen_to_vec3(geom.position[v]), 0.01, vec4(sol.velocity[v].length(), 0,0,1));
    // }

    std::cout << "Solution index: " << solution_index << "\n";
    
    if (world.input.keyboard.down(KEY_M)) {
        solution_index = (solution_index+1);
        load_solution(solution_index);
    }
    if (world.input.keyboard.down(KEY_N)) {
        solution_index = (solution_index-1);
        load_solution(solution_index);
    }
        
    // Draw boundaries.
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

    render_solution_texture();
}

void App::window_handler(WindowEvent e)
{
}

void App::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
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

    fclose(file);
    return true;
}


void App::render_solution_texture()
{
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
            scaled_pressure[v] = (sol.pressure[v] - min_pressure) / (max_pressure - min_pressure);
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
    GLuint vbos[5]; // position, velocity, pressure
    glGenBuffers(5, vbos);
    struct {
        const void *data;
        size_t data_size;
        size_t gl_data_number;
        GLenum gl_data_type;
    } data_to_upload[3] = {
        {&position_data[0], sizeof(vec2), 2, GL_FLOAT}, // position
        {&velocity_data[0], sizeof(vec2), 2, GL_FLOAT}, // velocity (P2)
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


int main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr, "give good args\n");
        exit(EXIT_FAILURE);
    }
    char *directory_path = argv[1];

    double theta0 = 0.13;
    vec2 obstruction_position = vec2(0.2,0.2);
    SurfaceGeometry *geom = square_minus_circle(0.25, theta0, 1, 1, 60, false, obstruction_position, false);

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
