#include "ReactionDiffusion/demo.h"

#include "mesh_generators.cpp"

Demo::Demo()
{
}

void Demo::take_screenshot()
{
    static int counter = 0;
    std::string pre = std::string(DATA) + "reaction_diffusion_" + std::to_string(counter);
    world->graphics.screenshot(pre + ".ppm",
                               screenshot_blx,
                               screenshot_bly,
                               screenshot_trx - screenshot_blx,
                               screenshot_try - screenshot_bly);
    counter += 1;
}

Aspect<Camera> main_camera;

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

    // Visualization
    solution_shader.add_shader(GLShader(VertexShader, SHADERS "reaction_diffusion/reaction_diffusion.vert"));
    solution_shader.add_shader(GLShader(FragmentShader, SHADERS "reaction_diffusion/reaction_diffusion.frag"));
    solution_shader.link();

    // geom = square_mesh(10);
    geom = torus_mesh(1.5, 0.4, 150);
    solver = new ReactorDiffuser(*geom, 0.02, [](vec3 pos, double u, double t)->double {
        return -10*u*(u-1)*(u+1);
    });
    solver->set_u([](vec3 pos)->double {
        return 0.;
    });

    show_wireframe = true;
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_P) solver->time_step(0.05);
        if (e.key.code == KEY_1) show_wireframe = !show_wireframe;
        if (e.key.code == KEY_R) {
            solver->set_u([](vec3 pos) {
                return 0.25 * sin(pos.x()*4) + 0.22*cos(pos.x()*pos.y()*10+1);
            });
        }
        // Take a screenshot.
        if (e.key.code == KEY_T) {
            take_screenshot();
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
}


void Demo::post_render_update()
{
    if (show_wireframe) {
        world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.002);
        return;
    }


    auto position_data = std::vector<vec3>();
    auto value_data = std::vector<float>();

    for (auto tri : geom->mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;

        // world->graphics.paint.sphere((1./3.)*(eigen_to_vec3(geom->position[he.vertex()])+eigen_to_vec3(geom->position[he.next().vertex()])+eigen_to_vec3(geom->position[he.next().next().vertex()])),
        //                             0.01, vec4(1,0,0,1));

        for (int i = 0; i < 3; i++) {
            Vertex v = he.vertex();
            position_data.push_back(eigen_to_vec3(geom->position[v]));
            value_data.push_back(solver->u_mesh[v]);
            he = he.next();
        }
    }
    assert(position_data.size() == value_data.size());

    GLuint vao;
    glCreateVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint vbos[2]; // position, value
    glGenBuffers(2, vbos);
    int data_num_vertices = position_data.size();
    struct {
        const void *data;
        size_t data_size;
        size_t gl_data_number;
        GLenum gl_data_type;
    } data_to_upload[2] = {
        {&position_data[0], sizeof(vec3), 3, GL_FLOAT},
        {&value_data[0], sizeof(float), 1, GL_FLOAT}
    };
    // Upload the data.
    for (int i = 0; i < 2; i++) {
        auto metadata = data_to_upload[i];
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        glBufferData(GL_ARRAY_BUFFER, data_num_vertices * metadata.data_size, metadata.data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(i, metadata.gl_data_number, metadata.gl_data_type, GL_FALSE, 0, (const void *) 0);
        glEnableVertexAttribArray(i);
    }

    glEnable(GL_DEPTH_TEST);
    solution_shader.bind();
    auto vp_matrix = main_camera->view_projection_matrix();
    glUniformMatrix4fv(solution_shader.uniform_location("mvp_matrix"), 1, GL_FALSE, (GLfloat *) &vp_matrix);
    glDrawArrays(GL_TRIANGLES, 0, data_num_vertices);
    solution_shader.unbind();


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
