#include "ReactionDiffusion/demo.h"

#include "mesh_generators.cpp"

Demo::Demo()
{
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

    // geom = square_mesh(10);
    geom = torus_mesh(1.5, 0.4, 50);
    solver = new ReactorDiffuser(*geom, 1, [](vec3 pos, double u, double t)->double {
        return 1;
    });
    solver->set_u([](vec3 pos) {
        return 0.;
    });
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_P) solver->time_step(0.05);
    }
}

void Demo::update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.01);
}


void Demo::post_render_update()
{
}


void Demo::mouse_handler(MouseEvent e)
{
}
