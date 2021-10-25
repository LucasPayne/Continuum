#include "NavierStokes/demo.h"
#include "NavierStokes/mesh_generators.h"


Demo::Demo()
{
    solver = nullptr;
    geom = nullptr;
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


    geom = square_mesh(10);
    double kinematic_viscosity = 1.;
    solver = new NavierStokesSolver(*geom, kinematic_viscosity);
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
    }
}

void Demo::update()
{
}


void Demo::post_render_update()
{
    double thickness = 0.005;
    world->graphics.paint.wireframe(*geom, mat4x4::translation(0,-0.01,0), thickness);
}


void Demo::mouse_handler(MouseEvent e)
{
}
