#include "core.h"
#include "SurfaceNavierStokes/demo.h"
#include "mesh_generators.cpp"
#include "mesh_processing/extensions/assimp_convert.h"

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

    geom = assimp_to_surface_geometry(std::string(MODELS) + "cylinder.stl");
    geom->mesh.lock();
    solver = new SurfaceNavierStokesSolver(*geom, 1);

}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_R) {
            solver->time_step(0.1);
        }
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
    if (e.action == MOUSE_BUTTON_PRESS) {
    }
}
