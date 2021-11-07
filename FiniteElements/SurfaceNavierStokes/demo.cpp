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

    // geom = assimp_to_surface_geometry(std::string(MODELS) + "tangram.stl");
    // geom = assimp_to_surface_geometry(std::string(MODELS) + "cylinder.stl");
    // geom->mesh.lock();
    
    geom = square_mesh(10);
    for (auto v : geom->mesh.vertices()) {
        vec3 p = eigen_to_vec3(geom->position[v]);
        geom->position[v] += Eigen::Vector3f(0, 0.1*sin(4*p.x()), 0);
        // geom->position[v] = Eigen::Vector3f(p.y(), p.x(), p.z());
    }
    
    solver = new SurfaceNavierStokesSolver(*geom, 1);

    solver->set_source([&](double x, double y, double z)->vec3 {
        double r = 0.25;
        if (x*x + z*z <= r*r) return vec3(1,0,0);
        // if (y*y + z*z <= r*r) return vec3(1,0,1);
        return vec3(0,0,0);
    });
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_R) {
            solver->time_step(0.1);
        }
        if (e.key.code == KEY_1) {
        }
    }
}

void Demo::update()
{
    // world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
    double velocity_mul = 1;
    for (auto v : geom->mesh.vertices()) {
        vec3 p = eigen_to_vec3(geom->position[v]);
        vec3 u = solver->velocity[v];
        world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(p, p + velocity_mul*u, 0.01, vec4(0,0,1,1));
    }
    for (auto e : geom->mesh.edges()) {
        vec3 a = eigen_to_vec3(geom->position[e.a().vertex()]);
        vec3 b = eigen_to_vec3(geom->position[e.b().vertex()]);
        world->graphics.paint.line(a,b,0.001,vec4(0,0,0,1));
        vec3 p = 0.5*a + 0.5*b;
        vec3 u = solver->velocity[e];
        world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(p, p + velocity_mul*u, 0.01, vec4(0,0,1,1));
    }

    for (auto tri : geom->mesh.faces()) {
        vec3 c = eigen_to_vec3(geom->barycenter(tri));
        vec3 n = solver->triangle_normal[tri];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
        // vec3 k = solver->triangle_projection_matrix[tri] * vec3(1,0,0);
        // world->graphics.paint.line(c,c+0.1*k,0.001,vec4(0.5,0,0.5,1));
    }
    for (auto e : geom->mesh.edges()) {
        vec3 c = eigen_to_vec3(geom->midpoint(e));
        vec3 n = solver->normal[e];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
    }
    for (auto v : geom->mesh.vertices()) {
        vec3 c = eigen_to_vec3(geom->position[v]);
        vec3 n = solver->normal[v];
        world->graphics.paint.line(c,c+0.1*n,0.001,vec4(0,1,0,1));
    }
}


void Demo::post_render_update()
{
}


void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
    }
}
