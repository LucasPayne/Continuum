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
    
    geom = square_mesh(60);
    for (auto v : geom->mesh.vertices()) {
        vec3 p = eigen_to_vec3(geom->position[v]);
        float theta = 0.95*(p.x()+1)*M_PI;
        // geom->position[v] += Eigen::Vector3f(0, 0.1*sin(4*p.x()), 0);
        // geom->position[v] = Eigen::Vector3f(cos(theta), sin(theta), p.z());
    }
    
    double viscosity = 0.001;
    solver = new SurfaceNavierStokesSolver(*geom, viscosity);

    solver->set_source([&](double x, double y, double z)->vec3 {
        #if 1
        double r = 0.125;
        if ((x+0.8)*(x+0.8) + z*z <= r*r) return vec3(50,0,0);
        // if (y*y + z*z <= r*r) return vec3(1,0,1);
        return vec3(0,0,0);
        #else
        double r = 0.25;
        if (x*x + z*z <= r*r) return vec3(1,-2,0);
        return vec3(0,0,0);
        #endif
    });
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) exit(EXIT_SUCCESS);
        if (e.key.code == KEY_R) {
            solver->time_step(0.01);
        }
        if (e.key.code == KEY_P) {
            solver->m_current_time_step_dt = 0.025;
            solver->explicit_advection();
        }
        if (e.key.code == KEY_1) {
            solver->set_velocity([&](double x, double y, double z)->vec3 {
                double r = 0.125;
                if (x*x + z*z <= r*r) return vec3(-2,0,-2);
                return vec3(0,0,0);
            });
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
        // world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(p, p + velocity_mul*u, 0.01, vec4(0,0,1,1));
    }
    for (auto e : geom->mesh.edges()) {
        vec3 a = eigen_to_vec3(geom->position[e.a().vertex()]);
        vec3 b = eigen_to_vec3(geom->position[e.b().vertex()]);
        world->graphics.paint.line(a,b,0.001,vec4(0,0,0,1));
        vec3 p = 0.5*a + 0.5*b;
        vec3 u = solver->velocity[e];
        // world->graphics.paint.sphere(p, 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(p, p + velocity_mul*u, 0.01, vec4(0,0,1,1));
    }
    #if 0
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
    #endif

    #if 0
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        vec3 c = eigen_to_vec3(geom->position[v]);
        vec3 r =  solver->_test_point_1[v];
        //world->graphics.paint.line(c,c+100*(r-c),0.005,vec4(0,0,0,1));
        world->graphics.paint.line(c,r,0.005,vec4(0,0,0,1));
        world->graphics.paint.sphere(r, 0.01, vec4(1,0,0,0));
        // world->graphics.paint.line(c,eigen_to_vec3(geom->barycenter(solver->_test_point_2[v])),0.005,vec4(1,0,0,1));
        // world->graphics.paint.line(c,solver->_test_point_2[v],0.005,vec4(1,0,0,1));
    }
    for (auto e : geom->mesh.edges()) {
        if (e.on_boundary()) continue;
        vec3 c = eigen_to_vec3(geom->midpoint(e));
        vec3 r =  solver->_test_point_1[e];
        // world->graphics.paint.line(c,r,0.005,vec4(0,0,0,1));
        // world->graphics.paint.line(c,solver->_test_point_2[e],0.005,vec4(1,0,0,1));
    }
    #endif

    solver->time_step(1./300.);
    static int counter = 0;
    save_solution(std::string(DATA) + "sns." + std::to_string(counter) + ".txt");
}


void Demo::post_render_update()
{
}


void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
    }
}

void Demo::save_solution(std::string filename)
{
    FILE *file = fopen(filename.c_str(), "w+");
    for (auto v : geom->mesh.vertices()) {
        fprintf(file, "%.7f %.7f\n", solver->velocity[v].x(), solver->velocity[v].y());
    }
    for (auto e : geom->mesh.edges()) {
        fprintf(file, "%.7f %.7f\n", solver->velocity[e].x(), solver->velocity[e].y());
    }
    for (auto v : geom->mesh.vertices()) {
        fprintf(file, "%.7g\n", solver->pressure[v]);
    }

    fclose(file);
}
