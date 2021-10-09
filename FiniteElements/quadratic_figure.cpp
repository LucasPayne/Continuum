#include "cg_sandbox.h"
#include "triangle_wrapper.h"
#include "CameraController.cpp"
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "mesh_generators.cpp"

using PlaneFunction = std::function<double(double x, double y)>; // Function of the XY plane.
using PlaneFunctionNL1 = std::function<double(double x, double y, double u)>; // First-order non-linear plane function.
using SparseMatrix = Eigen::SparseMatrix<double>;
using EigenTriplet = Eigen::Triplet<double>;

vec3 eigen_to_vec3(Eigen::Vector3f v)
{
    return vec3(v.x(), v.y(), v.z());
}
Eigen::Vector3f vec3_to_eigen(vec3 v)
{
    return Eigen::Vector3f(v.x(), v.y(), v.z());
}

#define num_funcs 6
std::function<double(double,double,double)> funcs[num_funcs] = {
    [](double x, double y, double z)->double {
        return x - 2*x*y - 2*x*z;
    },
    [](double x, double y, double z)->double {
        return y - 2*y*z - 2*y*x;
    },
    [](double x, double y, double z)->double {
        return z - 2*z*x - 2*z*y;
    },
    [](double x, double y, double z)->double {
        return 4*x*y;
    },
    [](double x, double y, double z)->double {
        return 4*y*z;
    },
    [](double x, double y, double z)->double {
        return 4*z*x;
    }
};
std::string func_names[num_funcs] = {
    "1_100",
    "2_010",
    "3_001",
    "4_110",
    "5_011",
    "6_101"
};


class App : public IGC::Callbacks {
public:
    World &world;
    App(World &world);

    void close();
    void loop();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void window_handler(WindowEvent e);

    Aspect<Camera> main_camera;


    vec3 tri_points[3];

    int func_index;

    int screenshot_blx;
    int screenshot_bly;
    int screenshot_trx;
    int screenshot_try;
    float f_screenshot_blx;
    float f_screenshot_bly;
    float f_screenshot_trx;
    float f_screenshot_try;
};



App::App(World &_world) : world{_world}
{
    // OpenGL ----------------------------------------------------
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //------------------------------------------------------------

    auto cameraman = world.entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    camera->background_color = vec4(1,1,1,1);
    auto t = cameraman.add<Transform>(0,0,2);
    main_camera = camera;
    auto controller = world.add<CameraController>(cameraman);
    t->position = vec3(0,1,3);
    controller->angle = -2;

    func_index = 0;

    screenshot_blx = 0;
    screenshot_bly = 0;
    screenshot_trx = 10;
    screenshot_try = 10;
}

void App::close()
{
}

void App::loop()
{
    auto tri_points_vec = std::vector<vec3>(4);
    for (int i = 0; i < 4; i++){
        tri_points_vec[i] = vec3(cos(i*2.f*M_PI/3.f), 0, sin(i*2.f*M_PI/3.f));
        if (i < 3) tri_points[i] = tri_points_vec[i];
    }
    world.graphics.paint.chain(tri_points_vec, 0.006, vec4(0,0,0,1));

    auto tri_point = [&](double x, double y, double z)->vec3 {
        return x*tri_points[0] + y*tri_points[1] + z*tri_points[2];
    };

    SurfaceMesh mesh;
    SurfaceGeometry geom(mesh);
    VertexAttachment<vec3> barycentric(mesh);
    int N = 100;
    auto vertex_grid = std::vector<Vertex>((N+1)*(N+1));
    for (int i = 0; i <= N; i++) {
        double x = i*1./N;
        for (int j = 0; j <= N-i; j++) {
            int k = N - i - j;
            float y = j*1./N;
            float z = k*1./N;
            // printf("%d %d %.2f %.2f %.2f\n", i,j,x,y,z);
            auto v = mesh.add_vertex();
            vec3 p = tri_point(x,y,z);
            geom.position[v] = Eigen::Vector3f(p.x(), p.y(), p.z());
            vertex_grid[i*(N+1) + j] = v;
            barycentric[v] = vec3(x,y,z);
        }
    }
    for (int i = 1; i <= N; i++) {
        for (int j = 0; j <= N-i; j++) {
            auto vgrid = [&](int a, int b) { return vertex_grid[(i+a)*(N+1) + j+b]; };
            mesh.add_triangle(vgrid(0,0), vgrid(-1,1), vgrid(-1,0));
            if (j > 0) {
                mesh.add_triangle(vgrid(0,0), vgrid(-1,0), vgrid(0,-1));
            }
        }
    }
    mesh.lock();
    // world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.001);
    for (auto tri : mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        auto ps = std::vector<vec3>();
        do {
            auto p = geom.position[he.vertex()];
            ps.push_back(vec3(p.x(), p.y(), p.z()));
            he = he.next();
        } while (he != start);
        ps.push_back(ps[0]);
        // world.graphics.paint.chain(ps, 0.0005, vec4(0,0,0,0.5));
    }

    auto func = funcs[func_index];
    SurfaceGeometry lifted_geom(mesh);
    for (auto v : mesh.vertices()) {
        auto b = barycentric[v];
        lifted_geom.position[v] = geom.position[v] + Eigen::Vector3f(0, func(b.x(), b.y(), b.z()), 0);
    }
    world.graphics.paint.wireframe(lifted_geom, mat4x4::identity(), 0.001);
    for (auto he : mesh.halfedges()) {
        // if (!edge.a().vertex().on_boundary() || !edge.b().vertex().on_boundary()) continue;
        if (!he.face().null()) continue;
        vec3 a = eigen_to_vec3(lifted_geom.position[he.vertex()]);
        vec3 b = eigen_to_vec3(lifted_geom.position[he.tip()]);
        world.graphics.paint.line(a, b, 0.005, vec4(0,0,0,1));
    }



    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2-i; j++) {
            double x = i*1./2;
            double y = j*1./2;
            vec3 p = tri_point(x,y, 1-x-y);
            world.graphics.paint.sphere(p, 0.03, vec4(1.4,0.7,0.7,1));
            world.graphics.paint.line(p, p + vec3(0,func(x,y,1-x-y),0), 0.005, vec4(0,0,0,1));
        }
    }

    // for (auto v : mesh.vertices()) {
    //     auto p = geom.position[v];
    //     vec3 pp = vec3(p.x(), p.y(), p.z());
    //     world.graphics.paint.sphere(pp, 0.02, vec4(1,0,0,1));
    // }
    //
    if (world.input.keyboard.down(KEY_R)) {
        // Draw the screenshot rectangle.
        std::vector<vec2> ps = {
            vec2(f_screenshot_blx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_bly)
        };
        world.graphics.paint.chain_2D(ps, 1, vec4(1,0,0,1));
    }

}

void App::window_handler(WindowEvent e)
{
}

void App::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) {
            exit(EXIT_SUCCESS);
        }
        if (e.key.code == KEY_M) {
            func_index = (func_index + 1) % num_funcs;
        }
        if (e.key.code == KEY_N) {
            func_index -= 1;
            if (func_index < 0) func_index = num_funcs - 1;
        }
        if (e.key.code == KEY_T) {
            world.graphics.screenshot("quadratic_basis_" + func_names[func_index] + ".ppm",
                                      screenshot_blx, screenshot_bly, screenshot_trx - screenshot_blx, screenshot_try - screenshot_bly);
        }
    }
}

void App::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
        if (e.button.code == MOUSE_LEFT) {
            screenshot_blx = int(world.graphics.window_viewport.w * e.cursor.x);
            screenshot_bly = int(world.graphics.window_viewport.h * e.cursor.y);
            f_screenshot_blx = e.cursor.x;
            f_screenshot_bly = e.cursor.y;
        }
        if (e.button.code == MOUSE_RIGHT) {
            screenshot_trx = int(world.graphics.window_viewport.w * e.cursor.x);
            screenshot_try = int(world.graphics.window_viewport.h * e.cursor.y);
            f_screenshot_trx = e.cursor.x;
            f_screenshot_try = e.cursor.y;
        }
    }
}



int main(int argc, char *argv[])
{
    printf("[main] Creating context...\n");
    IGC::Context context("A world");
    printf("[main] Creating world...\n");
    World world(context);

    printf("[main] Creating app...\n");
    App app(world);
    printf("[main] Adding app callbacks...\n");
    context.add_callbacks(&app);

    printf("[main] Entering loop...\n");
    context.enter_loop();
    context.close();
}
