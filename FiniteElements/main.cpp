#include "cg_sandbox.h"
typedef float REAL;
typedef void VOID;
extern "C" {
#include "triangle.h"
}

// #ifdef __cplusplus
// #endif


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
    std::vector<vec2> interior_points;
    std::vector<vec2> boundary_points;

};


App::App(World &_world) : world{_world}
{
    auto cameraman = world.entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    cameraman.add<Transform>(0,0,2);
    main_camera = camera;
    // world.add<CameraController>(cameraman);
    
    // Create a roughly regular point cloud for the circle, and points sampled on the boundary.
    int N = 10;
    for (int i = 0; i <= N; i++) {
        float x = -1+(i*2.f)/N;
        for (int j = 0; j <= N; j++) {
            float y = -1+(j*2.f)/N;
            if (x*x + y*y <= 1-1.f/N) interior_points.push_back(vec2(x,y));
        }
    }
    int num_angle_intervals = 2*N;
    for (int i = 0; i <= num_angle_intervals; i++) {
        float theta = (i*2.f*M_PI)/num_angle_intervals;
        boundary_points.push_back(vec2(cos(theta), sin(theta)));
    }
    for (auto *ps : {&interior_points, &boundary_points}) {
        for (int i = 0; i < ps->size(); i++) {
            (*ps)[i] = 0.2f*(*ps)[i] + vec2(0.5,0.5);
        }
    }
    
    std::vector<vec2> all_points;
    all_points.insert(all_points.end(), interior_points.begin(), interior_points.end());
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());

    // Triangulate this point cloud.
    std::string triswitches = "zQ";
    struct triangulateio in = {0};
    struct triangulateio out = {0};

    in.pointlist = (float *) &all_points[0];
    in.numberofpoints = all_points.size();
    in.numberofpointattributes = 0;
    triangulate(&triswitches[0], &in, &out, nullptr);

}

void App::close()
{
}

void App::loop()
{
    world.graphics.paint.circles(main_camera, boundary_points, 0.005, vec4(0,0,0,1));
    world.graphics.paint.circles(main_camera, interior_points, 0.005, vec4(0.5,0.5,0.5,1));
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
    }
}

void App::mouse_handler(MouseEvent e)
{
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
