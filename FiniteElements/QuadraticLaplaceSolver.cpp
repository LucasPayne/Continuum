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

Aspect<Camera> main_camera;


vec3 eigen_to_vec3(Eigen::Vector3f v)
{
    return vec3(v.x(), v.y(), v.z());
}
Eigen::Vector3f vec3_to_eigen(vec3 v)
{
    return Eigen::Vector3f(v.x(), v.y(), v.z());
}




class Solver {
public:
    Solver(SurfaceGeometry &_geom);

    // VertexAttachment<double> solve();

    void set_source(PlaneFunction func);
    void set_dirichlet_boundary(PlaneFunction func);
private:
    SurfaceGeometry *geom; // Assumes y is 0. (this is bad, maybe should template the SurfaceGeometry class...)
                           // The geometry does not change.
    int num_interior_vertices;
    int num_boundary_vertices;
    VertexAttachment<int> vertex_indices;
    PlaneFunction dirichlet_boundary_function;
    PlaneFunction source_function;
};

Solver::Solver(SurfaceGeometry &_geom) :
    geom{&_geom},
    vertex_indices(_geom.mesh),
    source_function([](double,double)->double { return 0.0; }),
    dirichlet_boundary_function([](double,double)->double { return 0.0; })
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
    }
}
void Solver::set_source(PlaneFunction func)
{
    source_function = func;
}
void Solver::set_dirichlet_boundary(PlaneFunction func)
{
    dirichlet_boundary_function = func;
}




struct Demo : public IBehaviour {
    Demo(PlaneFunction dirichlet_boundary_function);

    SurfaceGeometry *geom;
    PlaneFunction dirichlet_boundary_function;

    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void update();
    void post_render_update();

    int mesh_N;
    bool random;
    double source_force;

    // Toggleable options.
    bool render_exact_solution;
    bool linear_mode;

    
    GLShaderProgram quadratic_mesh_shader;


    int screenshot_blx;
    int screenshot_bly;
    int screenshot_trx;
    int screenshot_try;
    float f_screenshot_blx;
    float f_screenshot_bly;
    float f_screenshot_trx;
    float f_screenshot_try;
};


void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) {
            exit(EXIT_SUCCESS);
        }
        if (e.key.code == KEY_O) {
            mesh_N -= 1;
            if (mesh_N < 2) mesh_N = 2;
        }
        if (e.key.code == KEY_P) {
            mesh_N += 1;
        }
        if (e.key.code == KEY_R) {
            random = !random;
        }
        if (e.key.code == KEY_Y) {
            render_exact_solution = !render_exact_solution;
        }
        if (e.key.code == KEY_Z) {
            linear_mode = !linear_mode;
        }
        static int counter = 0;
        std::string pre = linear_mode ? "linear_approx_" : "quadratic_approx_";
        if (e.key.code == KEY_T) {
            world->graphics.screenshot(pre + std::to_string(mesh_N) + ".ppm",
        			      screenshot_blx, screenshot_bly, screenshot_trx - screenshot_blx, screenshot_try - screenshot_bly);
            counter += 1;
        }
    }
}
void Demo::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
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

Demo::Demo(PlaneFunction _dirichlet_boundary_function) :
    dirichlet_boundary_function{_dirichlet_boundary_function}
{
    mesh_N = 2;
    random = false;
    geom = nullptr;
    source_force = 0.0;
    render_exact_solution = false;
    linear_mode = false;

    quadratic_mesh_shader.add_shader(GLShader(VertexShader, SHADERS "quadratic_mesh/quadratic_mesh.vert"));
    quadratic_mesh_shader.add_shader(GLShader(TessControlShader, SHADERS "quadratic_mesh/quadratic_mesh.tcs"));
    quadratic_mesh_shader.add_shader(GLShader(TessEvaluationShader, SHADERS "quadratic_mesh/quadratic_mesh.tes"));
    quadratic_mesh_shader.add_shader(GLShader(FragmentShader, SHADERS "quadratic_mesh/quadratic_mesh.frag"));
    quadratic_mesh_shader.link();

    screenshot_blx = 0;
    screenshot_bly = 0;
    screenshot_trx = 10;
    screenshot_try = 10;
    f_screenshot_blx = 0;
    f_screenshot_bly = 0;
    f_screenshot_trx = 0.01;
    f_screenshot_try = 0.01;
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
    // Recreate the mesh.
    if (geom != nullptr) delete geom;
    srand(230192301);
    if (render_exact_solution) geom = circle_mesh(50, random);
    else {
        if (linear_mode) {
            geom = circle_mesh(mesh_N*2, random);
        } else {
            geom = circle_mesh(mesh_N, random);
        }
    }
    // geom = square_mesh(mesh_N);
    
    // Edge midpoints.
    auto midpoints = EdgeAttachment<Eigen::Vector3f>(geom->mesh);
    for (auto edge : geom->mesh.edges()) {
        auto pa = geom->position[edge.a().vertex()];
        auto pb = geom->position[edge.b().vertex()];
        midpoints[edge] = 0.5*pa + 0.5*pb;
    }

    // struct FaceData {
    // };
    // auto face_data = FaceAttachment<FaceData>(geom->mesh);

    if (!render_exact_solution) {
        world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
    } else {
        // draw boundary circle
        auto ps = std::vector<vec3>();
        auto ps_lift = std::vector<vec3>();
        int circ_N = 300;
        for (int i = 0; i <= circ_N; i++) {
            float theta = i*2.f*M_PI/circ_N;
            ps.push_back(vec3(cos(theta), 0, sin(theta)));
            ps_lift.push_back(vec3(cos(theta), dirichlet_boundary_function(cos(theta), sin(theta)), sin(theta)));
        }
        world->graphics.paint.chain(ps, 0.001, vec4(0,0,0,1));
        world->graphics.paint.chain(ps_lift, 0.001, vec4(0,0,0,1));
    }
    
    auto vertex_u = VertexAttachment<double>(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        auto p = geom->position[v];
        vertex_u[v] = dirichlet_boundary_function(p.x(), p.z());
    }
    auto edge_u = EdgeAttachment<double>(geom->mesh);
    for (auto e : geom->mesh.edges()) {
        auto p = midpoints[e];
        edge_u[e] = dirichlet_boundary_function(p.x(), p.z());
    }
    

if (!render_exact_solution) {
    float linewid = 0.00035;
    if (!linear_mode) {
        for (auto edge : geom->mesh.edges()) {
            auto p = eigen_to_vec3(midpoints[edge]);
            auto p_lift = p + vec3(0,edge_u[edge],0);
            // world->graphics.paint.sphere(p + vec3(0,edge_u[edge],0), 0.02, vec4(0,1,0,1));
            world->graphics.paint.sphere(p, 0.0235, vec4(0.6,0.6,1.4,1));
            world->graphics.paint.sphere(p_lift, 0.018, vec4(0.6,0.6,1.4,1));
            // world->graphics.paint.line(p, p_lift, linewid, vec4(0,0,0,1));
        }
    }
    for (auto vertex : geom->mesh.vertices()) {
        auto p = eigen_to_vec3(geom->position[vertex]);
        auto p_lift = p + vec3(0,vertex_u[vertex],0);
        world->graphics.paint.sphere(p, 0.0235, vec4(1,0.5,0,1));
        world->graphics.paint.sphere(p_lift, 0.018, vec4(1,0.5,0,1));
        // world->graphics.paint.line(p, p_lift, linewid, vec4(0,0,0,1));
    }
}

    // Upload and render the quadratic surface mesh.
    //--------------------------------------------------------------------------------
    GLuint q_vao;
    glCreateVertexArrays(1, &q_vao);
    glBindVertexArray(q_vao);

    auto q_position_data = std::vector<vec3>(6*geom->mesh.num_faces());
    auto q_value_data = std::vector<float>(6*geom->mesh.num_faces());
    int tri_index = 0;
    for (auto tri : geom->mesh.faces()) {
        // each vertex
        auto start = tri.halfedge();
        auto he = start;
        int tri_vertex_index = 0;
        do {
            // vertex data
            q_position_data[6*tri_index + tri_vertex_index] = eigen_to_vec3(geom->position[he.vertex()]);
            q_value_data[6*tri_index + tri_vertex_index] = vertex_u[he.vertex()];
            // midpoint data
            q_position_data[6*tri_index + 3 + tri_vertex_index] = eigen_to_vec3(midpoints[he.edge()]);
            q_value_data[6*tri_index + 3 + tri_vertex_index] = edge_u[he.edge()];
            
            he = he.next();
            tri_vertex_index += 1;
        } while (he != start);
        tri_index += 1;
    }
    for (auto v : q_position_data) std::cout << v.x() << ", " << v.y() << ", " << v.z() << "\n";

    // int i = 0;
    // for (auto v : q_position_data) {
    //     world->graphics.paint.sphere(v + vec3(0,q_value_data[i],0), 0.05, vec4(1,1,0,1));
    //     i += 1;
    // }

    GLuint q_positions;
    glGenBuffers(1, &q_positions);
    glBindBuffer(GL_ARRAY_BUFFER, q_positions);
    glBufferData(GL_ARRAY_BUFFER, q_position_data.size()*sizeof(float)*3, &q_position_data[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *) 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    GLuint q_values;
    glGenBuffers(1, &q_values);
    glBindBuffer(GL_ARRAY_BUFFER, q_values);
    glBufferData(GL_ARRAY_BUFFER, q_value_data.size()*sizeof(float), &q_value_data[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (const void *) 0);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Render the quadratic surface.
    quadratic_mesh_shader.bind();
        mat4x4 mvp_matrix = main_camera->view_projection_matrix();

        glClear(GL_DEPTH_BUFFER_BIT); //...
        glEnable(GL_DEPTH_TEST);
        glUniformMatrix4fv(quadratic_mesh_shader.uniform_location("mvp_matrix"), 1, GL_FALSE, (const GLfloat *) &mvp_matrix);
        glUniform1i(quadratic_mesh_shader.uniform_location("linear_mode"), linear_mode ? 1 : 0);
        glPatchParameteri(GL_PATCH_VERTICES, 6);
        glDrawArrays(GL_PATCHES, 0, q_position_data.size());

    quadratic_mesh_shader.unbind();

    // Clean up.
    glDeleteVertexArrays(1, &q_vao);
    glDeleteBuffers(1, &q_positions);
    glDeleteBuffers(1, &q_values);
    //--------------------------------------------------------------------------------
    
if (!render_exact_solution) {
    // Draw the lines of the wireframe.
    int line_N = 25;
    for (auto edge : geom->mesh.edges()) {
        auto p = geom->position[edge.a().vertex()];
        auto pp = geom->position[edge.b().vertex()];
        auto midpoint = midpoints[edge];
        auto ps = std::vector<vec3>(line_N);

        // Evaluate basis functions restricted to this edge.
        float val_a = dirichlet_boundary_function(p.x(), p.z());
        float val_b = dirichlet_boundary_function(midpoint.x(), midpoint.z());
        float val_c = dirichlet_boundary_function(pp.x(), pp.z());

        for (int i = 0; i < line_N; i++) {
            auto x = p.x() + (i*1.f/(line_N-1)) * (pp - p).x();
            auto z = p.z() + (i*1.f/(line_N-1)) * (pp - p).z();

            float u = i*1.f/(line_N-1);
        
            float val;
            if (linear_mode) {
                val = val_a * (1-u) + val_c * u;
            } else {
                val =   val_a * (1 - u - 2*u*(1-u))
                            + val_b * (4*u*(1-u))
                            + val_c * (u - 2*u*(1-u));
            }

            ps[i] = vec3(x, val + 0.003, z);
        }
        world->graphics.paint.chain(ps, 0.0021, vec4(0,0,0,1));
    }
}

}



class App : public IGC::Callbacks {
public:
    World &world;
    App(World &world);

    void close();
    void loop();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void window_handler(WindowEvent e);

    SurfaceGeometry *geom; // Delete mesh and geom then reset this pointer when changing mesh.

    Entity demo_e;
    Aspect<Behaviour> demo_b;
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
    
    // Create the demo.
    demo_e = world.entities.add();
    auto demo = world.add<Demo>(demo_e, [](double x, double y)->double {
        // return x*x - y*y + 1.13;
        // return 0.5*sin(5*x)*y + 1.13;
        // return 0.5*(x*x - y*y) + 0.2*sin(8*x)*y + 1.13;
        // return 0.5*(x*x - y*y)+1.13 + cos(8*x)*0.28;
        
        return 0.5*(x*x - y*y)+1.23 + cos(4*x)*0.3;
        
        // return 0.5*cos(8*x)*y + exp(-x*x)*0.5*y*y+1.3;
    });
    demo_b = demo_e.get<Behaviour>();

}

void App::close()
{
}

void App::loop()
{
}

void App::window_handler(WindowEvent e)
{
}

void App::keyboard_handler(KeyboardEvent e)
{
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
