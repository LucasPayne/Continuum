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


class FVSolver {
public:
    FVSolver(SurfaceGeometry &_geom);

    VertexAttachment<double> solve();

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

FVSolver::FVSolver(SurfaceGeometry &_geom) :
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
void FVSolver::set_source(PlaneFunction func)
{
    source_function = func;
}
void FVSolver::set_dirichlet_boundary(PlaneFunction func)
{
    dirichlet_boundary_function = func;
}
VertexAttachment<double> FVSolver::solve()
{
    auto rhs = Eigen::VectorXd(num_interior_vertices);
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        // Add the value to entry (i,j) in the resulting matrix.
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    /*--------------------------------------------------------------------------------
        Construct the system.
    --------------------------------------------------------------------------------*/
    for (Vertex v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom->position[v];

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];
            double C = 1.0/(4.0*geom->triangle_area(tri));

            // Diagonal term.
            double diagonal_term = C*(vpp_pos - vp_pos).dot(vpp_pos - vp_pos);
            add_entry(v_index, v_index, diagonal_term);

            // vp contribution.
            if (vp.on_boundary()) {
                double bv = dirichlet_boundary_function(vp_pos.x(), vp_pos.z());
                rhs[v_index] -= bv * C * (v_pos - vpp_pos).dot(vpp_pos - vp_pos);
            } else {
                double vp_term = C*(v_pos - vpp_pos).dot(vpp_pos - vp_pos);
                int vp_index = vertex_indices[vp];
                add_entry(v_index, vp_index, vp_term);
            }
            
            // vpp contribution.
            if (vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(vpp_pos.x(), vpp_pos.z());
                rhs[v_index] -= bv * C * (vp_pos - v_pos).dot(vpp_pos - vp_pos);
            } else {
                double vpp_term = C*(vp_pos - v_pos).dot(vpp_pos - vp_pos);
                int vpp_index = vertex_indices[vpp];
                add_entry(v_index, vpp_index, vpp_term);
            }

            he = he.twin().next();
        } while (he != start);
    }
    

    // Finalize the mass matrix.
    auto mass_matrix = SparseMatrix(num_interior_vertices, num_interior_vertices);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();

    /*--------------------------------------------------------------------------------
        Solve the system.
    --------------------------------------------------------------------------------*/
    //--------------------------------------------------------------------------------
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    //--------------------------------------------------------------------------------
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(mass_matrix);
    // Compute the numerical factorization 
    solver.factorize(mass_matrix);
    // Use the factors to solve the linear system 
    Eigen::VectorXd u = solver.solve(rhs);
    
    /*--------------------------------------------------------------------------------
    // Reassociate each coefficient (or boundary value) with the corresponding vertex of the mesh.
    --------------------------------------------------------------------------------*/
    auto mesh_u = VertexAttachment<double>(geom->mesh);
    int interior_vertex_index = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            auto p = geom->position[v];
            mesh_u[v] = dirichlet_boundary_function(p.x(), p.z());
        } else {
            mesh_u[v] = u[interior_vertex_index];
            interior_vertex_index += 1;
        } 
    }
    return mesh_u;
}





struct FVDemo : public IBehaviour {
    FVDemo(SurfaceGeometry &_geom, PlaneFunction dirichlet_boundary_function);

    SurfaceGeometry *geom;
    FVSolver fv_solver;

    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void update();
    void post_render_update();
};


void FVDemo::keyboard_handler(KeyboardEvent e)
{
}
void FVDemo::mouse_handler(MouseEvent e)
{
}

FVDemo::FVDemo(SurfaceGeometry &_geom, PlaneFunction dirichlet_boundary_function) :
    geom{&_geom},
    fv_solver(_geom)
{
    fv_solver.set_dirichlet_boundary(dirichlet_boundary_function);
}

void FVDemo::update()
{
}

void FVDemo::post_render_update()
{
    auto solver = FVSolver(*geom);
    solver.set_dirichlet_boundary([](double x, double y)->double {
        return exp(-1.6*x*x - 0.3*(y-1)*(y-1));
    });
    VertexAttachment<double> mesh_u = solver.solve();


    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
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

    Aspect<Camera> main_camera;

    SurfaceGeometry *geom; // Delete mesh and geom then reset this pointer when changing mesh.

    Entity demo_e;
    Aspect<Behaviour> demo_b;
    float source_force;
    int mesh_N;
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
    source_force = 20;
    mesh_N = 5;
    geom = circle_mesh(mesh_N, false);
    demo_e = world.entities.add();
    auto demo = world.add<FVDemo>(demo_e, *geom, [](double x, double y)->double {
        return 1+0.4*sin(8*x + total_time);
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
    auto reset = [&](bool random) {
        // note: The entity system is incomplete...
        for (auto b : world.entities.aspects<Behaviour>()) {
            if (b.entity() == demo_e) b->enabled = false;
        }
        geom = circle_mesh(mesh_N, random);
        world.add<FVDemo>(demo_e, *geom, [](double x, double y)->double {
            return 1+0.4*sin(8*x + total_time);
        });
    };
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) {
            exit(EXIT_SUCCESS);
        }
        if (e.key.code == KEY_P) {
            mesh_N += 5;
            reset(false);
        }
        if (e.key.code == KEY_R) {
            reset(true);
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
