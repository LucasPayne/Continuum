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


class HeatSolver {
public:
    HeatSolver(SurfaceGeometry &_geom, PlaneFunction u_initial, double _diffusivity);

    void step(double delta_time);

    void set_source(PlaneFunction func);

    VertexAttachment<double> u_prev;
private:
    SurfaceGeometry *geom; // Assumes y is 0. (this is bad, maybe should template the SurfaceGeometry class...)
                           // The geometry does not change.
    int num_interior_vertices;
    int num_boundary_vertices;
    VertexAttachment<int> vertex_indices;
    PlaneFunction dirichlet_boundary_function; // Interpolated from the initial condition.
    PlaneFunction source_function;

    double time;

    double diffusivity;
};

HeatSolver::HeatSolver(SurfaceGeometry &_geom, PlaneFunction u_initial, double _diffusivity) :
    geom{&_geom},
    vertex_indices(_geom.mesh),
    source_function([](double,double)->double { return 0.0; }),

    u_prev(_geom.mesh),
    time{0.0},

    diffusivity{_diffusivity}
{
    for (auto v : geom->mesh.vertices()) {
        auto p = geom->position[v];
        u_prev[v] = u_initial(p.x(), p.z());
    }
    dirichlet_boundary_function = u_initial; // Inherit the Dirichlet boundary condition from initial condition.

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
void HeatSolver::set_source(PlaneFunction func)
{
    source_function = func;
}

//note: "dt" is a global application time-step, do not use it here.
void HeatSolver::step(double delta_time)
{
    /*--------------------------------------------------------------------------------
        a(u, psi) =   < u/delta_time, psi >
		    + < diffusivity * grad(u), grad(psi) >
        f(psi)    =   < g + u_prev/delta_time, psi >
    --------------------------------------------------------------------------------*/
    std::vector<EigenTriplet> coefficients;
    auto rhs = Eigen::VectorXd(num_interior_vertices);

    // Precomputations.
    double inv_delta_time = 1.0 / delta_time;

    /*--------------------------------------------------------------------------------
        Construct the system.
    --------------------------------------------------------------------------------*/
    // Construct the RHS from interior source terms.
    // (Boundary conditions will later modify this RHS).
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;

        double total_triangle_area = 0.0;
        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            total_triangle_area += geom->triangle_area(tri);
            he = he.twin().next();
        } while (he != start);
        double val = source_function(geom->position[v].x(), geom->position[v].z());
        val += u_prev[v] * inv_delta_time;
        rhs[vertex_indices[v]] = (val/3) * total_triangle_area;
    }
    
    // Construct the mass matrix.
    for (auto tri : geom->mesh.faces()) {
        double triangle_area = geom->triangle_area(tri);

        // Get the triangle vertices (todo: provide mesh_processing function for this).
        Vertex verts[3];
        int vertex_index = 0;
        auto he = tri.halfedge();
        do {
            he = he.next();
            verts[vertex_index] = he.vertex();
        } while (++vertex_index < 3);

        // Compute the local stiffness matrix for the gradient inner product.
        Eigen::Matrix<double, 2,2> jacobian_matrix;
        auto A = geom->position[verts[0]];
        auto B = geom->position[verts[1]];
        auto C = geom->position[verts[2]];
        jacobian_matrix <<
            B.x() - A.x(),  C.x() - A.x(),
            B.z() - A.z(),  C.z() - A.z();
        auto grad_transform = jacobian_matrix.inverse().transpose();
        Eigen::Matrix<double, 3,2> gradients;
        gradients <<
            -1,-1,
            1,0,
            0,1;
        Eigen::Matrix<double, 3,3> local_stiffness_matrix;
        for (int i = 0; i < 3; i++) {
            auto grad_i = grad_transform * gradients.row(i).transpose();
            for (int j = 0; j < 3; j++) {
                auto grad_j = grad_transform * gradients.row(j).transpose();
                local_stiffness_matrix(i,j) = diffusivity * triangle_area * grad_i.dot(grad_j); //??????????????????????????????????????
            }
        }
        // Insert this stiffness matrix into the global mass matrix.
        // Boundary conditions:
        //     If the triangle is on a boundary, then a trial function centered on a boundary vertex has a predetermined (Dirichlet condition) coefficient,
        //     so the RHS is edited instead of editing the mass matrix.
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                bool b1 = verts[i].on_boundary();
                bool b2 = verts[j].on_boundary();
                if (b1) {
                    continue; // The trial function can't be centered on the boundary.
                } else if (b2) {
                    auto p = geom->position[verts[j]];
                    double boundary_val = dirichlet_boundary_function(p.x(), p.z());
                    rhs[vertex_indices[verts[i]]] -= boundary_val*local_stiffness_matrix(i,j);
                } else {
                    coefficients.push_back(EigenTriplet(vertex_indices[verts[i]], vertex_indices[verts[j]], local_stiffness_matrix(i,j)));
                }
            }
        }
    }
    // Modify the mass matrix to include the < u/delta_time, psi > term.
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;

        int j = vertex_indices[v];

        double total_triangle_area = 0.0;
        // For each adjacent triangle.
        auto start = v.halfedge(); //todo: provide this circulator in mesh_processing.
        auto he = start;
        do {
            auto tri = he.face();
            total_triangle_area += geom->triangle_area(tri);
            he = he.twin().next();
        } while (he != start);

        double val = total_triangle_area * inv_delta_time/3.0;
        coefficients.push_back(EigenTriplet(j, j, val));
    }
    // Finalize the mass matrix.
    auto mass_matrix = SparseMatrix(num_interior_vertices, num_interior_vertices);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();

    /*--------------------------------------------------------------------------------
        Solve the system.
    --------------------------------------------------------------------------------*/

    // std::cout << Eigen::MatrixXd(mass_matrix) << "\n";
    // std::cout << rhs << "\n";
    // getchar();
    //--------------------------------------------------------------------------------
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    //--------------------------------------------------------------------------------
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(mass_matrix);
    // Compute the numerical factorization 
    solver.factorize(mass_matrix);
    // Use the factors to solve the linear system 
    Eigen::VectorXd u = solver.solve(rhs);

    // Update u_prev.
    // Reassociate each coefficient (or boundary value) with the corresponding vertex of the mesh.
    int interior_vertex_index = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            auto p = geom->position[v];
            u_prev[v] = dirichlet_boundary_function(p.x(), p.z());
        } else {
            u_prev[v] = u[interior_vertex_index];
            interior_vertex_index += 1;
        } 
    }
}





struct HeatDemo : public IBehaviour {
    HeatDemo(SurfaceGeometry &_geom, PlaneFunction initial_u, double diffusivity);

    SurfaceGeometry *geom;
    HeatSolver heat_solver;

    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void update();
    void post_render_update();
};


void HeatDemo::keyboard_handler(KeyboardEvent e)
{
}
void HeatDemo::mouse_handler(MouseEvent e)
{
}

HeatDemo::HeatDemo(SurfaceGeometry &_geom, PlaneFunction initial_u, double diffusivity) :
    geom{&_geom},
    heat_solver(_geom, initial_u, diffusivity)
{
}

void HeatDemo::update()
{
    heat_solver.step(dt);
}

void HeatDemo::post_render_update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);

    auto lifted_geom = SurfaceGeometry(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        double u_val = heat_solver.u_prev[v];
        auto geom_pos = geom->position[v];
        lifted_geom.position[v] = Eigen::Vector3f(geom_pos.x(), u_val, geom_pos.z());
    }
    world->graphics.paint.wireframe(lifted_geom, mat4x4::identity(), 0.001);


    vec4 color = vec4(0.1,0.1,0.5,1);
    auto lifted_boundary_positions = std::vector<vec3>();
    for (auto v : geom->mesh.vertices()) {
        if (!v.on_boundary()) continue;
        auto p = geom->position[v];
        vec3 pp = *((vec3 *) &p[0]);
        world->graphics.paint.sphere(pp, 0.01, color);
        p = lifted_geom.position[v];
        pp = *((vec3 *) &p[0]);
        lifted_boundary_positions.push_back(pp);
    }
    lifted_boundary_positions.push_back(lifted_boundary_positions[0]);
    world->graphics.paint.chain(lifted_boundary_positions, 4, color);
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
    auto demo = world.add<HeatDemo>(demo_e, *geom, [](double x, double y)->double {
        // return 1+0.4*sin(8*x + total_time);
        float r = 0.5;
        if (x*x + y*y <= r*r) return 1.0;
        return 0.0;
    }, 0.01);
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

        // delete &geom->mesh;
        // delete geom;
        // geom = circle_mesh(mesh_N, false);

        // demo_e.remove<Behaviour>();
        // world.add<HeatDemo>(demo_e, *geom, [](double x, double y)->double {
        //     // return 1+0.4*sin(8*x + total_time);
        //     float r = 0.5;
        //     if (x*x + y*y <= r*r) return 1.0;
        //     return 0.0;
        // }, 0.01);

        // note: The entity system is incomplete...
        for (auto b : world.entities.aspects<Behaviour>()) {
            if (b.entity() == demo_e) b->enabled = false;
        }
        geom = circle_mesh(mesh_N, random);
        world.add<HeatDemo>(demo_e, *geom, [](double x, double y)->double {
            // return 1+0.4*sin(8*x + total_time);
            
            double t = 0.0;
            
            // t += 1+0.4*sin(8*x + total_time);
            t += 1+0.4*sin(8*x + y);

            double r = 0.5;
            if (x*x + y*y <= r*r) t = 2;
            return t;

        }, 0.01);
            

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
