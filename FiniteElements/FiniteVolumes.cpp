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



vec2 triangle_circumcenter(vec2 p, vec2 q, vec2 r)
{
    float pp = vec2::dot(p, p);
    float qq = vec2::dot(q, q);
    float rr = vec2::dot(r, r);
    vec2 pd = (rr - qq)*vec2(-p.y(), p.x());
    vec2 qd = (pp - rr)*vec2(-q.y(), q.x());
    vec2 rd = (qq - pp)*vec2(-r.y(), r.x());
    float denom = vec2::cross(p, q) + vec2::cross(q, r) + vec2::cross(r, p);
    return 0.5*(pd + qd + rd)/denom;
}

bool triangle_contains_point(vec2 p, vec2 q, vec2 r, vec2 point)
{
    bool p_side = vec2::dot(point - p, (p - q).perp()) > 0;
    bool q_side = vec2::dot(point - q, (q - r).perp()) > 0;
    bool r_side = vec2::dot(point - r, (r - p).perp()) > 0;
    return (p_side == q_side) && (q_side == r_side);
}



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
    
    /*--------------------------------------------------------------------------------
        Source term contribution.
    --------------------------------------------------------------------------------*/
    #if 0
    // Rectangle quadrature
    for (Vertex v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom->position[v];

        // Compute the control volume area.
        double control_volume_area = 0.0;
        
        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            // Get the triangle vertex positions.
            //     todo: This should be simple to do with mesh_processing.
            Eigen::Vector3f ps_3d[3];
            int vertex_index = 0;
            auto he_e = he;
            do {
                he_e = he_e.next();
                ps_3d[vertex_index] = geom->position[he_e.vertex()];
            } while (++vertex_index < 3);
            // Convert to 2D points.
            vec2 ps[3];
            for (int i = 0; i < 3; i++) {
                ps[i] = vec2(ps_3d[i].x(), ps_3d[i].z());
            }
            
            vec2 c = triangle_circumcenter(ps[0], ps[1], ps[2]);
            double segment_area = 0.0;
            if (triangle_contains_point(ps[0], ps[1], ps[2], c)) {
                vec2 midpoint1 = 0.5*(ps[0] + ps[1]);
                vec2 midpoint2 = 0.5*(ps[0] + ps[2]);
                // The circumcenter is inside the triangle.
                // The resulting segment consists of one quarter triangle,
                // and the triangle joining the midpoints with the circumcenter.
                segment_area = 0.25*geom->triangle_area(tri) + 0.5*vec2::cross(c - midpoint1, midpoint2 - midpoint1);
            } else {
                // The circumcenter is pulled back to the midpoint of the opposite edge.
                // The resulting segment consists of two quarter triangles.
                segment_area = 0.5*geom->triangle_area(tri);
            }
            control_volume_area += segment_area;

            he = he.twin().next();
        } while (he != start);

        // Rectangle-rule quadrature for source contribution.
        // (This takes one sample at the control volume center,
        //  then assumes the function is constant on the control volume.)
        rhs[v_index] += source_function(v_pos.x(), v_pos.z()) * control_volume_area;
    }
    #else
    // Linear quadrature.
    for (Vertex v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = vertex_indices[v]; // Global interior vertex index.
        auto v_pos = geom->position[v];

        double integral = 0.0;
        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            // Get the triangle vertex positions.
            //     todo: This should be simple to do with mesh_processing.
            Eigen::Vector3f ps_3d[3];
            int vertex_index = 0;
            auto he_e = he;
            do {
                he_e = he_e.next();
                ps_3d[vertex_index] = geom->position[he_e.vertex()];
            } while (++vertex_index < 3);
            // Convert to 2D points.
            vec2 ps[3];
            for (int i = 0; i < 3; i++) {
                ps[i] = vec2(ps_3d[i].x(), ps_3d[i].z());
            }
            // Compute the circumcenter.
            vec2 c = triangle_circumcenter(ps[0], ps[1], ps[2]);

	    vec2 midpoint1 = 0.5*(ps[0] + ps[1]);
	    vec2 midpoint2 = 0.5*(ps[0] + ps[2]);
            
            double f_v = source_function(v_pos.x(), v_pos.z());
            double f_midpoint1 = source_function(midpoint1.x(), midpoint1.y());
            double f_midpoint2 = source_function(midpoint2.x(), midpoint2.y());
            double f_c = source_function(c.x(), c.y());
            double tri_area = geom->triangle_area(tri);
            double c_tri_area = 0.5*vec2::cross(c - midpoint1, midpoint2 - midpoint1);

            integral += 0.5*tri_area*(f_v + f_midpoint1 + f_midpoint2)/6.0;
            integral += 2*c_tri_area*(f_midpoint1 + f_midpoint2 + f_c)/6.0;

            he = he.twin().next();
        } while (he != start);

        rhs[v_index] += integral;
    }
    #endif
    

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
    FVDemo(PlaneFunction dirichlet_boundary_function);

    SurfaceGeometry *geom;
    PlaneFunction dirichlet_boundary_function;

    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void update();
    void post_render_update();

    int mesh_N;
    bool random;
    double source_force;
};


void FVDemo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_Q) {
            exit(EXIT_SUCCESS);
        }
        if (e.key.code == KEY_O) {
            mesh_N -= 5;
            if (mesh_N < 2) mesh_N = 2;
        }
        if (e.key.code == KEY_P) {
            mesh_N += 5;
        }
        if (e.key.code == KEY_R) {
            random = !random;
        }
    }
}
void FVDemo::mouse_handler(MouseEvent e)
{
}

FVDemo::FVDemo(PlaneFunction _dirichlet_boundary_function) :
    dirichlet_boundary_function{_dirichlet_boundary_function}
{
    mesh_N = 5;
    random = false;
    geom = nullptr;
    source_force = 0.0;
}

void FVDemo::update()
{
    float source_force_change_speed = 20.f;
    if (world->input.keyboard.down(KEY_M)) {
        source_force += source_force_change_speed * dt;
    }
    if (world->input.keyboard.down(KEY_N)) {
        source_force -= source_force_change_speed * dt;
    }
}

void FVDemo::post_render_update()
{
    // Recreate the mesh.
    if (geom != nullptr) delete geom;
    geom = circle_mesh(mesh_N, random);

    // Solve the problem.
    auto solver = FVSolver(*geom);
    solver.set_dirichlet_boundary(dirichlet_boundary_function);
    solver.set_source([&](double x, double y)->double {
        float r = 0.3;
        if (x*x + y*y <= r*r) {
            return source_force;
        }
        return 0.0;
    });
    VertexAttachment<double> mesh_u = solver.solve();

    // Render the results.
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
    auto lifted_geom = SurfaceGeometry(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        double u_val = mesh_u[v];
        auto geom_pos = geom->position[v];
        lifted_geom.position[v] = Eigen::Vector3f(geom_pos.x(), u_val, geom_pos.z());
    }
    world->graphics.paint.wireframe(lifted_geom, mat4x4::identity(), 0.006f);
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
    auto demo = world.add<FVDemo>(demo_e, [](double x, double y)->double {
        return 1+0.4*sin(8*x);
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
