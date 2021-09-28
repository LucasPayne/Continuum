#include "cg_sandbox.h"
#include "triangle_wrapper.h"
#include "CameraController.cpp"


class LaplaceSolver {
    using PlaneFunction = std::function<double>(double, double);
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using EigenTriplet = Eigen::Triplet<double>;
public:
    LaplaceSolver(SurfaceGeometry &_geom);

    void set_dirichlet_boundary(PlaneFunction func);

    void solve();

private:
    SurfaceGeometry &geom; // Assumes y is 0. The geometry does not change.
    PlaneFunction dirichlet_boundary_function;
    int num_interior_vertices;
    int num_boundary_vertices;
    VertexAttachment<int> vertex_indices;
};

LaplaceSolver::LaplaceSolver(SurfaceGeometry &_geom) :
    geom{_geom}, vertex_indices(_geom.mesh), dirichlet_boundary_function([](double,double)-double> { return 0.0; })
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    {
        int full_index = 0;
        for (auto v : geom.mesh.vertices()) {
            if (v.on_boundary()) {
                num_boundary_vertices += 1;
                vertex_indices[full_index] = -1;
            } else {
                num_interior_vertices += 1;
                vertex_indices[full_index] = num_interior_vertices;
            }
            full_index += 1;
        }
    }
}
void LaplaceSolver::set_dirichlet_boundary(PlaneFunction func)
{
    dirichlet_boundary_function = func;
}


void LaplaceSolver::solve() const
{
    std::vector<EigenTriplet> coefficients;
    auto rhs = Eigen::VectorXd<double>(num_interior_vertices);

    for (auto tri : geom.mesh.faces()) {
        Vertex verts[3];
        int vertex_index = 0;
        auto he = tri.halfedge();
        do {
            he = he.next();
            verts[vertex_index] = he.vertex()
        } while (++vertex_index < 3);

        auto A = geom.positions[verts[0]];
        auto B = geom.positions[verts[1]];
        auto C = geom.positions[verts[2]];
        double triangle_area = 0.5*(B-A).cross(C-A).norm();
        
        Eigen::Matrix<double, 3,3> local_stiffness_matrix = {
            {1, -0.5, -0.5},
            {-0.5, 0.5, 0},
            {-0.5, 0, 0.5}
        };
        local_stiffness_matrix /= 2*triangle_area;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (verts[j].on_boundary()) {
                    rhs[i] -= dirichlet_boundary_function
                } else {
                    coefficients.pushback(EigenTriplet(vertex_indices[verts[i]], vertex_indices[verts[j]], local_stiffness_matrix[i,j]));
                }
            }
        }
    }

    auto mass_matrix = SparseMatrix(num_interior_vertices, num_interior_vertices);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
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
    std::vector<vec2> interior_points;
    std::vector<vec2> boundary_points;
    std::vector<vec2> all_points;

    SurfaceMesh mesh;
    SurfaceGeometry geom;
};


App::App(World &_world) : world{_world}, geom{mesh}
{
    auto cameraman = world.entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    auto t = cameraman.add<Transform>(0,0,2);
    main_camera = camera;
    auto controller = world.add<CameraController>(cameraman);
    t->position = vec3(0,1,3);
    controller->angle = -2;
    
    
    
    // Create a roughly regular point cloud for the circle, and points sampled on the boundary.
    int N = 30;
    for (int i = 0; i <= N; i++) {
        float x = -1+(i*2.f)/N;
        for (int j = 0; j <= N; j++) {
            float y = -1+(j*2.f)/N;
            if (x*x + y*y <= 1-1.1f/N) interior_points.push_back(vec2(x,y));
        }
    }
    int num_angle_intervals = 2*N;
    for (int i = 0; i < num_angle_intervals; i++) {
        float theta = (i*2.f*M_PI)/num_angle_intervals;
        boundary_points.push_back(vec2(cos(theta), sin(theta)));
    }
    // for (auto *ps : {&interior_points, &boundary_points}) {
    //     for (int i = 0; i < ps->size(); i++) {
    //         (*ps)[i] = 0.2f*(*ps)[i] + vec2(0.5,0.5);
    //     }
    // }
    all_points.insert(all_points.end(), interior_points.begin(), interior_points.end());
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());

    // Write to a .node file.
    FILE *node_file = fopen(DATA "circle_cloud.node", "w+");
    fprintf(node_file, "%zu 2 0 0\n", all_points.size()); //num_vertices, num_dimensions=2, num_attributes, num_boundary_markers=0,1
    for (int i = 0; i < all_points.size(); i++) {
        fprintf(node_file, "%d %.6f %.6f\n", i, all_points[i].x(), all_points[i].y());
    }
    fclose(node_file);

    // Triangulate this point cloud.
    std::string triswitches = "zV";
    double *p_mem = (double *) malloc(2*sizeof(double)*all_points.size());
    for (int i = 0; i < all_points.size(); i++) {
        p_mem[2*i] = all_points[i].x();
        p_mem[2*i+1] = all_points[i].y();
    }
    struct triangulateio in = {0};
    in.pointlist = p_mem;
    in.numberofpoints = all_points.size();
    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;

    struct triangulateio out = {0};
    out.pointlist = nullptr;

    in.numberofpointattributes = 0;
    triangulate(&triswitches[0], &in, &out, nullptr);
    

    auto vertices = std::vector<Vertex>(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++) {
        auto v = geom.mesh.add_vertex();
        geom.position[v] = Eigen::Vector3f(out.pointlist[2*i+0], 0, out.pointlist[2*i+1]);
        vertices[i] = v;
    }
    for (int i = 0; i < out.numberoftriangles; i++) {
        int a = out.trianglelist[3*i+0];
        int b = out.trianglelist[3*i+1];
        int c = out.trianglelist[3*i+2];
        geom.mesh.add_triangle(vertices[a], vertices[b], vertices[c]);
    }
    geom.mesh.lock();
}

void App::close()
{
}

void App::loop()
{
    world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.001);

    auto solver = LaplaceSolver(geom);
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
