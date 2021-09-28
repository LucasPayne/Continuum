#include "cg_sandbox.h"
#include "triangle_wrapper.h"
#include "CameraController.cpp"
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>



class LaplaceSolver {
    using PlaneFunction = std::function<double(double, double)>;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using EigenTriplet = Eigen::Triplet<double>;
public:
    LaplaceSolver(SurfaceGeometry &_geom);

    void set_dirichlet_boundary(PlaneFunction func);

    VertexAttachment<double> solve();

private:
    SurfaceGeometry &geom; // Assumes y is 0. (this is bad, maybe should template the SurfaceGeometry class...)
        //The geometry does not change.
    int num_interior_vertices;
    int num_boundary_vertices;
    VertexAttachment<int> vertex_indices;
    PlaneFunction dirichlet_boundary_function;
};

LaplaceSolver::LaplaceSolver(SurfaceGeometry &_geom) :
    geom{_geom}, vertex_indices(_geom.mesh), dirichlet_boundary_function([](double,double)->double { return 0.0; })
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
    }
}
void LaplaceSolver::set_dirichlet_boundary(PlaneFunction func)
{
    dirichlet_boundary_function = func;
}


VertexAttachment<double> LaplaceSolver::solve()
{
    std::vector<EigenTriplet> coefficients;
    auto rhs = Eigen::VectorXd(num_interior_vertices);
    for (int i = 0; i < num_interior_vertices; i++) rhs[i] = 0; //---zero initialize?

    for (auto tri : geom.mesh.faces()) {
        Vertex verts[3];
        int vertex_index = 0;
        auto he = tri.halfedge();
        do {
            he = he.next();
            verts[vertex_index] = he.vertex();
        } while (++vertex_index < 3);

        auto A = geom.position[verts[0]];
        auto B = geom.position[verts[1]];
        auto C = geom.position[verts[2]];
        double triangle_area = 0.5*(B-A).cross(C-A).norm();

        Eigen::Matrix<double, 2,2> jacobian_matrix;
        jacobian_matrix <<
            B.x() - A.x(),  C.x() - A.x(),
            B.z() - A.z(),  C.z() - A.z();
        // auto jacobian_matrix_inverse = jacobian_matrix.inverse();
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
                local_stiffness_matrix(i,j) = grad_i.dot(grad_j);
            }
        }

        
        // Eigen::Matrix<double, 3,3> local_stiffness_matrix;
        // local_stiffness_matrix << 
        //     1, -0.5, -0.5,
        //     -0.5, 0.5, 0,
        //     -0.5, 0, 0.5;
        // local_stiffness_matrix *= 2*triangle_area;

        // i: Row, trial function.
        // j: Column, test function.

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                bool b1 = verts[i].on_boundary();
                bool b2 = verts[j].on_boundary();
                if (b1) {
                    continue; // The trial function can't be centered on the boundary.
                } else if (b2) {
                    auto p = geom.position[verts[j]];
                    double boundary_val = dirichlet_boundary_function(p.x(), p.z());
                    rhs[vertex_indices[verts[i]]] -= boundary_val*local_stiffness_matrix(i,j);
                } else {
                    coefficients.push_back(EigenTriplet(vertex_indices[verts[i]], vertex_indices[verts[j]], local_stiffness_matrix(i,j)));
                }
            }
        }

        // std::cout << local_stiffness_matrix << "\n";
        // for (int i = 0; i < 3; i++) std::cout << (verts[i].on_boundary() ? "b " : "_ ");
        // std::cout << "\n";
        // std::cout << rhs << "\n";
        // getchar();
    }
    auto mass_matrix = SparseMatrix(num_interior_vertices, num_interior_vertices);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();
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

    auto mesh_u = VertexAttachment<double>(geom.mesh);
    // Reassociate each coefficient (or boundary value) with the corresponding vertex of the mesh.
    int interior_vertex_index = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            auto p = geom.position[v];
            mesh_u[v] = dirichlet_boundary_function(p.x(), p.z());
        } else {
            mesh_u[v] = u[interior_vertex_index];
            interior_vertex_index += 1;
        } 
    }
    // std::cout << u << "\n";
    return mesh_u;
}


struct SolvingMesh {
    SolvingMesh(int N);

    std::vector<vec2> interior_points;
    std::vector<vec2> boundary_points;
    std::vector<vec2> all_points;

    SurfaceMesh mesh;
    SurfaceGeometry geom;
};
SolvingMesh::SolvingMesh(int N) : geom(mesh)
{
    // Create a roughly regular point cloud for the circle, and points sampled on the boundary.
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
    all_points.insert(all_points.end(), interior_points.begin(), interior_points.end());
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());

    // Write to a .node file.
    // FILE *node_file = fopen(DATA "circle_cloud.node", "w+");
    // fprintf(node_file, "%zu 2 0 0\n", all_points.size()); //num_vertices, num_dimensions=2, num_attributes, num_boundary_markers=0,1
    // for (int i = 0; i < all_points.size(); i++) {
    //     fprintf(node_file, "%d %.6f %.6f\n", i, all_points[i].x(), all_points[i].y());
    // }
    // fclose(node_file);

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
    free(p_mem);

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

    int mesh_N;
    SolvingMesh *solving_mesh;
};



App::App(World &_world) : world{_world}
{
    auto cameraman = world.entities.add();
    auto camera = cameraman.add<Camera>(0.1, 300, 0.1, 0.566);
    auto t = cameraman.add<Transform>(0,0,2);
    main_camera = camera;
    auto controller = world.add<CameraController>(cameraman);
    t->position = vec3(0,1,3);
    controller->angle = -2;
    
    mesh_N = 5;
    solving_mesh = new SolvingMesh(mesh_N);
}

void App::close()
{
}

void App::loop()
{
    world.graphics.paint.wireframe(solving_mesh->geom, mat4x4::identity(), 0.001);

    auto solver = LaplaceSolver(solving_mesh->geom);
    solver.set_dirichlet_boundary([](double x, double y)->double {
        // return 1+0.4*sin(8*x);
        return x < 0 ? 1.5 : 0.5;
        // return 1;
    });
    VertexAttachment<double> mesh_u = solver.solve();
    auto lifted_geom = SurfaceGeometry(solving_mesh->geom.mesh);
    for (auto v : solving_mesh->geom.mesh.vertices()) {
        double u_val = mesh_u[v];
        auto geom_pos = solving_mesh->geom.position[v];
        lifted_geom.position[v] = Eigen::Vector3f(geom_pos.x(), u_val, geom_pos.z());
    }
    world.graphics.paint.wireframe(lifted_geom, mat4x4::identity(), 0.001);

    for (auto v : solving_mesh->geom.mesh.vertices()) {
        auto p = solving_mesh->geom.position[v];
        vec3 pp = *((vec3 *) &p[0]);
        if (v.on_boundary()) world.graphics.paint.sphere(pp, 0.01, vec4(1,0,0,1));
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
        if (e.key.code == KEY_P) {
            mesh_N += 5;
            delete solving_mesh;
            solving_mesh = new SolvingMesh(mesh_N);
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
