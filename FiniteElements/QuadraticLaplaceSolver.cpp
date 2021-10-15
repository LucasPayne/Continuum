#include "cg_sandbox.h"
#include "triangle_wrapper.h"
#include "CameraController.cpp"
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "mesh_generators.cpp"

FILE *error_metric_file = nullptr;



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


struct ManufacturedSolution {
    PlaneFunction dirichlet_boundary_function;
    PlaneFunction source_function;
};


class Solver {
public:
    Solver(SurfaceGeometry &_geom);

    // std::pair<VertexAttachment<double>, EdgeAttachment<double>> solve();
    void solve();

    void set_source(PlaneFunction func);
    void set_dirichlet_boundary(PlaneFunction func);

    
    // Solution (initializes to 0).
    VertexAttachment<double> vertex_u;
    EdgeAttachment<double> midpoint_u;

// private:
    SurfaceGeometry *geom; // Assumes y is 0. (this is bad, maybe should template the SurfaceGeometry class...)
                           // The geometry does not change.
    int num_interior_vertices;
    int num_boundary_vertices;
    int num_interior_midpoints;
    int num_boundary_midpoints;
    VertexAttachment<int> vertex_indices;
    EdgeAttachment<int> midpoint_indices;
    EdgeAttachment<Eigen::Vector3f> midpoints;

    PlaneFunction dirichlet_boundary_function;
    PlaneFunction source_function;
};

Solver::Solver(SurfaceGeometry &_geom) :
    // Solution.
    vertex_u(_geom.mesh),
    midpoint_u(_geom.mesh),

    geom{&_geom},
    vertex_indices(_geom.mesh),
    midpoint_indices(_geom.mesh),
    midpoints(_geom.mesh),
    source_function([](double,double)->double { return 0.0; }),
    dirichlet_boundary_function([](double,double)->double { return 0.0; })
{
    assert(geom->mesh.locked());

    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    num_boundary_midpoints = 0;
    num_interior_midpoints = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
    }
    for (auto edge : geom->mesh.edges()) {
        if (edge.on_boundary()) {
            midpoint_indices[edge] = -1;
            num_boundary_midpoints += 1;
        } else {
            midpoint_indices[edge] = num_interior_midpoints;
            num_interior_midpoints += 1;
        }
    }
    // Compute midpoints.
    for (auto edge : geom->mesh.edges()) {
        auto p = geom->position[edge.a().vertex()];
        auto pp = geom->position[edge.b().vertex()];
        midpoints[edge] = 0.5*p + 0.5*pp;
    }
    // printf("%d %d %d %d\n", num_boundary_vertices, num_interior_vertices,  
    //                          num_boundary_midpoints, num_interior_midpoints);
    // getchar();
}
void Solver::set_source(PlaneFunction func)
{
    source_function = func;
}
void Solver::set_dirichlet_boundary(PlaneFunction func)
{
    dirichlet_boundary_function = func;
}

// std::pair<VertexAttachment<double>, EdgeAttachment<double>> Solver::solve()
void Solver::solve()
{
    // The indices of interior vertices and interior midpoints are concatenated.
    //     note: When accessing a midpoint, the index must be shifted to start at num_interior_vertices.
    int N = num_interior_vertices + num_interior_midpoints;
    auto rhs = Eigen::VectorXd(N);
    for (int i = 0; i < N; i++) rhs[i] = 0.;
    std::vector<EigenTriplet> coefficients;
    auto add_entry = [&](int i, int j, double value) {
        // Add the value to entry (i,j) in the resulting matrix.
        // printf("%d %d %.5g\n", i,j,value);
        coefficients.push_back(EigenTriplet(i, j, value));
    };

    /*--------------------------------------------------------------------------------
        Compute the integrals for interior vertices.
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
            Edge vp_edge = he.edge(); // contains midpoint_vp
            Vertex vpp = he.next().tip();
            Edge vpp_edge = he.next().next().edge(); // contains midpoint_vpp
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];
            double C = 1.0/(4.0*geom->triangle_area(tri));
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            double val = 0.;
            // Diagonal term.
            val = C * 0.5 * K3.dot(K3);
            add_entry(v_index, v_index, val);

            // vp contribution.
            val = -(1./6.) * C * K1.dot(K3);
            if (vp.on_boundary()) {
                double bv = dirichlet_boundary_function(vp_pos.x(), vp_pos.z());
                rhs[v_index] -= bv * val;
            } else {
                int vp_index = vertex_indices[vp];
                add_entry(v_index, vp_index, val);
            }
            
            // vpp contribution.
            val = -(1./6.) * C * K2.dot(K3);
            if (vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(vpp_pos.x(), vpp_pos.z());
                rhs[v_index] -= bv * val;
            } else {
                int vpp_index = vertex_indices[vpp];
                add_entry(v_index, vpp_index, val);
            }

            // midpoint_vp contribution.
            val = (2./3.)*C*K1.dot(K3);
            int midpoint_vp_index = midpoint_indices[vp_edge];
            add_entry(v_index, num_interior_vertices + midpoint_vp_index, val);
            
            // midpoint_vpp contribution.
            val = (2./3.)*C*K2.dot(K3);
            int midpoint_vpp_index = midpoint_indices[vpp_edge];
            add_entry(v_index, num_interior_vertices + midpoint_vpp_index, val);

            he = he.twin().next();
        } while (he != start);
    }

    /*--------------------------------------------------------------------------------
        Compute the integrals for interior midpoints.
    --------------------------------------------------------------------------------*/
    for (auto edge : geom->mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};
        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            double C = 1.0/(4.0*geom->triangle_area(tri));
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            int v_index = vertex_indices[v];
            int vp_index = vertex_indices[vp];
            int vpp_index = vertex_indices[vpp];
            auto v_pos = geom->position[v];
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();
            auto midpoint_vp_index = midpoint_indices[midpoint_vp];
            auto midpoint_vpp_index = midpoint_indices[midpoint_vpp];
            auto midpoint_vp_pos = midpoints[midpoint_vp];
            auto midpoint_vpp_pos = midpoints[midpoint_vpp];

            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;

            // printf("self midpoint\n");
            double val = 0.;
            // 110, 110
            val = 4.*C/3. * (K3.dot(K3) - K1.dot(K2));
            // val = 4.*C/3. * (K1.dot(K1) + K1.dot(K2) + K2.dot(K2));
            add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_index, val);

            // printf("boundary midpoints\n");
            // 110, 011
            val = 4.*C/3. * (K1.dot(K3));
            if (midpoint_vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(midpoint_vpp_pos.x(), midpoint_vpp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vpp_index, val);
            }
            
            // 110, 101
            val = 4.*C/3. * (K2.dot(K3));
            if (midpoint_vp.on_boundary()) {
                double bv = dirichlet_boundary_function(midpoint_vp_pos.x(), midpoint_vp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, num_interior_vertices + midpoint_vp_index, val);
            }
            
            // printf("vertices\n");
            // 110, 200
            val = 2.*C/3. * (K1.dot(K2));
            if (vp.on_boundary()) {
                double bv = dirichlet_boundary_function(vp_pos.x(), vp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, vp_index, val);
            }
            // 110, 020
            if (vpp.on_boundary()) {
                double bv = dirichlet_boundary_function(vpp_pos.x(), vpp_pos.z());
                rhs[num_interior_vertices + midpoint_index] -= bv * val;
            } else {
                add_entry(num_interior_vertices + midpoint_index, vpp_index, val);
            }
        }
    }
    /*--------------------------------------------------------------------------------
        Source term integration.
    --------------------------------------------------------------------------------*/
    // Sample the source function at each vertex and midpoint.
    // This determines a piecewise quadratic interpolation.
    VertexAttachment<double> g_vertex(geom->mesh);
    EdgeAttachment<double> g_edge(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        auto v_pos = geom->position[v];
        g_vertex[v] = source_function(v_pos.x(), v_pos.z());
    }
    for (auto edge : geom->mesh.edges()) {
        auto midpoint_pos = midpoints[edge];
        g_edge[edge] = source_function(midpoint_pos.x(), midpoint_pos.z());
    }

    // Vertex trial integrals.
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) continue;
        int v_index = vertex_indices[v];
        auto v_pos = geom->position[v];

        double integral = 0.;

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            double tri_area = geom->triangle_area(tri);
            double J = tri_area;
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom->position[vp];
            auto vpp_pos = geom->position[vpp];

            auto opposite_edge = he.next().edge();

            // v
            integral += (1./60.) * J * g_vertex[v];
            // vp
            integral += (-1./360.) * J * g_vertex[vp];
            // vpp
            integral += (-1./360.) * J * g_vertex[vpp];
            // opposite_edge midpoint
            integral += (-1./90.) * J * g_edge[opposite_edge];

            he = he.next();
        } while (he != start);
        rhs[v_index] += integral;
    }

    // Midpoint trial integrals.
    for (auto edge : geom->mesh.edges()) {
        if (edge.on_boundary()) continue;
        int midpoint_index = midpoint_indices[edge];
        Halfedge hes[2] = {edge.a(), edge.b()};

        double integral = 0.;

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            double tri_area = geom->triangle_area(tri);
            double J = tri_area;
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            // Opposite triangle midpoints.
            auto midpoint_vp = he.next().next().edge();
            auto midpoint_vpp = he.next().edge();
            
            // edge
            integral += (4./45.) * J * g_edge[edge];
            // midpoint_vp
            integral += (2./45.) * J * g_edge[midpoint_vp];
            // midpoint_vpp
            integral += (2./45.) * J * g_edge[midpoint_vpp];
            // v
            integral += (-1./90.) * J * g_vertex[v];
        }
        rhs[num_interior_vertices + midpoint_index] += integral;
    }


    // Finalize the mass matrix.
    auto mass_matrix = SparseMatrix(N, N);
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    mass_matrix.makeCompressed();

    // std::cout << "mass matrix:\n" << Eigen::MatrixXd(mass_matrix) << "\n";
    // printf("rhs: ");
    // for (int i = 0; i < N; i++) printf("%.5g, ", rhs[i]);
    // printf("\n");
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
    // Reassociate each coefficient (or boundary value) with the corresponding vertex or edge of the mesh.
    --------------------------------------------------------------------------------*/

    int interior_vertex_index = 0;
    for (auto v : geom->mesh.vertices()) {
        if (v.on_boundary()) {
            auto p = geom->position[v];
            vertex_u[v] = dirichlet_boundary_function(p.x(), p.z());
        } else {
            vertex_u[v] = u[interior_vertex_index];
            interior_vertex_index += 1;
        }
    }

    int interior_midpoint_index = 0;
    for (auto edge : geom->mesh.edges()) {
        if (edge.on_boundary()) {
            auto p = midpoints[edge];
            midpoint_u[edge] = dirichlet_boundary_function(p.x(), p.z());
        } else {
            midpoint_u[edge] = u[num_interior_vertices + interior_midpoint_index];
            interior_midpoint_index += 1;
        }
    }

    // return {vertex_u, midpoint_u};
    // return std::pair<VertexAttachment<double>, EdgeAttachment<double>>(vertex_u, midpoint_u);
}




struct Demo : public IBehaviour {
    Demo();

    SurfaceGeometry *geom;

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
    GLShaderProgram quadratic_function_shader;


    int screenshot_blx;
    int screenshot_bly;
    int screenshot_trx;
    int screenshot_try;
    float f_screenshot_blx;
    float f_screenshot_bly;
    float f_screenshot_trx;
    float f_screenshot_try;

    // Test manufactured solutions.
    std::vector<ManufacturedSolution> manufactured_solutions;
    int sol_index;

    // Error metrics.
    GLuint sample_geometry_fbo;
    GLuint sample_geometry_texture;
    GLuint sample_geometry_depth_texture;
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
        if (e.key.code == KEY_V) {
            sol_index = (sol_index + 1) % manufactured_solutions.size();
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

Demo::Demo()
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

    quadratic_function_shader.add_shader(GLShader(VertexShader, SHADERS "quadratic_function/quadratic_function.vert"));
    quadratic_function_shader.add_shader(GLShader(TessControlShader, SHADERS "quadratic_function/quadratic_function.tcs"));
    quadratic_function_shader.add_shader(GLShader(TessEvaluationShader, SHADERS "quadratic_function/quadratic_function.tes"));
    quadratic_function_shader.add_shader(GLShader(FragmentShader, SHADERS "quadratic_function/quadratic_function.frag"));
    quadratic_function_shader.link();


    screenshot_blx = 0;
    screenshot_bly = 0;
    screenshot_trx = 10;
    screenshot_try = 10;
    f_screenshot_blx = 0;
    f_screenshot_bly = 0;
    f_screenshot_trx = 0.01;
    f_screenshot_try = 0.01;

    sol_index = 0;

    manufactured_solutions.push_back({
        [](double x, double y)->double {
            return exp(-2*(x*x + y*y));
        },
        [](double x, double y)->double {
            return (8 - 16*x*x - 16*y*y) * exp(-2*(x*x + y*y));
        }
    });
    manufactured_solutions.push_back({
        [](double x, double y)->double {
            return 0.5*cos(8*x)*y + exp(-x*x)*0.5*y*y+1.3;
        },
        [](double x, double y)->double {
            return -y*(2.0*x*x*y*exp(-x*x) - 1.0*y*exp(-x*x) - 32.0*cos(8*x)) - 1.0*exp(-x*x);
        }
    });
    manufactured_solutions.push_back({
        [](double x, double y)->double {
            return x*x - y*y + 1.3;
        },
        [](double x, double y)->double {
            return 0.;
        }
    });


    glGenFramebuffers(1, &sample_geometry_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, sample_geometry_fbo);

    glGenTextures(1, &sample_geometry_texture);
    glBindTexture(GL_TEXTURE_2D, sample_geometry_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 1024, 1024, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sample_geometry_texture, 0);

    glGenTextures(1, &sample_geometry_depth_texture);
    glBindTexture(GL_TEXTURE_2D, sample_geometry_depth_texture);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, 1024, 1024, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, sample_geometry_depth_texture, 0);

    GLenum framebuffer_status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if(framebuffer_status != GL_FRAMEBUFFER_COMPLETE) {
        fprintf(stderr, "Framebuffer incomplete.\n");
        exit(EXIT_FAILURE);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
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
    if (linear_mode) {
        // geom = circle_mesh(mesh_N*2, random);
        geom = square_mesh(mesh_N*2);
    } else {
        // geom = circle_mesh(mesh_N, random);
        geom = square_mesh(mesh_N);
    }
    // Compute Edge midpoints.
    auto midpoints = EdgeAttachment<Eigen::Vector3f>(geom->mesh);
    for (auto edge : geom->mesh.edges()) {
        auto pa = geom->position[edge.a().vertex()];
        auto pb = geom->position[edge.b().vertex()];
        midpoints[edge] = 0.5*pa + 0.5*pb;
    }

    // Solve the system.
    auto solver = Solver(*geom);
    solver.set_dirichlet_boundary(manufactured_solutions[sol_index].dirichlet_boundary_function);
    solver.set_source(manufactured_solutions[sol_index].source_function);
    solver.solve();

    // for (auto vertex : geom->mesh.vertices()) {
    //     printf("vertex: %.5g\n", solver.vertex_u[vertex]);
    // }
    // for (auto midpoint : geom->mesh.edges()) {
    //     printf("midpoint: %.5g\n", solver.midpoint_u[midpoint]);
    // }

    
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
            ps_lift.push_back(vec3(cos(theta), solver.dirichlet_boundary_function(cos(theta), sin(theta)), sin(theta)));
        }
        world->graphics.paint.chain(ps, 0.001, vec4(0,0,0,1));
        world->graphics.paint.chain(ps_lift, 0.001, vec4(0,0,0,1));
    }
    
    // Exact solution assumes that the Dirichlet boundary function is of a manufactured solution.
    auto exact_vertex_u = VertexAttachment<double>(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        auto p = geom->position[v];
        exact_vertex_u[v] = solver.dirichlet_boundary_function(p.x(), p.z());
    }
    auto exact_edge_u = EdgeAttachment<double>(geom->mesh);
    for (auto e : geom->mesh.edges()) {
        auto p = midpoints[e];
        exact_edge_u[e] = solver.dirichlet_boundary_function(p.x(), p.z());
    }
    auto &vertex_u = render_exact_solution ? exact_vertex_u : solver.vertex_u;
    auto &edge_u = render_exact_solution ? exact_edge_u : solver.midpoint_u;


    // // Compute the error.
    // // Compute square norms of each basis function.
    // VertexAttachment<double> vertex_square_norms(geom->mesh);
    // VertexAttachment<double> vertex_square_norms(geom->mesh);


    
#if 0
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
#endif

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
    glClear(GL_DEPTH_BUFFER_BIT); //...
    glEnable(GL_DEPTH_TEST);
    #if 1
    GLShaderProgram *quadratic_shaders[2] = {&quadratic_mesh_shader, &quadratic_function_shader};
    #else
    GLShaderProgram *quadratic_shaders[1] = {&quadratic_mesh_shader};
    #endif
    mat4x4 mvp_matrix = main_camera->view_projection_matrix();
    for (auto *shader : quadratic_shaders) {
        shader->bind();

            glUniformMatrix4fv(shader->uniform_location("mvp_matrix"), 1, GL_FALSE, (const GLfloat *) &mvp_matrix);
            glUniform1i(shader->uniform_location("linear_mode"), linear_mode ? 1 : 0);
            glPatchParameteri(GL_PATCH_VERTICES, 6);
            glDrawArrays(GL_PATCHES, 0, q_position_data.size());
        shader->unbind();
    }

    // // Draw into sample texture, for error calculation.
    quadratic_function_shader.bind();
    glUniformMatrix4fv(quadratic_function_shader.uniform_location("mvp_matrix"), 1, GL_FALSE, (const GLfloat *) &mvp_matrix);
    glUniform1i(quadratic_function_shader.uniform_location("linear_mode"), linear_mode ? 1 : 0);
    glBindFramebuffer(GL_FRAMEBUFFER, sample_geometry_fbo);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);
    glPatchParameteri(GL_PATCH_VERTICES, 6);
    glDrawArrays(GL_PATCHES, 0, q_position_data.size());
    glEnable(GL_DEPTH_TEST);
    glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);
    quadratic_function_shader.unbind();

    // Read in the result.
    glBindFramebuffer(GL_FRAMEBUFFER, sample_geometry_fbo);
    auto geometry_sample_pixels = std::vector<float>(4*1024*1024);
    glReadPixels(0,0,1024,1024, GL_RGBA, GL_FLOAT, &geometry_sample_pixels[0]);
    glBindFramebuffer(GL_FRAMEBUFFER, world->graphics.screen_buffer.id);

    float coeff = 1.f/(1024*1024); //----------should be for [-1,1]^2
    float total = 0.f;
    for (int i = 0; i < 1024; i++) {
        for (int j = 0; j < 1024; j++) {
            float val = geometry_sample_pixels[4*(1024*j + i)];
            total += val * coeff;
        }
    }
    // Compute the mesh parameter.
    float h = 0.f;
    float h_coeff = 1.f / geom->mesh.num_edges();
    for (auto edge : geom->mesh.edges()) {
        auto a = geom->position[edge.a().vertex()];
        auto b = geom->position[edge.b().vertex()];
        h += h_coeff * (a - b).norm();
    }

    // float error = sqrt(total);
    // fprintf(error_metric_file, "%d %.10f %.10f\n", mesh_N, h, error);
    // printf("%d %.10f %.10f\n", mesh_N, h, error);




    // Clean up.
    glDeleteVertexArrays(1, &q_vao);
    glDeleteBuffers(1, &q_positions);
    glDeleteBuffers(1, &q_values);
    //--------------------------------------------------------------------------------
    
    // Draw the lines of the wireframe.
    int line_N = 25;
    for (auto edge : geom->mesh.edges()) {
        auto p = geom->position[edge.a().vertex()];
        auto pp = geom->position[edge.b().vertex()];
        auto midpoint = midpoints[edge];
        auto ps = std::vector<vec3>(line_N);

        // Evaluate basis functions restricted to this edge.
        float val_a = solver.vertex_u[edge.a().vertex()];
        float val_b = solver.midpoint_u[edge];
        float val_c = solver.vertex_u[edge.b().vertex()];

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

    // Draw boundary condition
    // Square
    int num = 300;
    for (int t = -1; t <= 1; t += 2) {
        for (int axis = 0; axis <= 1; axis++) {
            auto bc = std::vector<vec3>();
            for (int i = 0; i <= num; i++) {
                double val = -1+i*2.f/num;
                vec2 v;
                if (axis == 0) v = vec2(t, val);
                else v = vec2(val, t);
                bc.push_back(vec3(v.x(), solver.dirichlet_boundary_function(v.x(), v.y()), v.y()));
            }
            world->graphics.paint.chain(bc, 0.01, vec4(0,0,0,1));
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
    auto demo = world.add<Demo>(demo_e);
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
    error_metric_file = fopen(DATA "error_metric_file.txt", "w+");

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

    fclose(error_metric_file);
}
