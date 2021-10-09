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

double alpha = 2;



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
    #if 1
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

        // Apparently...
        control_volume_area *= 2.f/3.f; //????????????????????????????????????????????????????????????

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

    bool render_exact_solution; // Toggleable option. This uses the same mesh as the approximate solution.



    // Plotting.
    GLShaderProgram sprite_shader;
    GLuint sprite_vao;

    GLShaderProgram sample_geometry_shader;

    GLuint sample_geometry_fbo;
    GLuint sample_geometry_texture;
    GLuint sample_geometry_depth_texture;

    void error_metrics(VertexAttachment<double> &mesh_u);
};


void FVDemo::keyboard_handler(KeyboardEvent e)
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
        if (e.key.code == KEY_T) {
            render_exact_solution = !render_exact_solution;
        }
    }
}
void FVDemo::mouse_handler(MouseEvent e)
{
}

FVDemo::FVDemo(PlaneFunction _dirichlet_boundary_function) :
    dirichlet_boundary_function{_dirichlet_boundary_function}
{
    mesh_N = 2;
    random = false;
    geom = nullptr;
    source_force = 0.0;
    render_exact_solution = false;

    // Geometry sampling.
    sample_geometry_shader.add_shader(GLShader(VertexShader, SHADERS "sample_geometry/sample_geometry.vert"));
    sample_geometry_shader.add_shader(GLShader(FragmentShader, SHADERS "sample_geometry/sample_geometry.frag"));
    sample_geometry_shader.link();

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
        

    // Plotting
    sprite_shader.add_shader(GLShader(VertexShader, SHADERS "plot/plot.vert"));
    sprite_shader.add_shader(GLShader(FragmentShader, SHADERS "plot/plot.frag"));
    sprite_shader.link();
    glGenVertexArrays(1, &sprite_vao);

    GLuint sprite_vbo;
    glGenBuffers(1, &sprite_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sprite_vbo);
    float sprite_uvs[2*4] = {
        0,0, 1,0, 1,1, 0,1
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4 * 2, sprite_uvs, GL_STATIC_DRAW);

    glGenVertexArrays(1, &sprite_vao);
    glBindVertexArray(sprite_vao);
    glBindBuffer(GL_ARRAY_BUFFER, sprite_vbo);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);


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


void FVDemo::error_metrics(VertexAttachment<double> &mesh_u)
{
    // Quality metric.
    auto sample_positions = std::vector<vec2>();
    auto sample_values = std::vector<float>();
    for (auto tri : geom->mesh.faces()) {
        auto start = tri.halfedge();
        auto he = start;
        do {
            auto v = he.vertex();
            auto p = geom->position[v];
            sample_positions.push_back(vec2(p.x(), p.z()));
            // Reconstruct from hat functions. This hat function has a coefficient in the interior, and a boundary value on the boundary.
            if (v.on_boundary()) {
                sample_values.push_back(dirichlet_boundary_function(p.x(), p.z()));
            } else {
                sample_values.push_back(mesh_u[v]);
            }
            he = he.next();
        } while (he != start);
    }
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLuint vbo_positions;
    glGenBuffers(1, &vbo_positions);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_positions);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * sample_positions.size(), &sample_positions[0], GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

    GLuint vbo_values;
    glGenBuffers(1, &vbo_values);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_values);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * sample_values.size(), &sample_values[0], GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(1);

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport); // save the viewport to restore after rendering sample geometry.
    glDisable(GL_SCISSOR_TEST);
    glViewport(0,0,1024,1024);

    sample_geometry_shader.bind();
    glBindFramebuffer(GL_FRAMEBUFFER, sample_geometry_fbo);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_TRIANGLES, 0, sample_positions.size());
    glEnable(GL_DEPTH_TEST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    sample_geometry_shader.unbind();

    glEnable(GL_SCISSOR_TEST);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo_positions);
    glDeleteBuffers(1, &vbo_values);

    // Read in the result.
    glBindFramebuffer(GL_FRAMEBUFFER, sample_geometry_fbo);
    auto geometry_sample_pixels = std::vector<float>(4*1024*1024);
    glReadPixels(0,0,1024,1024, GL_RGBA, GL_FLOAT, &geometry_sample_pixels[0]);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    float coeff = 1.f/(1024*1024);
    float total = 0.f;
    for (int i = 0; i < 1024; i++) {
        for (int j = 0; j < 1024; j++) {
            float val = geometry_sample_pixels[4*(1024*j + i)];
            total += val * coeff;
        }
    }
    float h = 0.f;
    float h_coeff = 1.f / geom->mesh.num_edges();
    for (auto edge : geom->mesh.edges()) {
        auto a = geom->position[edge.a().vertex()];
        auto b = geom->position[edge.b().vertex()];
        h += h_coeff * (a - b).norm();
    }

    float error = sqrt(total);
    fprintf(error_metric_file, "%d %.10f %.10f\n", mesh_N, h, error);
    printf("%d %.10f %.10f\n", mesh_N, h, error);
}

void FVDemo::post_render_update()
{
    // Recreate the mesh.
    if (geom != nullptr) delete geom;
    srand(230192301);
    geom = circle_mesh(mesh_N, random);
    // geom = square_mesh(mesh_N);

    // Solve the problem.
    auto solver = FVSolver(*geom);
    solver.set_dirichlet_boundary(dirichlet_boundary_function);
    solver.set_source([&](double x, double y)->double {
        // return -(4*pow(alpha, 2)*pow(x, 2)*exp(-alpha*pow(x, 2))*exp(-alpha*pow(y, 2)) + 4*pow(alpha, 2)*pow(y, 2)*exp(-alpha*pow(x, 2))*exp(-alpha*pow(y, 2)) - 4*alpha*exp(-alpha*pow(x, 2))*exp(-alpha*pow(y, 2)));
        
        return 0.0;


        // float r = 0.3;
        // if (x*x + y*y <= r*r) {
        //     return 10;
        // }
        // return 0.0;
    });
    VertexAttachment<double> mesh_u = solver.solve();

    // Render the results.
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
    auto lifted_geom = SurfaceGeometry(geom->mesh);
    for (auto v : geom->mesh.vertices()) {
        double u_val = mesh_u[v];
        auto geom_pos = geom->position[v];
        lifted_geom.position[v] = Eigen::Vector3f(geom_pos.x(), render_exact_solution ? dirichlet_boundary_function(geom_pos.x(), geom_pos.z()) : u_val, geom_pos.z());
    }
    world->graphics.paint.wireframe(lifted_geom, mat4x4::identity(), 0.006f);

    // Plot the solution.
    sprite_shader.bind();
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glUniform2f(sprite_shader.uniform_location("bottom_left"), 0,0);
    float wid = 0.4;
    glUniform2f(sprite_shader.uniform_location("top_right"), wid, wid/0.566);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, sample_geometry_texture);
    glUniform1i(sprite_shader.uniform_location("tex"), 0);

    glBindVertexArray(sprite_vao);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    sprite_shader.unbind();

    if (1) {
        // Square
        // Draw the exact boundary condition.
        #if 0
        int num = 300;
        for (int t = -1; t <= 1; t += 2) {
            for (int axis = 0; axis <= 1; axis++) {
	        auto bc = std::vector<vec3>();
                for (int i = 0; i <= num; i++) {
                    double val = -1+i*2.f/num;
                    vec2 v;
                    if (axis == 0) v = vec2(t, val);
                    else v = vec2(val, t);
                    bc.push_back(vec3(v.x(), dirichlet_boundary_function(v.x(), v.y()), v.y()));
                }
	        world->graphics.paint.chain(bc, 0.01, vec4(0,0,0,1));
            }
        }
        #else
        // Circle
        int num = 300;
        auto boundary_condition_loop = std::vector<vec3>(num+1);
        for (int i = 0; i <= num; i++) {
            float theta = 2*i*M_PI/num;
            float c = cos(theta);
            float s = sin(theta);
            boundary_condition_loop[i] = vec3(c, dirichlet_boundary_function(c, s), s);
        }
        world->graphics.paint.chain(boundary_condition_loop, 0.01, vec4(0,0,0,1));
        #endif
    }


    error_metrics(mesh_u);
    // mesh_N += 1;
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
        // return 1+0.4*sin(8*x);
        // return cos(x) + y*y;
        // return exp(-alpha*(x*x + y*y));
        return x*x - y*y + 1.13;
        
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
