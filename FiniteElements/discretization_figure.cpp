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



vec2 curve(float t)
{
    // return vec2(0.5+0.2*cos(t)+0.1*t*(2*M_PI-t), 0.5+0.4*sin(t));
    return vec2(0.5+0.2*cos(t) + 0.1*sin(t), 0.5+0.2*sin(t));
}

enum CenterModes {
    CenterModeBarycentric,
    CenterModeVoronoi,
    CenterModeMixedVoronoi,
    NUM_CENTER_MODES
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

    std::vector<Eigen::Vector3f> ps;
    std::vector<vec2> ps_draw;
    SurfaceMesh mesh;
    SurfaceGeometry geom;

    int center_mode;
};



App::App(World &_world) :
    world{_world},
    geom(mesh)
{
    center_mode = CenterModeBarycentric;
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

    //--------------------------------------------------------------------------------
    // Figure
    //--------------------------------------------------------------------------------
    // Compute curve boundary points.
    int curve_N = 360;
    ps_draw = std::vector<vec2>(curve_N+1);
    for (int i = 0; i <= curve_N; i++) {
        float t = i*2*M_PI/curve_N;
        vec2 c = curve(t);
        ps_draw[i] = c;
        // aspect ratio...
        ps_draw[i].y() *= 1/0.566;
        ps_draw[i] *= 0.5;
    }

    // int boundary_N = 18;
    int boundary_N = 10;
    auto boundary_ps = std::vector<Eigen::Vector3f>(boundary_N);
    for (int i = 0; i < boundary_N; i++) {
        float t = i*2*M_PI/boundary_N;
        vec2 c = curve(t);
        boundary_ps[i] = Eigen::Vector3f(c.x(), 0, c.y());
    }

    // Add interior points.
    // int interior_N = 5;
    int interior_N = 2;
    float from = 0.43f;
    float to = 0.57f;
    for (int i = 0; i <= interior_N; i++) {
        float x = from + i*(to-from)/interior_N;
        for (int j = 0; j <= interior_N; j++) {
            float y = from + j*(to-from)/interior_N;
            auto p = Eigen::Vector3f(x,0,y);
            bool exterior = false;
            for (int k = 0; k < boundary_ps.size()-1; k++) {
                auto A = boundary_ps[k];
                auto B = boundary_ps[k+1];
                float det = (A.x()-p.x())*(B.z()-p.z()) - (B.x()-p.x())*(A.z()-p.z());
                if (det < 0) {
                    exterior = true;
                    break;
                }
            }
            if (!exterior) ps.push_back(p);
        }
    }
    for (auto bp : boundary_ps) ps.push_back(bp);
    
    // Triangulate.
    std::string triswitches = "zV";
    double *p_mem = (double *) malloc(2*sizeof(double)*ps.size());
    for (int i = 0; i < ps.size(); i++) {
        p_mem[2*i] = ps[i].x();
        p_mem[2*i+1] = ps[i].z();
    }
    struct triangulateio in = {0};
    in.pointlist = p_mem;
    in.numberofpoints = ps.size();
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
        // Move some vertices to visualize hat functions.
        if (i == 7) geom.position[v].y() += 0.15;
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
    world.graphics.paint.chain_2D(ps_draw, 4, vec4(0,0,0,1));
    // world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.003);
    world.graphics.paint.wireframe(geom, mat4x4::identity(), 0.003);
    

    auto get_center = [&](Face tri)->Eigen::Vector3f {
        if (center_mode == CenterModeBarycentric) {
            return geom.barycenter(tri);
        } else if (center_mode == CenterModeVoronoi || center_mode == CenterModeMixedVoronoi) {
            Eigen::Vector3f positions[3];
            int vertex_index = 0;
            auto he = tri.halfedge();
            do {
                he = he.next();
                positions[vertex_index] = geom.position[he.vertex()];
            } while (++vertex_index < 3);

            // https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
            Eigen::Vector3f a = positions[0];
            Eigen::Vector3f b = positions[1];
            Eigen::Vector3f c = positions[2];
            Eigen::Vector3f ac = c - a ;
            Eigen::Vector3f ab = b - a ;
            Eigen::Vector3f abXac = ab.cross(ac);

            // this is the vector from a TO the circumsphere center
            Eigen::Vector3f toCircumsphereCenter = (abXac.cross( ab )*ac.dot(ac) + ac.cross( abXac )*ab.dot(ab)) / (2.f*abXac.dot(abXac)) ;
            float circumsphereRadius = toCircumsphereCenter.norm();
            // The 3 space coords of the circumsphere center then:
            Eigen::Vector3f ccs = a  +  toCircumsphereCenter ; // now this is the actual 3space location

            // if (center_mode == CenterModeMixedVoronoi) {
            //     if (ccs - a
            // }

            return ccs;
        }
        assert(0);
    };

    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        // For each adjacent triangle.
        auto start = v.halfedge(); //todo: provide this circulator in mesh_processing.
        auto he = start;
        auto cell_boundary = std::vector<vec3>();

        vec3 vp = *((vec3 *) &geom.position[v][0]);
        //world.graphics.paint.sphere(vp, 0.02, vec4(1,0,0,1));
        do {
            auto tri = he.face();
            auto barycenter = get_center(tri);

            auto next_tri = he.twin().next().face();
            auto next_barycenter = get_center(next_tri);
            vec3 a = *((vec3 *) &barycenter[0]);
            vec3 b = *((vec3 *) &next_barycenter[0]);
            vec3 tip = *((vec3 *) &geom.position[he.tip()][0]);
            // Intersect a<->b with vp<->tip.
            auto get_line = [](vec3 A, vec3 B)->vec3 {
                return vec3(
                    A.z() - B.z(),
                    B.x() - A.x(),
                    A.x()*(B.z() - A.z()) - A.z()*(B.x() - A.x())
                );
            };
            vec3 intersect = vec3::cross(get_line(a,b), get_line(vp,tip));
            intersect /= intersect.z();
            float dist = sqrt((intersect.x()-vp.x())*(intersect.x()-vp.x())
                                +(intersect.y()-vp.z())*(intersect.y()-vp.z()));
            float mid_length = sqrt((tip.x()-vp.x())*(tip.x()-vp.x())
                                    +(tip.z()-vp.z())*(tip.z()-vp.z()));
            float t = dist / mid_length;
            vec3 midpoint = vp*(1-t) + tip*t;

            cell_boundary.push_back(midpoint);
            cell_boundary.push_back(b);
            
            he = he.twin().next();
        } while (he != start);
        if (cell_boundary.size() > 1) {
            cell_boundary.push_back(cell_boundary[0]);
            cell_boundary.push_back(cell_boundary[1]);
        }

        // push it up a little bit
        for (int i = 0; i < cell_boundary.size(); i++) {
            cell_boundary[i].y() += 0.001;
        }
        // world.graphics.paint.chain(cell_boundary, 15, vec4(1,0.6,0.6,1));
        world.graphics.paint.chain(cell_boundary, 0.005, vec4(1,0,0,1));
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
            center_mode = (center_mode + 1) % NUM_CENTER_MODES;
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
