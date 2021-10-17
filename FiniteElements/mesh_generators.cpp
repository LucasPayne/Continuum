
// Create a mesh.
SurfaceGeometry *circle_mesh(int N, bool random)
{
    SurfaceMesh *mesh = new SurfaceMesh();
    SurfaceGeometry *geom = new SurfaceGeometry(*mesh);

    std::vector<vec2> interior_points;
    std::vector<vec2> boundary_points;
    std::vector<vec2> all_points;


    // Create a roughly regular point cloud for the circle, and points sampled on the boundary.
    if (random) {
        for (int i = 0; i <= N*N; i++) {
            vec2 rand_vec = vec2::random(-1,1);
            if (vec2::dot(rand_vec, rand_vec) <= 1-1.1/N) {
                interior_points.push_back(rand_vec);
            }
        }
    } else {
        for (int i = 0; i <= N; i++) {
            float x = -1+(i*2.f)/N;
            for (int j = 0; j <= N; j++) {
                float y = -1+(j*2.f)/N;
                if (x*x + y*y <= 1-1.1f/N) {
                    interior_points.push_back(vec2(x,y));
                }
            }
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
    // std::string triswitches = "zV";
    std::string triswitches = "z";
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
        auto v = geom->mesh.add_vertex();
        geom->position[v] = Eigen::Vector3f(out.pointlist[2*i+0], 0, out.pointlist[2*i+1]);
        vertices[i] = v;
    }
    for (int i = 0; i < out.numberoftriangles; i++) {
        int a = out.trianglelist[3*i+0];
        int b = out.trianglelist[3*i+1];
        int c = out.trianglelist[3*i+2];
        geom->mesh.add_triangle(vertices[a], vertices[b], vertices[c]);
    }
    geom->mesh.lock();
    return geom;
}

// [-1,1]^2
SurfaceGeometry *square_mesh(int N)
{
    SurfaceMesh *mesh = new SurfaceMesh();
    SurfaceGeometry *geom = new SurfaceGeometry(*mesh);
    
    auto vertices = std::vector<Vertex>((N+1)*(N+1));
    // Create vertices.
    for (int i = 0; i <= N; i++) {
        double x = -1 + i*2.f/N;
        for (int j = 0; j <= N; j++) {
            double y = -1 + j*2.f/N;
            auto v = mesh->add_vertex();
            geom->position[v] = Eigen::Vector3f(x, 0, y);
            vertices[i*(N+1) + j] = v;
        }
    }
    // Create triangles.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Vertex a = vertices[i*(N+1) + j];
            Vertex b = vertices[(i+1)*(N+1) + j];
            Vertex c = vertices[(i+1)*(N+1) + j+1];
            Vertex d = vertices[i*(N+1) + j+1];
            mesh->add_triangle(a,b,d);
            mesh->add_triangle(b,c,d);
        }
    }
    mesh->lock();
    return geom;
}

SurfaceGeometry *square_mesh_with_square_hole(int N)
{
    N = std::max(N, 10);
    SurfaceMesh *mesh = new SurfaceMesh();
    SurfaceGeometry *geom = new SurfaceGeometry(*mesh);

    auto vertices = std::vector<Vertex>((N+1)*(N+1));
    // Create vertices.
    for (int i = 0; i <= N; i++) {
        double x = -1 + i*2.f/N;
        for (int j = 0; j <= N; j++) {
            double y = -1 + j*2.f/N;
            auto v = mesh->add_vertex();
            geom->position[v] = Eigen::Vector3f(x, 0, y);
            vertices[i*(N+1) + j] = v;
        }
    }
    // Create triangles.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Vertex a = vertices[i*(N+1) + j];
            Vertex b = vertices[(i+1)*(N+1) + j];
            Vertex c = vertices[(i+1)*(N+1) + j+1];
            Vertex d = vertices[i*(N+1) + j+1];
            if (i != N/2 || j != N/2) {
                mesh->add_triangle(a,b,d);
                mesh->add_triangle(b,c,d);
            }
        }
    }
    mesh->lock();
    return geom;
}

// Return with a data structure which can be used to access the square mesh with rectangular indexing.
class SquareMesh
{
public:
    SquareMesh(int _N) :
        N{_N}
    {
        SurfaceMesh *mesh = new SurfaceMesh();
        geom = new SurfaceGeometry(*mesh);
        
        vertices = std::vector<Vertex>((N+1)*(N+1));
        // Create vertices.
        for (int i = 0; i <= N; i++) {
            double x = -1 + i*2.f/N;
            for (int j = 0; j <= N; j++) {
                double y = -1 + j*2.f/N;
                auto v = mesh->add_vertex();
                geom->position[v] = Eigen::Vector3f(x, 0, y);
                vertices[i*(N+1) + j] = v;
            }
        }
        // Create triangles.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Vertex a = vertices[i*(N+1) + j];
                Vertex b = vertices[(i+1)*(N+1) + j];
                Vertex c = vertices[(i+1)*(N+1) + j+1];
                Vertex d = vertices[i*(N+1) + j+1];
                mesh->add_triangle(a,b,d);
                mesh->add_triangle(b,c,d);
            }
        }
        mesh->lock();
    }

    inline Vertex vertex(int i, int j) {
        return vertices[i*(N+1) + j];
    }

    SurfaceGeometry *geom;
private:
    int N;
    std::vector<Vertex> vertices;
};



SurfaceGeometry *square_minus_circle(float r, float theta0, float a, float b, int square_N, bool misc_curve=false, vec2 obstruction_position=vec2(0,0),
                                     bool many_sample_curve=false,
                                     float rect_x_scale = 1.f,
                                     float rect_y_scale = 1.f
) {
    SurfaceMesh *mesh = new SurfaceMesh();
    SurfaceGeometry *geom = new SurfaceGeometry(*mesh);

    #if 0
    int N = 10;
    std::vector<vec2> interior_points;
    std::vector<vec2> boundary_points;
    std::vector<vec2> all_points;

    // Create a roughly regular point cloud for the circle, and points sampled on the boundary.
    if (random) {
        for (int i = 0; i <= N*N; i++) {
            vec2 rand_vec = vec2::random(-1,1);
            if (vec2::dot(rand_vec, rand_vec) <= 1-1.1/N) {
                interior_points.push_back(rand_vec);
            }
        }
    } else {
        for (int i = 0; i <= N; i++) {
            float x = -1+(i*2.f)/N;
            for (int j = 0; j <= N; j++) {
                float y = -1+(j*2.f)/N;
                if (x*x + y*y <= 1-1.1f/N) {
                    interior_points.push_back(vec2(x,y));
                }
            }
        }
    }
    int num_angle_intervals = 2*N;
    for (int i = 0; i < num_angle_intervals; i++) {
        float theta = (i*2.f*M_PI)/num_angle_intervals;
        boundary_points.push_back(vec2(cos(theta), sin(theta)));
    }
    all_points.insert(all_points.end(), interior_points.begin(), interior_points.end());
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());
    #else
    std::vector<vec2> all_points;

    int circle_N = many_sample_curve ? 150 : 4*square_N;
    std::vector<int> segment_pieces;

    // float theta0 = 1.2;
    // float a = 0.7;
    // float b = 1.34;
    float cos_theta0 = cos(theta0);
    float sin_theta0 = sin(theta0);

    for (int i = 0; i < circle_N; i++) {
        float theta = (i*2.f*M_PI)/circle_N;
        // all_points.push_back(r*vec2(cos(theta), sin(theta)));
        if (misc_curve) {
            all_points.push_back(obstruction_position + r*vec2(a*cos_theta0*cos(theta) - b*sin_theta0*sin(theta), a*sin_theta0*cos(theta) + b*cos_theta0*sin(theta)));
        } else {
            double c = 1 + 0.7*pow(sin(theta*2), 4);
            all_points.push_back(obstruction_position + c*r*vec2(a*cos_theta0*cos(theta) - b*sin_theta0*sin(theta), a*sin_theta0*cos(theta) + b*cos_theta0*sin(theta)));
        }
        segment_pieces.push_back(i);
        segment_pieces.push_back((i+1)%circle_N);
    }

    int square_N_x = 4*1.6108*square_N;
    int square_N_y = 1*square_N;
    for (int i = 0; i <= square_N_x; i++) {
        float x = -1+(i*2.f)/square_N_x;
        for (int j = 0; j <= square_N_y; j++) {
            float y = -1+(j*2.f)/square_N_y;
	    all_points.push_back(vec2(rect_x_scale*x,rect_y_scale*y));
        }
    }


    #endif

    // Write to a .node file.
    // FILE *node_file = fopen(DATA "circle_cloud.node", "w+");
    // fprintf(node_file, "%zu 2 0 0\n", all_points.size()); //num_vertices, num_dimensions=2, num_attributes, num_boundary_markers=0,1
    // for (int i = 0; i < all_points.size(); i++) {
    //     fprintf(node_file, "%d %.6f %.6f\n", i, all_points[i].x(), all_points[i].y());
    // }
    // fclose(node_file);

    // Triangulate this point cloud.
    std::string triswitches = "zpcj";
    double *p_mem = (double *) malloc(2*sizeof(double)*all_points.size());
    for (int i = 0; i < all_points.size(); i++) {
        p_mem[2*i] = all_points[i].x();
        p_mem[2*i+1] = all_points[i].y();
    }
    int *segment_mem = (int *) malloc(sizeof(int)*segment_pieces.size());
    for (int i = 0; i < segment_pieces.size(); i++) segment_mem[i] = segment_pieces[i];

    struct triangulateio in = {0};
    memset(&in, 0, sizeof(struct triangulateio));
    in.pointlist = p_mem;
    in.numberofpoints = all_points.size();
    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;

    in.numberofsegments = segment_pieces.size()/2;
    in.segmentlist = segment_mem;
    in.segmentmarkerlist = nullptr;
    // in.numberofholes = 0;
    double *hole_mem = (double *) malloc(sizeof(double)*2);
    in.numberofholes = 1;
    in.holelist = hole_mem;
    in.holelist[0] = obstruction_position.x();
    in.holelist[1] = obstruction_position.y();
    
    in.numberofregions = 0;

    struct triangulateio out = {0};
    memset(&out, 0, sizeof(struct triangulateio));

    triangulate(&triswitches[0], &in, &out, nullptr);

    

    printf("%d\n", out.numberofpoints);
    printf("%d\n", out.numberoftriangles);
    auto vertices = std::vector<Vertex>(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++) {
        printf("%.6f %.6f\n", out.pointlist[2*i+0], out.pointlist[2*i+1]);
        auto v = geom->mesh.add_vertex();
        geom->position[v] = Eigen::Vector3f(out.pointlist[2*i+0], 0, out.pointlist[2*i+1]);
        vertices[i] = v;
    }
    for (int i = 0; i < out.numberoftriangles; i++) {
        int a = out.trianglelist[3*i+0];
        int b = out.trianglelist[3*i+1];
        int c = out.trianglelist[3*i+2];
        printf("%d %d %d\n", a,b,c);
        geom->mesh.add_triangle(vertices[a], vertices[b], vertices[c]);
    }
    // getchar();

    geom->mesh.lock();
    free(p_mem);
    free(segment_mem);
    free(hole_mem);
    return geom;
}
