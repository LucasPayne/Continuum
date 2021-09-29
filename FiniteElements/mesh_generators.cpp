
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


