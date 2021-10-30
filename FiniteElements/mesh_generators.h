#ifndef HEADER_DEFINED_MESH_GENERATORS
#define HEADER_DEFINED_MESH_GENERATORS

// Create a mesh.
SurfaceGeometry *circle_mesh(int N, bool random);

// [-1,1]^2
SurfaceGeometry *square_mesh(int N);

SurfaceGeometry *square_mesh_with_square_hole(int N);

// Return with a data structure which can be used to access the square mesh with rectangular indexing.
class SquareMesh
{
public:
    SquareMesh(int _N);

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
);

#endif // HEADER_DEFINED_MESH_GENERATORS
