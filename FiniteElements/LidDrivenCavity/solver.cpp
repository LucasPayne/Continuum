
// A P2 triangle mesh discretization has nodal points at vertices and edge midpoints.
// Functions on these nodes are aggregated into a P2Attachment.
template <typename T>
class P2Attachment {
public:
    P2Attachment(SurfaceMesh &_mesh) :
        vertex_attachment(_mesh),
        edge_attachment(_mesh),
        mesh{_mesh}
    {}
    inline T operator[](Vertex v) {
        return vertex_attachment[v];
    }
    inline T operator[](Edge e) {
        return edge_attachment[e];
    }
    VertexAttachment<T> vertex_attachment;
    EdgeAttachment<T> edge_attachment;
private:
    SurfaceMesh &mesh;
};


struct Solver {
    Solver(SurfaceGeometry &_geom);
    
    // Solution (initializes to 0).
    //------------------------------------------------------------
    // P2 coefficients for velocity u.
    P2Attachment<double> u;
    // P1 coefficients for pressure p.
    VertexAttachment<double> p;

    // Dirichlet boundary condition.
    P2Attachment<vec2> u_boundary;
    // Set boundary condition from a function.
    void set_u_boundary(PlaneVectorField _u_boundary);

    // Additional mesh data.
    //------------------------------------------------------------
    // Flat index ordering of vertices and midpoints.
    VertexAttachment<int> vertex_indices;
    EdgeAttachment<int> midpoint_indices;
    // Store precomputed midpoints for convenience.
    EdgeAttachment<Eigen::Vector3f> midpoints;

    // Mesh properties.
    int num_boundary_vertices;
    int num_interior_vertices;
    int num_boundary_edges;
    int num_interior_edges;

    SurfaceGeometry &geom;
};

Solver::Solver(SurfaceGeometry &_geom) :
    u(_geom.mesh),
    p(_geom.mesh),
    u_boundary(_geom.mesh),
    vertex_indices(_geom.mesh),
    midpoint_indices(_geom.mesh),
    midpoints(_geom.mesh),

    geom{_geom}
{
    num_boundary_vertices = 0;
    num_interior_vertices = 0;
    num_boundary_edges = 0;
    num_interior_edges = 0;
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) {
            vertex_indices[v] = -1;
            num_boundary_vertices += 1;
        } else {
            vertex_indices[v] = num_interior_vertices;
            num_interior_vertices += 1;
        }
    }
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) {
            midpoint_indices[edge] = -1;
            num_boundary_edges += 1;
        } else {
            midpoint_indices[edge] = num_interior_edges;
            num_interior_edges += 1;
        }
    }
    // Compute midpoints.
    for (auto edge : geom.mesh.edges()) {
        auto p = geom.position[edge.a().vertex()];
        auto pp = geom.position[edge.b().vertex()];
        midpoints[edge] = 0.5*p + 0.5*pp;
    }
}


void Solver::set_u_boundary(PlaneVectorField vf)
{
    for (auto v : geom.mesh.vertices()) {
        auto pos = geom.position[v];
        u_boundary[v] = vf(pos.x(), pos.z());
    }
    for (auto e : geom.mesh.edges()) {
        auto pos = midpoints[e];
        u_boundary[e] = vf(pos.x(), pos.z());
    }
}
