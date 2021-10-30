#ifndef HEADER_DEFINED_P2_P1
#define HEADER_DEFINED_P2_P1
// #include "NavierStokes/core.h"


// A P2 triangle mesh discretization has nodal points at vertices and edge midpoints.
// Functions on these nodes are aggregated into a P2Attachment.
// Wrapper for both Vertex and Edge, for accessing a P2Attachment.
template <typename T>
class P2Attachment;
class P2Element {
public:
    P2Element(Vertex v) :
        vertex{v}, m_is_vertex{true}
    {}
    P2Element(Edge e) :
        edge{e}, m_is_vertex{false}
    {}
    inline bool is_vertex() const { return m_is_vertex; }
    inline bool on_boundary() {
        if (is_vertex()) return vertex.on_boundary();
        return edge.on_boundary();
    }
// private: //how can a template class be made a friend?
    Vertex vertex;
    Edge edge;
    bool m_is_vertex;
};
template <typename T>
class P2Attachment {
public:
    P2Attachment(SurfaceMesh &_mesh) :
        vertex_attachment(_mesh),
        edge_attachment(_mesh),
        mesh{_mesh}
    {}
    inline T operator[](Vertex v) const {
        return vertex_attachment[v];
    }
    inline T operator[](Edge e) const {
        return edge_attachment[e];
    }
    inline T &operator[](P2Element element) const {
        if (element.is_vertex()) {
            return vertex_attachment[element.vertex];
        } else {
            return edge_attachment[element.edge];
        }
    }
    inline T &operator[](Vertex v) {
        return vertex_attachment[v];
    }
    inline T &operator[](Edge e) {
        return edge_attachment[e];
    }
    inline T &operator[](P2Element element) {
        if (element.is_vertex()) {
            return vertex_attachment[element.vertex];
        } else {
            return edge_attachment[element.edge];
        }
    }
    VertexAttachment<T> vertex_attachment;
    EdgeAttachment<T> edge_attachment;
private:
    SurfaceMesh &mesh;
};


// A P1 triangle mesh discretization has nodal points at vertices.
// Functions on these nodes are aggregated into a P1Attachment.
template <typename T>
using P1Attachment = VertexAttachment<T>;

#endif // HEADER_DEFINED_P2_P1
