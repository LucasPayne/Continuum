#include "NavierStokes/NavierStokesSolver.h"
#include "NavierStokes/core.h"

void NavierStokesSolver::compute_pressure_residual(P2Attachment<double> &pressure_residual)
{
    // Make sure that the residual begins at zero.
    for (auto v : geom.mesh.vertices()) {
        pressure_residual[v] = 0.;
    }

    // For each trial pressure basis function at a vertex node.
    for (auto v : geom.mesh.vertices()) {
        auto v_pos = geom.position[v];

        double integral = 0.;

        // For each adjacent triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            // Define terms.
            Face tri = he.face();
            Vertex vp = he.next().vertex();
            Vertex vpp = he.next().tip();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            Edge edge_101 = he.edge();
            Edge edge_110 = he.next().edge();
            Edge edge_011 = he.next().next().edge();
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);

            integral += vec2::dot(velocity[v], K1/6 + K2/6);
            if (!edge_110.on_boundary()) {
                integral += vec2::dot(velocity[edge_110], -K1/6 - K2/6);
            }
            if (!edge_011.on_boundary()) {
                integral += vec2::dot(velocity[edge_110], K1/6 - K2/6);
            }
            if (!edge_101.on_boundary()) {
                integral += vec2::dot(velocity[edge_110], -K1/6 + K2/6);
            }

            he = he.twin().next();
        } while (he != start);

        pressure_residual[v] += integral;
    }
}
