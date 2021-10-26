#include "NavierStokes/NavierStokesSolver.h"

std::vector<TopLeftEntry> NavierStokesSolver::compute_linear_term_matrix_top_left()
{
    auto coeffs = std::vector<TopLeftEntry>();

    double inv_dt = 1./m_current_time_step_dt;

    // For each velocity trial basis function psi^u
    //------------------------------------------------------------
    // For each psi^u based on a vertex.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

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
            Edge edge_110 = he.next().edge();
            Edge edge_011 = he.next().next().edge();
            Edge edge_101 = he.edge();
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;
            double tri_area = geom.triangle_area(tri);

            /*--------------------------------------------------------------------------------
                Time step term
            --------------------------------------------------------------------------------*/
            P2Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            element_weights = {
                1./60., -1./360., -1./360., -1./90., 0, 0
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                double val = 2 * tri_area * element_weights[i] * inv_dt;
                coeffs.emplace_back(P2Element(v), 0, P2Element(elements[i]), 0, val);
                coeffs.emplace_back(P2Element(v), 1, P2Element(elements[i]), 1, val);
            }

            /*--------------------------------------------------------------------------------
                Viscosity term
            --------------------------------------------------------------------------------*/
            double C = m_kinematic_viscosity * 1.0/(4.0*geom.triangle_area(tri)); //---------half?

            double val = 0.;
            // Diagonal term.
            val = C * 0.5 * K3.dot(K3);
            coeffs.emplace_back(P2Element(v), 0, P2Element(v), 0, val);
            coeffs.emplace_back(P2Element(v), 1, P2Element(v), 1, val);

            // vp contribution.
            val = -(1./6.) * C * K1.dot(K3);
            if (!vp.on_boundary()) {
                coeffs.emplace_back(P2Element(v), 0, P2Element(vp), 0, val);
                coeffs.emplace_back(P2Element(v), 1, P2Element(vp), 1, val);
            }
            
            // vpp contribution.
            val = -(1./6.) * C * K2.dot(K3);
            if (!vpp.on_boundary()) {
                coeffs.emplace_back(P2Element(v), 0, P2Element(vpp), 0, val);
                coeffs.emplace_back(P2Element(v), 1, P2Element(vpp), 1, val);
            }

            // edge_101 contribution.
            val = (2./3.)*C*K1.dot(K3);
            if (!edge_101.on_boundary()) {
                coeffs.emplace_back(P2Element(v), 0, P2Element(edge_101), 0, val);
                coeffs.emplace_back(P2Element(v), 1, P2Element(edge_101), 1, val);
            }
            
            // edge_011 contribution.
            val = (2./3.)*C*K2.dot(K3);
            if (!edge_011.on_boundary()) {
                coeffs.emplace_back(P2Element(v), 0, P2Element(edge_011), 0, val);
                coeffs.emplace_back(P2Element(v), 1, P2Element(edge_011), 1, val);
            }

            he = he.twin().next();
        } while (he != start);
    }
    // For each psi^u based at an edge midpoint.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.vertex();
            auto vpp = he.tip();
            auto edge_110 = he.edge();
            auto edge_011 = he.next().edge();
            auto edge_101 = he.next().next().edge();
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto K1 = v_pos - vpp_pos;
            auto K2 = vp_pos - v_pos;
            auto K3 = vpp_pos - vp_pos;
            double tri_area = geom.triangle_area(tri);

            /*--------------------------------------------------------------------------------
                Time step term
            --------------------------------------------------------------------------------*/
            P2Element elements[6] = {
                v, vp, vpp, edge_110, edge_011, edge_101
            };
            auto element_weights = std::array<double,6>();
            element_weights = {
                -1./90., 0, 0, 4./45., 2./45., 2./45.
            };
            for (int i = 0; i < 6; i++) {
                if (elements[i].on_boundary()) continue;
                double val = 2 * tri_area * element_weights[i] * inv_dt;
                coeffs.emplace_back(P2Element(edge), 0, P2Element(elements[i]), 0, val);
                coeffs.emplace_back(P2Element(edge), 1, P2Element(elements[i]), 1, val);
            }

            /*--------------------------------------------------------------------------------
                Viscosity term
            --------------------------------------------------------------------------------*/
            double C = m_kinematic_viscosity * 1.0/(4.0*geom.triangle_area(tri)); //-----half?

            double val = 0.;
            // 110, 110
            val = 4.*C/3. * (K3.dot(K3) - K1.dot(K2));
            coeffs.emplace_back(P2Element(edge), 0, P2Element(edge_110), 0, val);
            coeffs.emplace_back(P2Element(edge), 1, P2Element(edge_110), 1, val);

            // 110, 011
            val = 4.*C/3. * (K1.dot(K3));
            if (!edge_011.on_boundary()) {
                coeffs.emplace_back(P2Element(edge), 0, P2Element(edge_011), 0, val);
                coeffs.emplace_back(P2Element(edge), 1, P2Element(edge_011), 1, val);
            }
            
            // 110, 101
            val = 4.*C/3. * (K2.dot(K3));
            if (!edge_101.on_boundary()) {
                coeffs.emplace_back(P2Element(edge), 0, P2Element(edge_101), 0, val);
                coeffs.emplace_back(P2Element(edge), 1, P2Element(edge_101), 1, val);
            }
            
            // 110, 200
            val = 2.*C/3. * (K1.dot(K2));
            if (!vp.on_boundary()) {
                coeffs.emplace_back(P2Element(edge), 0, P2Element(vp), 0, val);
                coeffs.emplace_back(P2Element(edge), 1, P2Element(vp), 1, val);
            }
            // 110, 020
            if (!vpp.on_boundary()) {
                coeffs.emplace_back(P2Element(edge), 0, P2Element(vpp), 0, val);
                coeffs.emplace_back(P2Element(edge), 1, P2Element(vpp), 1, val);
            }
        }
    }

    return coeffs;
}
