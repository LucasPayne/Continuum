

Eigen::VectorXd Solver::pressure_gradient_source()
{
    auto source = Eigen::VectorXd(2*N_u);
    for (int i = 0; i < 2*N_u; i++) source[i] = 0.;

    // Trial psi^u at vertices.
    for (auto v : geom.mesh.vertices()) {
        if (v.on_boundary()) continue;
        auto v_pos = geom.position[v];

        double integral_x = 0.;
        double integral_y = 0.;
        // For each triangle.
        auto start = v.halfedge();
        auto he = start;
        do {
            auto tri = he.face();
            auto vp = he.next().vertex();
            auto vpp = he.next().next().vertex();
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            double R = -0.5/geom.triangle_area(tri);

            vec2 val_v = vec2(0.,0.);

            // Integrate psi_u at v against p[v]*phi_p at v.
            // zero
            
            // val_v = p[v] * R * ((-1./6.)*K1 + (-1./6.)*K2); ---incorrect
            // integral_x += vec2::dot(val_v, K1);
            // integral_y += vec2::dot(val_v, K2);
            
            // Integrate psi_u at v against phi_p at vp.
            // zero
            
            // Integrate psi_u at v against phi_p at vpp.
            // zero
        
            he = he.twin().next();
        } while (he != start);

        int v_index = interior_vertex_indices[v];
        source[2*v_index + 0] = integral_x;
        source[2*v_index + 1] = integral_y;
    }
    // Trial psi^u at midpoints.
    for (auto edge : geom.mesh.edges()) {
        if (edge.on_boundary()) continue;
        Halfedge hes[2] = {edge.a(), edge.b()};

        double integral_x = 0.;
        double integral_y = 0.;
        // For the two incident triangles.
        for (int t = 0; t < 2; t++) {
            // printf("t = %d\n------------------------------------------------------------\n", t);
            // Define terms.
            auto he = hes[t];
            auto tri = he.face();
            // Triangle vertices.
            auto v = he.next().tip(); // v is the opposite vertex.
            auto vp = he.tip();
            auto vpp = he.vertex();

            int v_index = interior_vertex_indices[v];
            int vp_index = interior_vertex_indices[vp];
            int vpp_index = interior_vertex_indices[vpp];
            auto v_pos = geom.position[v];
            auto vp_pos = geom.position[vp];
            auto vpp_pos = geom.position[vpp];
            // Triangle side vectors.
            auto vec2_extract = [](Eigen::Vector3f evec) { return vec2(evec.x(), evec.z()); };
            vec2 K1 = vec2_extract(v_pos - vpp_pos);
            vec2 K2 = vec2_extract(vp_pos - v_pos);
            vec2 K3 = vec2_extract(vpp_pos - vp_pos);
            //double R = -1./(12*geom.triangle_area(tri));
            double R = 1./(12*geom.triangle_area(tri)); // ?--------------------------------------------------------------------------------

            vec2 val_v = vec2(0.,0.);
            // Integrate psi_u at edge_110 against phi_p at 002 (v).
            val_v = p[v] * R * K3.perp();
            integral_x += val_v.x();
            integral_y += val_v.y();

            // Integrate psi_u at edge_110 against phi_p at 200 (vp).
            val_v = p[vp] * R * K1.perp();
            integral_x += val_v.x();
            integral_y += val_v.y();
            
            // Integrate psi_u at edge_110 against phi_p at 002 (vpp).
            val_v = p[vpp] * R * K2.perp();
            integral_x += val_v.x();
            integral_y += val_v.y();
        }
        int index = num_interior_vertices + interior_midpoint_indices[edge];
        source[2*index + 0] = integral_x;
        source[2*index + 1] = integral_y;
    }
    

    return source;
}
