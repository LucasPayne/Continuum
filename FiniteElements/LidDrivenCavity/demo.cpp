


struct Demo : public IBehaviour {
    Demo();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);

    void recreate_solver();

    SurfaceGeometry *geom;
    Solver *solver;

    int mesh_N;
    double mu;
};


void Demo::recreate_solver()
{
    if (geom != nullptr) delete geom;
    // geom = circle_mesh(mesh_N, false);
    auto sq_mesh = SquareMesh(mesh_N);
    geom = sq_mesh.geom;

    if (solver != nullptr) delete solver;
    solver = new Solver(*geom, mu);

    // Set the lid boundary condition explicitly, on the vertex and midpoint sample points.
    // (This is to avoid possible errors at corners if the boundary condition was specified with a function.)
    for (int i = 0; i < mesh_N+1; i++) {
        solver->u_boundary[sq_mesh.vertex(i,0)] = vec2(1,0);
    }
    for (int i = 0; i < mesh_N; i++) {
        auto v1 = sq_mesh.vertex(i,0);
        auto v2 = sq_mesh.vertex(i+1,0);
        solver->u_boundary[geom->mesh.vertices_to_edge(v1, v2)] = vec2(1,0);
    }
}


Demo::Demo()
{
    mesh_N = 4;
    mu = 1.0;
    recreate_solver();
}

void Demo::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_P) {
            mesh_N += 1;
            recreate_solver();
        }
        if (e.key.code == KEY_O) {
            mesh_N -= 1;
            if (mesh_N < 2) mesh_N = 2;
            recreate_solver();
        }
        if (e.key.code == KEY_R) {
            solver->solve();
        }
        if (e.key.code == KEY_T) {
            solver->write_sparsity_pattern = true;
            solver->solve();
        }
    }
}


void Demo::post_render_update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);

    double u_multiplier = 0.1;

    auto draw_u_vec = [&](Eigen::Vector3f pos, vec2 u_val) {
        world->graphics.paint.sphere(eigen_to_vec3(pos), 0.01, vec4(0,0,1,1));
        world->graphics.paint.line(eigen_to_vec3(pos), eigen_to_vec3(pos) + u_multiplier*vec3(u_val.x(), 0.05, u_val.y()), 0.005, vec4(0,0,0,1));
    };

    // // Draw boundary velocity.
    // for (auto v : geom->mesh.vertices()) {
    //     vec3 pos = eigen_to_vec3(geom->position[v]);
    //     vec2 vec = solver->u_boundary[v];
    // }
    // for (auto e : geom->mesh.edges()) {
    //     vec3 pos = eigen_to_vec3(solver->midpoints[e]);
    //     vec2 vec = solver->u_boundary[e];
    //     world->graphics.paint.line(pos, pos + 0.06*vec3(vec.x(), 0, vec.y()), 0.01, vec4(1,0,0,1));
    // }

    // Draw the solution.
    for (auto v : geom->mesh.vertices()) {
        draw_u_vec(geom->position[v], solver->u[v]);
    }
    for (auto e : geom->mesh.edges()) {
        draw_u_vec(solver->midpoints[e], solver->u[e]);
    }
}
