


struct Demo : public IBehaviour {
    Demo();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);

    void recreate_solver();

    SurfaceGeometry *geom;
    Solver *solver;

    int mesh_N;
};


void Demo::recreate_solver()
{
    if (geom != nullptr) delete geom;
    // geom = circle_mesh(mesh_N, false);
    auto sq_mesh = SquareMesh(mesh_N);
    geom = sq_mesh.geom;

    if (solver != nullptr) delete solver;
    solver = new Solver(*geom);

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
    }
}


void Demo::post_render_update()
{
    world->graphics.paint.wireframe(*geom, mat4x4::identity(), 0.001);
}
