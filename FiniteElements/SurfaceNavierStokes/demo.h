#ifndef HEADER_DEFINED_DEMO
#define HEADER_DEFINED_DEMO
#include "core.h"
#include "SurfaceNavierStokes/SurfaceNavierStokesSolver.h"
#include "CameraController.h"

struct Demo : public IBehaviour {
    Demo();
    void init();

    CameraController *controller;

    void update();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);

    SurfaceNavierStokesSolver *solver;
    SurfaceGeometry *geom;

    void save_solution(std::string filename);

    std::tuple<Face, vec3> traverse(Face tri, vec3 origin, vec3 shift);
};

#endif // HEADER_DEFINED_DEMO
