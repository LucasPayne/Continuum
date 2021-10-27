#ifndef HEADER_DEFINED_DEMO
#define HEADER_DEFINED_DEMO
#include "ReactionDiffusion/core.h"
#include "ReactionDiffusion/ReactionDiffusion.h"

/*
 * CameraController behaviour.
 * Attach this to an entity (with a transform) to give it mouse view and movement controls.
 * This is intended to be used alongside an attached Camera to create a camera man.
 */
struct CameraController : public IBehaviour {
    float azimuth;
    float angle;

    float strafe_speed;
    float forward_speed;
    float lift_speed;

    bool view_with_mouse;
    float key_view_speed_horizontal;
    float key_view_speed_vertical;

    float min_angle;
    float max_angle;

    #define BASE_MOUSE_SENSITIVITY 1.22
    // mouse_sensitivity multiplies the base sensitivity.
    float mouse_sensitivity;

    CameraController() {}
    inline void lock_angle() {
        if (angle < min_angle) angle = min_angle;
        if (angle > max_angle) angle = max_angle;
    }

    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);
    void init();
    void update();
};

struct Demo : public IBehaviour {
    Demo();
    void init();

    CameraController *controller;

    void update();
    void post_render_update();
    void keyboard_handler(KeyboardEvent e);
    void mouse_handler(MouseEvent e);

    ReactorDiffuser *solver;
    SurfaceGeometry *geom;

    // Visualization
    GLShaderProgram solution_shader;
    bool show_wireframe;
    
    // Screenshots
    int screenshot_blx;
    int screenshot_bly;
    int screenshot_trx;
    int screenshot_try;
    float f_screenshot_blx;
    float f_screenshot_bly;
    float f_screenshot_trx;
    float f_screenshot_try;
    void take_screenshot();
};

#endif // HEADER_DEFINED_DEMO
