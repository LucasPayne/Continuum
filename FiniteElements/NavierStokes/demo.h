#ifndef HEADER_DEFINED_DEMO
#define HEADER_DEFINED_DEMO
#include "NavierStokes/core.h"
#include "NavierStokes/NavierStokesSolver.h"

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

    NavierStokesSolver *solver;
    SurfaceGeometry *geom;

    // Mesh selection
    int mesh_mode;
    void recreate_solver();

    // Visualization
    void render_solution_texture();

    bool show_wireframe;
    bool show_vector_field;
    bool show_source;
    
    GLShaderProgram solution_shader; // Render the solution (velocity and pressure) to textures.
    GLuint solution_fbo; // For render-to-texture.
    GLuint solution_texture; // 3 components: r:velocity_x, g:velocity_y, b:pressure.
    GLuint solution_depth_texture; //--- Is this needed for completeness?
    
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
    // Films
    double film_seconds;
    double film_dt;
    double film_num_frames;
    int film_frame;
    int filming;

    // For interactive fluid.
    vec2 source_position;

    void save_solution(std::string filename);
};

#endif // HEADER_DEFINED_DEMO
