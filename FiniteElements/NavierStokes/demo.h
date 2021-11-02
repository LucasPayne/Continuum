#ifndef HEADER_DEFINED_DEMO
#define HEADER_DEFINED_DEMO
#include "core.h"
#include "NavierStokes/NavierStokesSolver.h"
#include "CameraController.h"

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
