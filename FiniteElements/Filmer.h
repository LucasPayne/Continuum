#ifndef FILMER_H
#define FILMER_H
#include "core.h"

struct Filmer : public IBehaviour
{
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
    int start_frame;
    int film_frame;
    int filming;
    std::string base_path;

    std::function<void(int)> post_frame_callback;

    Filmer() {}
    Filmer(std::string _base_path, std::function<void(int)> _post_frame_callback, int _start_frame=0);
    void init();
    void update();
    void post_render_update();
    void mouse_handler(MouseEvent e);
    void keyboard_handler(KeyboardEvent e);
};

#endif // FILMER_H
