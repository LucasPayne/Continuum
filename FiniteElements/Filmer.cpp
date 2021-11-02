#include "Filmer.h"

Filmer::Filmer(std::string _base_path, std::function<void(int)> _post_frame_callback, int _start_frame) :
    base_path{_base_path},
    post_frame_callback{_post_frame_callback},
    start_frame{_start_frame}
{
}

void Filmer::init()
{
    filming = false;
}

void Filmer::update()
{
    if (world->input.keyboard.down(KEY_G)) {
        // Draw the screenshot rectangle.
        std::vector<vec2> ps = {
            vec2(f_screenshot_blx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_bly),
            vec2(f_screenshot_trx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_try),
            vec2(f_screenshot_blx, f_screenshot_bly)
        };
        world->graphics.paint.chain_2D(ps, 1, vec4(1,0,0,1));
    }
}
void Filmer::post_render_update()
{
    if (filming) {
        post_frame_callback(film_frame);
        take_screenshot();
        film_frame += 1;
    }
}

void Filmer::keyboard_handler(KeyboardEvent e)
{
    if (e.action == KEYBOARD_PRESS) {
        if (e.key.code == KEY_9) {
            filming = true;
            film_frame = start_frame;
        }
    }
}


void Filmer::mouse_handler(MouseEvent e)
{
    if (e.action == MOUSE_BUTTON_PRESS) {
        if (!filming ) {
            // Set the screenshot rectangle.
            if (e.button.code == MOUSE_LEFT) {
                screenshot_blx = int(world->graphics.window_viewport.w * e.cursor.x);
                screenshot_bly = int(world->graphics.window_viewport.h * e.cursor.y);
                f_screenshot_blx = e.cursor.x;
                f_screenshot_bly = e.cursor.y;
            }
            if (e.button.code == MOUSE_RIGHT) {
                screenshot_trx = int(world->graphics.window_viewport.w * e.cursor.x);
                screenshot_try = int(world->graphics.window_viewport.h * e.cursor.y);
                f_screenshot_trx = e.cursor.x;
                f_screenshot_try = e.cursor.y;
            }
        }
    }
}

void Filmer::take_screenshot()
{
    std::string pre = base_path + "." + std::to_string(film_frame);
    world->graphics.screenshot(pre + ".ppm",
                               screenshot_blx,
                               screenshot_bly,
                               screenshot_trx - screenshot_blx,
                               screenshot_try - screenshot_bly);
}

