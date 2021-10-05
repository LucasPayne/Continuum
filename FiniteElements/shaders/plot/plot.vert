#version 420

in vec2 uv;

uniform vec2 bottom_left;
uniform vec2 top_right;

out vec2 f_uv;

void main(void)
{
    vec2 position = vec2(-1,-1) + 2*bottom_left + 2*(top_right - bottom_left)*uv;
    gl_Position = vec4(position, 0, 1);
    f_uv = uv;
}
