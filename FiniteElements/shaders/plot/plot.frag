#version 420

in vec2 f_uv;

out vec4 color;


uniform sampler2D tex;


void main(void)
{
    // color = vec4(1,0,0,1);
    color = vec4(texture(tex, f_uv).rgb, 1);
}
