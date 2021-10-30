#version 420

in VS_OUT {
    vec3 position;
    vec2 uv;
} fs_in;

uniform sampler2D solution;
uniform int mode;

out vec4 color;

void main(void)
{
    vec4 val = texture(solution, fs_in.uv);
    // color = vec4(vec3(val[mode]), 1);
    color = vec4(vec3(val.z), 1);
}
