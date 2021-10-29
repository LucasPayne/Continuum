#version 420

layout (location = 0) in vec3 position;
layout (location = 1) in float u;
layout (location = 2) in float v;

uniform mat4x4 mvp_matrix;

out vec3 f_position;
out float f_u;
out float f_v;

void main(void)
{
    gl_Position = mvp_matrix * vec4(position, 1);
    f_position = position;
    f_u = u;
    f_v = v;
}
