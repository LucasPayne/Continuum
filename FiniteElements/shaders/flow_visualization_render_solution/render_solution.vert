#version 420

layout (location = 0) in vec3 position;
layout (location = 1) in vec2 uv;

out VS_OUT {
    vec3 position;
    vec2 uv;
} vs_out;

uniform mat4x4 mvp_matrix;

void main(void)
{
    gl_Position = mvp_matrix * vec4(position, 1);

    vs_out.position = position;
    vs_out.uv = uv;
}

