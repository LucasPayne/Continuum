#version 420

layout (location = 0) in vec3 position;
layout (location = 1) in float value;

out VS_OUT {
    vec3 position;
    float value;
} vs_out;

// uniform mat4x4 mvp_matrix;

void main(void)
{
    vs_out.position = position;
    vs_out.value = value;
    
    // gl_Position = mvp_matrix * vec4(position, 1);
}

