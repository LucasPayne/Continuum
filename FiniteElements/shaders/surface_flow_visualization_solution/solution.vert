#version 420

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 velocity;
layout (location = 2) in float pressure;

out VS_OUT {
    vec3 position;
    vec3 velocity;
    float pressure;
} vs_out;

void main(void)
{
    vs_out.position = position;
    vs_out.velocity = velocity;
    vs_out.pressure = pressure;
}

