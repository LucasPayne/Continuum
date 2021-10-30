#version 420

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 velocity;
layout (location = 2) in float pressure;

out VS_OUT {
    vec2 position;
    vec2 velocity;
    float pressure;
} vs_out;

void main(void)
{
    vs_out.position = position;
    vs_out.velocity = velocity;
    vs_out.pressure = pressure;
}

