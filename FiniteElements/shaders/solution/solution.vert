#version 420

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 velocity;
layout (location = 2) in float pressure;
layout (location = 3) in float divergence;
layout (location = 4) in float divergence_P2;

out VS_OUT {
    vec2 position;
    vec2 velocity;
    float pressure;
    float divergence;
    float divergence_P2;
} vs_out;

void main(void)
{
    vs_out.position = position;
    vs_out.velocity = velocity;
    vs_out.pressure = pressure;
    vs_out.divergence = divergence;
    vs_out.divergence_P2 = divergence_P2;
}

