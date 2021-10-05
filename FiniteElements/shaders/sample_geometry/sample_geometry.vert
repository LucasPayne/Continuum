#version 420

//Position is normalized to [-1,1]^2.
layout (location = 0) in vec2 position;
layout (location = 1) in float value;

out float f;

void main(void)
{
    gl_Position = vec4(position, 0, 1);
    f = value;
}
