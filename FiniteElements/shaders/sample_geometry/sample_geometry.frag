#version 420

in float f;

out vec4 color;

void main(void)
{
    color = vec4(f, 0,0,1);
}
