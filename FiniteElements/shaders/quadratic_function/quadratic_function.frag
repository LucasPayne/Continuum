#version 420

in TES_OUT {
    vec2 plane_position;
    vec3 barycentric;
    float value;
} fs_in;

out vec4 color;


// Basis functions
float p0(in float x, in float y, in float z)
{
    return x - 2*x*y - 2*x*z;
}
float p1(in float x, in float y, in float z)
{
    return y - 2*y*z - 2*y*x;
}
float p2(in float x, in float y, in float z)
{
    return z - 2*z*x - 2*z*y;
}
float p3(in float x, in float y, in float z)
{
    return 4*x*y;
}
float p4(in float x, in float y, in float z)
{
    return 4*y*z;
}
float p5(in float x, in float y, in float z)
{
    return 4*z*x;
}

void main(void)
{
    // Draw manufactured solution on the left.
    float x = fs_in.plane_position.x;
    float y = fs_in.plane_position.y;
    float sol = exp(-2*(x*x + y*y));

#if 1
    // Visualizing error.
    color = vec4(vec3(fs_in.value), 1);
    if (x < 0) color = vec4(vec3(sol), 1);
#else
    // Computing error.
    color = vec4(vec3((sol - fs_in.value)*(sol - fs_in.value)), 1);
#endif
}

