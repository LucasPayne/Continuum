#version 420

in TES_OUT {
    vec3 barycentric;
    flat vec2 velocities[6];
    float pressure;
} fs_in;

out vec4 color; // render-to-texture, (velocity_x, velocity_y, pressure, 1)


// Quadratic basis functions
float q0(in float x, in float y, in float z)
{
    return x - 2*x*y - 2*x*z;
}
float q1(in float x, in float y, in float z)
{
    return 4*x*y;
}
float q2(in float x, in float y, in float z)
{
    return y - 2*y*z - 2*y*x;
}
float q3(in float x, in float y, in float z)
{
    return 4*y*z;
}
float q4(in float x, in float y, in float z)
{
    return z - 2*z*x - 2*z*y;
}
float q5(in float x, in float y, in float z)
{
    return 4*z*x;
}

void main(void)
{
    float x = fs_in.barycentric.x;
    float y = fs_in.barycentric.y;
    float z = fs_in.barycentric.z;

    float q[6];
    q[0] = q0(x,y,z);
    q[1] = q1(x,y,z);
    q[2] = q2(x,y,z);
    q[3] = q3(x,y,z);
    q[4] = q4(x,y,z);
    q[5] = q5(x,y,z);

    // Evaluate velocity.
    vec2 velocity = vec2(0,0);
    for (int i = 0; i < 6; i++) {
        velocity += q[i] * fs_in.velocities[i];
    }
    color = vec4(velocity, fs_in.pressure, 1);
}
