#version 420

in TES_OUT {
    vec3 barycentric;
    flat vec2 velocities[6];
    float pressure;
    float divergence;
    float divergence_P2[6];
} fs_in;

out vec4 color; // render-to-texture, (velocity_x, velocity_y, pressure, 1)

uniform int show_div_P2;


void main(void)
{
    float x = fs_in.barycentric.x;
    float y = fs_in.barycentric.y;
    float z = fs_in.barycentric.z;

    // Evaluate quadratic basis functions
    float q[6] = {
        x - 2*x*y - 2*x*z,
        4*x*y,
        y - 2*y*z - 2*y*x,
        4*y*z,
        z - 2*z*x - 2*z*y,
        4*z*x,
    };

    // Evaluate velocity.
    vec2 velocity = vec2(0,0);
    for (int i = 0; i < 6; i++) {
        velocity += q[i] * fs_in.velocities[i];
    }
    color = vec4(velocity, fs_in.pressure, fs_in.divergence);

    // Evaluate P2 divergence.
    float divergence_P2 = 0.;
    for (int i = 0; i < 6; i++) {
        divergence_P2 += q[i] * fs_in.divergence_P2[i];
    }
    if (show_div_P2 == 1) {
        color = vec4(velocity, fs_in.pressure, divergence_P2);
    } else {
        color = vec4(velocity, fs_in.pressure, fs_in.divergence);
    }
}
