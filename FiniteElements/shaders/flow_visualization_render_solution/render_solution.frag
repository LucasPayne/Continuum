#version 420

in VS_OUT {
    vec3 position;
    vec2 uv;
} fs_in;

uniform sampler2D solution;
uniform int mode;

out vec4 color;

void main(void)
{
    float velocity_multiplier = 1;

    vec4 val = texture(solution, fs_in.uv);
    if (mode == 0 || mode == 1) {
        // Velocity component
        float v = val[mode];
        if (v > 0) {
            color = vec4(v, 0,0,1);
        } else {
            color = vec4(0,0,abs(v),1);
        }
        color = vec4(vec3(1-exp(-0.1*length(velocity_multiplier * val.xy))), 1);
    } else if (mode == 2) {
        // Pressure
        color = vec4(vec3(val[mode]), 1);
    } else if (mode == 3) {
        // Divergence
        color = vec4(vec3(val[mode]), 1);
    }
}
