#version 420

in vec2 f_uv;

out vec4 color;


uniform sampler2D tex;
uniform int mode;


void main(void)
{
    // color = vec4(1,0,0,1);
    // color = vec4(texture(tex, f_uv).rgb, 1);

    // vec3 val = texture(tex, f_uv).rgb;
    vec4 val = texture(tex, vec2(1-f_uv.x, f_uv.y)).rgba;
    vec2 velocity = val.rg;
    float pressure = val.b;
    float div_u = val.a;

    const float p_mul = 1.f;
    const float v_mul = 10.f;
    const float div_mul = 4.f;

    if (mode == 3) { // Pressure
        color = vec4(vec3(p_mul * pressure), 1);
    } else if (mode == 2) { // Velocity x
        if (velocity.x < 0) {
            color = vec4(v_mul * abs(velocity.x), 0,0,1);
        } else {
            color = vec4(0,0, v_mul * abs(velocity.x), 1);
        }
    } else if (mode == 1) { // Velocity y
        if (velocity.y < 0) {
            color = vec4(v_mul * abs(velocity.y), 0,0,1);
        } else {
            color = vec4(0,0, v_mul * abs(velocity.y), 1);
        }
    } else if (mode == 0) { // Velocity divergence
        color = vec4(0,0, v_mul * abs(velocity.y), 1);
        if (div_u < 0) {
            color = vec4(div_mul*vec2(abs(div_u)),0,1);
        } else {
            color = vec4(0,div_mul*vec2(abs(div_u)),1);
        }
    }
    // if (val.a == 0) color = vec4(1); //...
}
