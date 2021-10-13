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
    vec3 val = texture(tex, vec2(1-f_uv.x, f_uv.y)).rgb; // swap x
    vec2 velocity = val.rg;
    float pressure = val.b;

    const float p_mul = 1.f;
    const float v_mul = 4.f;

    if (mode == 0) {
        color = vec4(vec3(p_mul * pressure), 1);
    } else if (mode == 1) {
        if (velocity.x < 0) {
            color = vec4(v_mul * abs(velocity.x), 0,0,1);
        } else {
            color = vec4(0,0, v_mul * abs(velocity.x), 1);
        }
    } else if (mode == 2) {
        if (velocity.y < 0) {
            color = vec4(v_mul * abs(velocity.y), 0,0,1);
        } else {
            color = vec4(0,0, v_mul * abs(velocity.y), 1);
        }
    }
}
