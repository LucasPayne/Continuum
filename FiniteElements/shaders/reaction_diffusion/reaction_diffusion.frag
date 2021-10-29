#version 420
in vec3 f_position;
in float f_u;
in float f_v;

out vec4 color;

void main(void)
{
    const float mul = 3000;
    if (f_position.x > 0) {
        vec3 up_color = vec3(0,1,0);
        vec3 down_color = vec3(0,0,1);
        if (f_u > 0) {
            color = vec4(mul*up_color*abs(f_u), 1);
        } else {
            color = vec4(mul*down_color*abs(f_v), 1);
        }
    } else {
        vec3 up_color = vec3(1,1,0);
        vec3 down_color = vec3(0,1,1);
        if (f_v > 0) {
            color = vec4(mul*up_color*abs(f_v), 1);
        } else {
            color = vec4(mul*down_color*abs(f_v), 1);
        }
    }
}
