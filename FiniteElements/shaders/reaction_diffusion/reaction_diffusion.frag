#version 420
in vec3 f_position;
in float f_value;

out vec4 color;

void main(void)
{
    vec3 up_color = vec3(0,1,0);
    vec3 down_color = vec3(0,0,1);

    if (f_value > 0) {
        color = vec4(up_color*abs(f_value), 1);
    } else {
        color = vec4(down_color*abs(f_value), 1);
    }
}
