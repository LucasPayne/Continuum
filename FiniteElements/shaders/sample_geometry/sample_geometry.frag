#version 420

in vec2 f_position;
in float f;

out vec4 color;

void main(void)
{
    float x = f_position.x;
    float y = f_position.y;

    float solution = exp(-2*(x*x + y*y));
    // float solution = x*x - y*y;

    float diff = f - solution;


    // if (f_position.x > 0) {
    //     color = vec4(vec3(f),1);
    // } else {
    //     color = vec4(vec3(solution),1);
    // }

    // color = vec4(30*abs(diff), 0, 0,1);
    // color = vec4(15*abs(diff), 0, 0,1);

    // color = vec4(30*diff*diff, 0, 0,1);
    color = vec4(diff*diff, 0, 0,1);
}
