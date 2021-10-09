#version 420

in TES_OUT {
    vec3 barycentric;
    vec3 normal;
} fs_in;

out vec4 color;

void main(void)
{
    vec3 light_dir = normalize(vec3(0,-1,0));
    vec3 n = fs_in.normal;

    vec3 face_diffuse = vec3(0.93);
    const float lightness = 0.5;
    vec4 face_color = vec4(face_diffuse * (lightness + (1-lightness)*max(0, dot(light_dir, n))), 1);

    vec4 wireframe_color = vec4(0,0,0,1);

    float x = fs_in.barycentric.x;
    float y = fs_in.barycentric.y;
    float z = fs_in.barycentric.z;

    color = face_color;
    
    float b = 0.01;
    if (x < b || y < b || z < b) {
        color = wireframe_color;
    }
    
}

