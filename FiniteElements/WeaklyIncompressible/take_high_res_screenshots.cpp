
void Demo::take_high_res_screenshot(int n)
{
    float x_half_wid = 1.25;
    float extent = 1.05*x_half_wid/4;

    main_camera->projection_matrix = mat4x4( // column major
        1.f/extent,0,0,0,
        0,1.f/(0.566*extent),0,0,
        0,0,1,0,
        0,0,0,1
    );
    float py = 0.5;
    vec3 positions[4] = {
        vec3(-5*x_half_wid/8, py, 0),
        vec3(-x_half_wid/4, py, 0),
        vec3(x_half_wid/4, py, 0),
        vec3(5*x_half_wid/8, py, 0)
    };
    controller->azimuth = M_PI;
    controller->angle = -M_PI/2;
    controller->entity.get<Transform>()->position = positions[n];

}
