#version 330 core

in vec4 vertex_color;
uniform vec4 model_color = vec4(1.0, 1.0, 1.0, 1.0);

out vec4 fragment_color;

void main() {
    // Use vertex color if it has meaningful values, otherwise use uniform color
    if (vertex_color.rgb != vec3(1.0, 1.0, 1.0) || vertex_color.a != 1.0) {
        fragment_color = vertex_color;
    }
    else {
        fragment_color = model_color;
    }
}
