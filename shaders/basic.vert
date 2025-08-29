#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec4 color;

uniform mat4 mvp_matrix;
uniform mat4 model_matrix;

out vec4 vertex_color;
out vec3 world_position;

void main() {
    // Transform position to world space
    vec4 world_pos = model_matrix * vec4(position, 1.0);
    world_position = world_pos.xyz;
    
    // Transform to clip space
    gl_Position = mvp_matrix * vec4(position, 1.0);
    
    // Pass color to fragment shader
    vertex_color = color;
}
