#version 330 core

in vec4 vertex_color;
in vec3 world_position;

uniform vec4 model_color = vec4(1.0, 1.0, 1.0, 1.0);
uniform vec3 light_direction = vec3(0.0, 1.0, 0.5);
uniform float ambient_strength = 0.3;
uniform float diffuse_strength = 0.7;

out vec4 fragment_color;

void main() {
    // Use vertex color if available, otherwise use uniform color
    vec4 base_color = vertex_color;
    if (length(vertex_color.rgb) < 0.01) {
        base_color = model_color;
    }
    
    // Simple lighting calculation
    vec3 light_dir = normalize(light_direction);
    float diff = max(dot(normalize(vec3(0.0, 1.0, 0.0)), light_dir), 0.0);
    
    vec3 ambient = ambient_strength * base_color.rgb;
    vec3 diffuse = diffuse_strength * diff * base_color.rgb;
    
    fragment_color = vec4(ambient + diffuse, base_color.a);
}
