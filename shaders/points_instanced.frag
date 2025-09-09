#version 330 core

// Inputs from vertex shader
in vec4 fragColor;

// Output
out vec4 FragColor;

void main() {
    // Create circular points by discarding fragments outside circle
    vec2 coord = gl_PointCoord - vec2(0.5);
    if (dot(coord, coord) > 0.25) {
        discard;
    }
    
    // Output the color directly without lighting
    FragColor = fragColor;
}
