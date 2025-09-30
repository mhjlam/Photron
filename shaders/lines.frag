#version 450 core

// Input from vertex shader
in vec4 vertexColor;

// Output color
out vec4 FragColor;

void main() {
    // Simple passthrough with line color
    FragColor = vertexColor;
}
