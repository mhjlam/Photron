#version 330 core

// Inputs from vertex shader
in vec4 vertexColor;

// Output
out vec4 FragColor;

void main()
{
    // Use raw color to match original renderer
    FragColor = vertexColor;
    
    // Discard very transparent fragments to improve performance
    if (FragColor.a < 0.01) {
        discard;
    }
}
