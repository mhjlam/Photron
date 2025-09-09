#version 330 core

// Inputs from vertex shader
in vec4 vColor;

// Output
out vec4 FragColor;

void main()
{
    // Use raw color for photon paths
    FragColor = vColor;
}
