#version 330 core

// Vertex attributes (cube geometry)
layout(location = 0) in vec3 aPosition;

// Instance attributes (per voxel)
layout(location = 2) in vec3 aInstancePosition;
layout(location = 3) in vec4 aInstanceColor;
layout(location = 4) in float aInstanceScale;

// Uniforms
uniform mat4 uMVP;

// Outputs to fragment shader
out vec4 vColor;

void main()
{
    // Transform vertex position by instance data
    vec3 worldPos = aPosition * aInstanceScale + aInstancePosition;
    gl_Position = uMVP * vec4(worldPos, 1.0);
    
    // Pass color to fragment shader
    vColor = aInstanceColor;
}
