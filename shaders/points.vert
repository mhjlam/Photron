#version 330 core

// Vertex attributes
layout (location = 0) in vec3 aPosition;        // Point position (always at origin)

// Instance attributes
layout (location = 1) in vec3 aInstancePos;     // Instance position
layout (location = 2) in vec4 aInstanceColor;   // Instance color
layout (location = 3) in float aInstanceSize;   // Instance size

// Uniforms
uniform mat4 uMVP;

// Outputs to fragment shader
out vec4 vertexColor;

void main() {
    // Transform the point position by the instance position
    vec3 worldPos = aPosition + aInstancePos;
    
    // Apply MVP transformation
    gl_Position = uMVP * vec4(worldPos, 1.0);
    
    // Set point size based on instance size
    gl_PointSize = aInstanceSize;
    
    // Pass color to fragment shader
    vertexColor = aInstanceColor;
}
