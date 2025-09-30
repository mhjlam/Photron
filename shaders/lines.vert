#version 450 core

// Base line vertex (local coordinates)
layout(location = 0) in vec3 aPosition;

// Per-instance attributes (supports both single color and gradient lines)
layout(location = 1) in vec3 aInstanceStart;
layout(location = 2) in vec3 aInstanceEnd;
layout(location = 3) in vec4 aInstanceStartColor;
layout(location = 4) in vec4 aInstanceEndColor;  // For gradient support; same as start color for solid lines

// Uniforms
uniform mat4 uMVP;

// Output to fragment shader
out vec4 vertexColor;

void main() {
    // Transform line vertex based on instance data
    // aPosition.x determines which endpoint: 0 = start, 1 = end
    vec3 worldPosition = mix(aInstanceStart, aInstanceEnd, aPosition.x);
    gl_Position = uMVP * vec4(worldPosition, 1.0);
    
    // Interpolate color based on position along line
    // For solid color lines, both start and end colors will be the same
    vertexColor = mix(aInstanceStartColor, aInstanceEndColor, aPosition.x);
}
