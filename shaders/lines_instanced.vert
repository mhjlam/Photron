#version 330 core

// Vertex attributes (line geometry - just two endpoints)
layout(location = 0) in vec3 aPosition;

// Instance attributes (per line segment)
layout(location = 1) in vec3 aInstanceStart;
layout(location = 2) in vec3 aInstanceEnd;
layout(location = 3) in vec4 aInstanceStartColor;
layout(location = 4) in vec4 aInstanceEndColor;

// Uniforms
uniform mat4 uMVP;

// Outputs to fragment shader
out vec4 vColor;

void main()
{
    // Transform line vertex based on instance data
    // aPosition.x determines which endpoint: 0 = start, 1 = end
    vec3 worldPos = mix(aInstanceStart, aInstanceEnd, aPosition.x);
    gl_Position = uMVP * vec4(worldPos, 1.0);
    
    // Interpolate color based on position along line
    vColor = mix(aInstanceStartColor, aInstanceEndColor, aPosition.x);
}
