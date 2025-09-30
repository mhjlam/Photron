#version 450 core

// Base triangle vertex (local coordinates)
layout (location = 0) in vec3 aPosition;

// Per-instance attributes  
layout (location = 1) in vec3 aInstanceV0;
layout (location = 2) in vec3 aInstanceV1;
layout (location = 3) in vec3 aInstanceV2;
layout (location = 4) in vec4 aInstanceColor;
layout (location = 5) in vec3 aInstanceNormal;

// Uniforms
uniform mat4 uMVP;

// Output to fragment shader
out vec4 vertexColor;
out vec3 worldPos;
out vec3 worldNormal;

void main() {
    // Transform base vertex using barycentric coordinates
    vec3 worldPosition;
    
    // Map base triangle (0,0,0), (1,0,0), (0,1,0) to instance triangle
    if (aPosition.x == 0.0 && aPosition.y == 0.0) {
        // First vertex maps to v0
        worldPosition = aInstanceV0;
    }
    else if (aPosition.x == 1.0 && aPosition.y == 0.0) {
        // Second vertex maps to v1
        worldPosition = aInstanceV1;
    }
    else {
        // Third vertex maps to v2
        worldPosition = aInstanceV2;
    }
    
    // Transform to clip space
    gl_Position = uMVP * vec4(worldPosition, 1.0);
    
    // Pass through instance data
    vertexColor = aInstanceColor;
    worldPos = worldPosition;
    worldNormal = aInstanceNormal;
}
