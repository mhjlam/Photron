#version 450 core

layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec4 aColor;

uniform mat4 uMVP;
uniform vec4 uClipPlanes[6]; // Support up to 6 clipping planes
uniform int uNumClipPlanes;  // Number of active clipping planes

out vec4 vColor;
out float gl_ClipDistance[6]; // Built-in clipping distances

void main() {
    gl_Position = uMVP * vec4(aPosition, 1.0);
    vColor = aColor;
    
    // Calculate clip distances for each active plane
    for (int i = 0; i < uNumClipPlanes && i < 6; i++) {
        // Plane equation: ax + by + cz + d = 0
        // Distance from point to plane: ax + by + cz + d
        gl_ClipDistance[i] = dot(aPosition, uClipPlanes[i].xyz) + uClipPlanes[i].w;
    }
    
    // Disable unused clip planes
    for (int i = uNumClipPlanes; i < 6; i++) {
        gl_ClipDistance[i] = 1.0; // Always positive = not clipped
    }
}
