#version 450 core

// Input from vertex shader
in vec4 vertexColor;
in vec3 worldPos;
in vec3 worldNormal;

// Output color
out vec4 FragColor;

void main() {
    // Basic lighting calculation
    vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
    vec3 normal = normalize(worldNormal);
    
    // Simple diffuse lighting
    float lightIntensity = max(dot(normal, lightDir), 0.2); // Minimum ambient
    
    // Apply lighting to color, preserving alpha for transparency
    vec3 litColor = vertexColor.rgb * lightIntensity;
    FragColor = vec4(litColor, vertexColor.a);
}
