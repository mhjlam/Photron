#version 450 core

in vec4 vColor;
out vec4 FragColor;

void main() {
    // Make points circular instead of square
    vec2 coord = gl_PointCoord - vec2(0.5);
    if (length(coord) > 0.5) {
        discard;
    }
    
    // Add some shading for better visibility
    float dist = length(coord) * 2.0;
    float alpha = 1.0 - (dist * dist * 0.3);
    
    FragColor = vec4(vColor.rgb, vColor.a * alpha);
}
