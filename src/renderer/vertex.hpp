#pragma once

#include <glm/glm.hpp>

/**
 * @struct RenderVertex
 * @brief Represents a single vertex with position, normal, texture coordinates, and color.
 */
struct RenderVertex
{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 tex_coords;
    glm::vec4 color;

    // Constructors
    RenderVertex() : position(0.0f), normal(0.0f, 1.0f, 0.0f), tex_coords(0.0f), color(1.0f) {}
    
    explicit RenderVertex(const glm::vec3& pos)
        : position(pos), normal(0.0f, 1.0f, 0.0f), tex_coords(0.0f), color(1.0f) {}
    
    RenderVertex(const glm::vec3& pos, const glm::vec4& col)
        : position(pos), normal(0.0f, 1.0f, 0.0f), tex_coords(0.0f), color(col) {}
    
    RenderVertex(const glm::vec3& pos, const glm::vec3& norm, const glm::vec2& tex, const glm::vec4& col)
        : position(pos), normal(norm), tex_coords(tex), color(col) {}
};
