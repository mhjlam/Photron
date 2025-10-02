/**
 * @file vertex.hpp
 * @brief Vertex data structures for OpenGL rendering pipeline
 *
 * Defines vertex formats and data structures used throughout the rendering
 * system for efficient GPU data transfer and visualization.
 */

#pragma once

#include <glm/glm.hpp>

/**
 * @struct RenderVertex
 * @brief Standard vertex data structure for OpenGL rendering pipeline
 *
 * Contains all common vertex attributes needed for 3D rendering including
 * position, normal, texture coordinates, and color. Uses GLM types for
 * efficient GPU transfer and mathematical operations.
 */
struct RenderVertex
{
	glm::vec3 position {0.0f};           ///< 3D world space position (x, y, z)
	glm::vec3 normal {0.0f, 1.0f, 0.0f}; ///< Surface normal vector (normalized)
	glm::vec2 tex_coords {0.0f};         ///< Texture coordinates (u, v)
	glm::vec4 color {1.0f};              ///< RGBA color (0.0-1.0 range)

	/**
	 * @brief Default constructor using default member initialization
	 */
	RenderVertex() = default;

	/**
	 * @brief Construct vertex with position only
	 * @param pos 3D position vector
	 */
	explicit RenderVertex(const glm::vec3& pos) noexcept : position(pos) {}

	/**
	 * @brief Construct vertex with position and color
	 * @param pos 3D position vector
	 * @param col RGBA color vector
	 */
	RenderVertex(const glm::vec3& pos, const glm::vec4& col) noexcept : position(pos), color(col) {}

	/**
	 * @brief Construct vertex with all attributes
	 * @param pos 3D position vector
	 * @param norm Surface normal vector
	 * @param tex Texture coordinates
	 * @param col RGBA color vector
	 */
	RenderVertex(const glm::vec3& pos, const glm::vec3& norm, const glm::vec2& tex, const glm::vec4& col) noexcept :
		position(pos), normal(norm), tex_coords(tex), color(col) {}
};
