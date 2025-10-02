#pragma once

#include <glm/glm.hpp>

/**
 * @file path_instance.hpp
 * @brief Shared instance data structures for photon path rendering
 *
 * This file defines the instance data structures used for photon path
 * visualization, shared between Renderer and PathRenderer to maintain
 * type compatibility and prevent circular dependencies.
 */

/**
 * @struct LineInstance
 * @brief Instance data structure for line rendering with gradient support
 *
 * Contains all necessary data for rendering a single line segment with
 * gradient coloring between start and end points.
 */
struct LineInstance
{
	glm::vec3 start {};           ///< Line start position in world coordinates
	glm::vec3 end {};             ///< Line end position in world coordinates
	glm::vec4 start_color {1.0f}; ///< Color at line start (RGBA)
	glm::vec4 end_color {1.0f};   ///< Color at line end (RGBA)
};

/**
 * @struct PointInstance
 * @brief Instance data structure for point rendering
 *
 * Contains all necessary data for rendering a single point with
 * position, color, and size information.
 */
struct PointInstance
{
	glm::vec3 position {};  ///< Point position in world coordinates
	glm::vec4 color {1.0f}; ///< Point color (RGBA)
	float size {1.0f};      ///< Point size in pixels
};
