/**
 * @file path_renderer.hpp
 * @brief Specialized photon path visualization system
 *
 * Handles rendering of photon transport paths, scattering points, emitter vectors,
 * and energy flow visualization. Uses instanced rendering for performance with
 * large numbers of photons.
 */

#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "path_instance.hpp"
#include "renderer/settings.hpp"
#include "simulator/layer.hpp"

// Forward declarations
class Simulator;
class Camera;

/**
 * @class PathRenderer
 * @brief Self-contained photon path rendering system
 *
 * Handles all aspects of photon path visualization including line segments,
 * scatter points, emitter vectors, and surface interactions. Maintains its
 * own OpenGL resources and sophisticated caching system for optimal performance.
 *
 * Features:
 * - Incremental photon path caching for real-time performance
 * - Energy-based adaptive color mapping with logarithmic scaling
 * - Surface interaction detection and visualization
 * - Emitter vector grouping and deduplication
 * - Specular reflection visualization
 * - Complex path traversal with emitter connections
 * - Independent OpenGL resource management
 */
class PathRenderer
{
public:
	/**
	 * @brief Constructor - initializes OpenGL resources
	 */
	PathRenderer();

	/**
	 * @brief Destructor - cleans up OpenGL resources
	 */
	~PathRenderer();

	// Delete copy constructor and assignment operator
	PathRenderer(const PathRenderer&) = delete;
	PathRenderer& operator=(const PathRenderer&) = delete;

	/**
	 * @brief Initialize OpenGL resources for path rendering
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize();

	/**
	 * @brief Main photon path rendering entry point
	 *
	 * Executes complete photon path visualization pipeline including:
	 * - Incremental cache updates for new photons
	 * - Energy-based color mapping
	 * - Surface interaction calculations
	 * - Emitter vector processing and deduplication
	 * - Line and point instance rendering
	 *
	 * @param settings Current rendering settings for visualization control
	 * @param simulator Reference to simulation data for photon paths
	 * @param camera Reference to camera system for MVP calculations
	 */
	void render_paths(const Settings& settings, const Simulator& simulator, const Camera& camera);

	/**
	 * @brief Force invalidation of all cached path data
	 */
	void invalidate_cache();

private:
	/**
	 * @brief Generate adaptive energy-based color for photon visualization
	 *
	 * Creates enhanced non-linear color mapping optimized for photon energy
	 * visualization with expanded low-energy range and clear high-energy distinction.
	 *
	 * @param energy Current energy value to map
	 * @param min_energy Minimum energy in dataset for normalization
	 * @param max_energy Maximum energy in dataset for normalization
	 * @return glm::vec4 RGBA color vector for energy visualization
	 */
	glm::vec4 get_adaptive_energy_color(float energy, float min_energy, float max_energy) const;
	// ========================================
	// Internal Methods
	// ========================================

	/**
	 * @brief Update cached energy range for color mapping
	 */
	void update_cached_energy_range(const Simulator& simulator);

	/**
	 * @brief Update cached surface Y-coordinate
	 */
	void update_cached_surface(const Simulator& simulator);

	/**
	 * @brief Collect photon path line instances
	 */
	void collect_path_line_instances(const Simulator& simulator,
									 const std::function<glm::vec4(float)>& adaptive_log_color);

	/**
	 * @brief Collect emitter direction vectors
	 */
	void collect_emitter_vectors(const Simulator& simulator, const std::function<glm::vec4(float)>& adaptive_log_color);

	/**
	 * @brief Collect scatter and emitter points
	 */
	void collect_path_point_instances(const Simulator& simulator,
									  const std::function<glm::vec4(float)>& adaptive_log_color);

	/**
	 * @brief Render cached line instances using OpenGL
	 */
	void render_line_instances(const Camera& camera);

	/**
	 * @brief Render cached point instances using OpenGL
	 */
	void render_point_instances(const Camera& camera);

	/**
	 * @brief Setup line instanced rendering VAO and attributes
	 */
	bool setup_line_instanced_rendering();

	/**
	 * @brief Setup point instanced rendering VAO and attributes
	 */
	bool setup_point_instanced_rendering();

	// OpenGL state management
	void enable_blending() const;
	void disable_blending() const;

	// ========================================
	// OpenGL Resources
	// ========================================

	// Line rendering resources
	GLuint lines_vao_ {0};                           ///< VAO for line geometry
	GLuint lines_vbo_ {0};                           ///< VBO for line vertex data
	GLuint lines_instance_vbo_ {0};                  ///< VBO for line instance data
	GLuint lines_shader_ {0};                        ///< Shader program for line rendering
	GLint line_instanced_mvp_uniform_location_ {-1}; ///< Cached MVP uniform location

	// Point rendering resources
	GLuint points_vao_ {0};                           ///< VAO for point geometry
	GLuint points_vbo_ {0};                           ///< VBO for point vertex data
	GLuint points_instance_vbo_ {0};                  ///< VBO for point instance data
	GLuint points_shader_ {0};                        ///< Shader program for point rendering
	GLint point_instanced_mvp_uniform_location_ {-1}; ///< Cached MVP uniform location

	// ========================================
	// Caching System
	// ========================================

	// Instance data caches
	mutable std::vector<LineInstance> cached_line_instances_;   ///< Cached photon path line instances
	mutable std::vector<PointInstance> cached_point_instances_; ///< Cached photon scattering point instances
	mutable bool path_instances_cached_ {false};                ///< Flag indicating if path instances are cached
	mutable size_t cached_photon_count_ {0};                    ///< Number of photons already processed in cache

	// Energy range cache
	mutable float cached_min_energy_ {1.0f};   ///< Cached minimum energy value for color mapping
	mutable float cached_max_energy_ {0.0f};   ///< Cached maximum energy value for color mapping
	mutable bool energy_range_cached_ {false}; ///< Flag indicating if energy range cache is valid

	// Surface geometry cache
	mutable float cached_surface_y_ {0.1f}; ///< Cached surface Y-coordinate for geometry calculations
	mutable bool surface_cached_ {false};   ///< Flag indicating if surface geometry cache is valid

	// Buffer upload tracking
	mutable bool line_buffer_uploaded_ {false};  ///< Flag indicating if line buffer is uploaded to GPU
	mutable bool point_buffer_uploaded_ {false}; ///< Flag indicating if point buffer is uploaded to GPU

	// OpenGL state management
	mutable GLuint current_shader_program_ {0}; ///< Currently bound shader program for state caching
};
