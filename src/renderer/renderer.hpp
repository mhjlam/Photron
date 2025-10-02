/**
 * @file renderer.hpp
 * @brief OpenGL-based 3D visualization renderer for Monte Carlo simulation results
 *
 * The Renderer class provides real-time 3D visualization of photon transport
 * simulation results using modern OpenGL techniques. It renders voxelized geometry,
 * photon paths, energy deposition, and provides interactive camera controls through
 * coordinated specialized rendering components.
 */

#pragma once

#include <concepts>
#include <functional>
#include <optional>
#include <string>

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "renderer/camera.hpp"
#include "renderer/geometry_renderer.hpp"
#include "renderer/label_renderer.hpp"
#include "renderer/path_renderer.hpp"
#include "renderer/settings.hpp"
#include "renderer/shader.hpp"
#include "renderer/voxel_renderer.hpp"

// Forward declarations
class Simulator;

/**
 * @class Renderer
 * @brief OpenGL-based 3D visualization system for Monte Carlo photon transport results
 *
 * The Renderer class provides comprehensive 3D visualization capabilities for
 * Monte Carlo photon transport simulation results. Key features include:
 *
 * - **Voxel Rendering**: Visualizes 3D voxelized geometry with energy-based color mapping
 * - **Photon Path Visualization**: Renders photon trajectories and scattering events
 * - **Energy Deposition Display**: Color-coded visualization of absorbed energy
 * - **Interactive Camera**: Arc-ball and free-flight camera modes with smooth controls
 * - **Performance Optimization**: Instanced rendering, frustum culling, and caching
 * - **Real-time Updates**: Dynamic visualization of ongoing simulations
 *
 * The renderer uses modern OpenGL with shader-based pipeline for efficient
 * rendering of potentially millions of voxels and photon path segments.
 */
class Renderer
{
public:
	/**
	 * @brief Construct a new Renderer object with default settings
	 */
	Renderer() : voxel_renderer(), path_renderer(), geometry_renderer(), label_renderer() {}

	/**
	 * @brief Destroy the Renderer object and clean up OpenGL resources
	 */
	~Renderer() = default;

	/**
	 * @brief Initialize OpenGL context and prepare rendering resources
	 *
	 * Sets up shaders, vertex buffers, and OpenGL state for rendering.
	 * Must be called after valid OpenGL context is established.
	 *
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize();

	/**
	 * @brief Update per-frame renderer state
	 *
	 * Performs any necessary updates before rendering, such as
	 * camera matrix calculations and uniform buffer updates.
	 */
	void update();

	/**
	 * @brief Render complete 3D scene with simulation results
	 *
	 * Main rendering entry point that coordinates all drawing operations
	 * including voxel geometry, photon paths, and energy visualization.
	 *
	 * @param simulator Reference to simulator containing data to visualize
	 */
	void render(Simulator& simulator);

	/**
	 * @brief Update viewport dimensions for window resize events
	 *
	 * @param width New viewport width in pixels
	 * @param height New viewport height in pixels
	 */
	void set_viewport(int width, int height);

	/**
	 * @brief Get mutable reference to camera for direct manipulation
	 *
	 * @return Camera& Reference to internal camera object
	 */
	Camera& get_camera() { return camera; }

	/**
	 * @brief Update rendering and visualization settings
	 *
	 * Applies new rendering configuration including visualization flags,
	 * camera mode, and voxel display options.
	 *
	 * @param settings New settings configuration to apply
	 */
	void set_settings(const Settings& settings);

	/**
	 * @brief Automatically manage energy label visibility based on photon count
	 *
	 * Disables energy labels when photon count exceeds threshold to prevent
	 * visual clutter while maintaining interactive control.
	 *
	 * @param settings Settings object to modify for label visibility
	 */
	void auto_manage_energy_labels(Settings& settings);

	/**
	 * @brief Invalidate all cached rendering data
	 *
	 * Forces regeneration of all cached geometry, energy ranges, and
	 * visualization data when simulation parameters change.
	 */
	void invalidate_all_caches();

	/**
	 * @brief Convert 3D world coordinates to 2D screen coordinates
	 *
	 * Projects world-space position through view and projection matrices
	 * to determine screen-space pixel coordinates for UI positioning.
	 *
	 * @param world_pos 3D position in world coordinate system
	 * @param screen_width Current viewport width in pixels
	 * @param screen_height Current viewport height in pixels
	 * @return glm::vec2 Screen coordinates (pixels from top-left corner)
	 */
	glm::vec2 world_to_screen(const glm::vec3& world_pos, int screen_width, int screen_height) const;

	/**
	 * @brief Set callback function for text rendering integration
	 *
	 * Registers a callback for rendering text labels through external systems
	 * like ImGui, allowing flexible text display without OpenGL dependencies.
	 *
	 * @tparam Callable Function object type for text rendering
	 * @param callback Function accepting (text, x, y, color) parameters for text display
	 */
	template<typename Callable>
		requires std::invocable<Callable, const std::string&, float, float, const glm::vec4&>
	void set_text_render_callback(Callable&& callback) {
		text_render_callback_ = std::forward<Callable>(callback);
		label_renderer.set_text_render_callback(std::forward<Callable>(callback));
	}

	// ========================================
	// Core Renderer State (Public for InputHandler access)
	// ========================================

	Camera camera;                 ///< 3D camera system with orbital and free-flight modes (public for InputHandler)
	bool orbit_camera_mode {true}; ///< Camera mode flag: true for Orbit, false for Free (public for InputHandler)

private:
	// ========================================
	// Core Renderer Implementation
	// ========================================

	void setup_opengl();
	void update_camera();
	void update_camera_target(const Simulator& simulator);

	Settings settings_; ///< Current rendering and visualization settings

	// These renderers handle specific aspects of the visualization pipeline.
	// Made mutable to allow const methods to update rendering state.

	mutable VoxelRenderer voxel_renderer;                        ///< Voxel rendering and energy visualization
	mutable PathRenderer path_renderer;                          ///< Photon path visualization and line rendering
	mutable GeometryRenderer geometry_renderer;                  ///< Medium geometry and wireframe rendering
	mutable LabelRenderer label_renderer;                        ///< Energy label rendering and text overlays

	int viewport_width_ {800};                                   ///< Current viewport width in pixels
	int viewport_height_ {600};                                  ///< Current viewport height in pixels

	std::function<void(const std::string&, float, float, const glm::vec4&)>
		text_render_callback_;                                   ///< Text rendering callback for ImGui integration

	mutable uint64_t last_simulation_version_ {0};               ///< Last simulation version for cache invalidation
	std::optional<std::reference_wrapper<Simulator>> simulator_; ///< Non-owning reference to current simulator instance
};
