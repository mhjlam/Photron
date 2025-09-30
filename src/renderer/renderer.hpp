/**
 * @file renderer.hpp
 * @brief OpenGL-based 3D visualization renderer for Monte Carlo simulation results
 *
 * The Renderer class provides real-time 3D visualization of photon transport
 * simulation results using modern OpenGL techniques. It renders voxelized geometry,
 * photon paths, energy deposition, and provides interactive camera controls.
 */

#pragma once

#include <atomic>
#include <concepts>
#include <functional>
#include <future>
#include <memory>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "renderer/camera.hpp"
#include "renderer/settings.hpp"

// Forward declarations
class Simulator;

// Forward declare OpenGL types
using GLuint = unsigned int;
using GLint = int;

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
	Renderer();

	/**
	 * @brief Destroy the Renderer object and clean up OpenGL resources
	 */
	~Renderer();

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

	// Input handling methods

	/**
	 * @brief Handle keyboard input events for camera and rendering control
	 *
	 * @param key GLFW key code
	 * @param scancode Platform-specific scan code
	 * @param action GLFW action (press, release, repeat)
	 * @param mods Modifier key flags (Shift, Ctrl, Alt, etc.)
	 */
	void handle_key_input(int key, int scancode, int action, int mods);

	/**
	 * @brief Handle mouse cursor movement for camera control
	 *
	 * @param xpos Current cursor X coordinate in screen space
	 * @param ypos Current cursor Y coordinate in screen space
	 */
	void handle_mouse_move(float xpos, float ypos);

	/**
	 * @brief Handle mouse button press/release events
	 *
	 * @param button Mouse button code (left, right, middle, etc.)
	 * @param action GLFW action (press or release)
	 * @param mods Modifier key flags
	 */
	void handle_mouse_button(int button, int action, int mods);

	/**
	 * @brief Handle mouse scroll events for camera zoom control
	 *
	 * @param xoffset Horizontal scroll offset (typically 0)
	 * @param yoffset Vertical scroll offset (mouse wheel)
	 */
	void handle_mouse_scroll(float xoffset, float yoffset);

	// Camera control methods

	/**
	 * @brief Reset camera to default position and orientation
	 *
	 * Returns camera to initial view position with appropriate zoom
	 * to frame the simulation geometry.
	 */
	void reset_camera();

	/**
	 * @brief Set camera interaction mode
	 *
	 * Switches between arc-ball camera (orbits around target) and
	 * free-flight camera (first-person style movement).
	 *
	 * @param is_arc_mode True for arc-ball mode, false for free-flight mode
	 */
	void set_camera_mode(bool is_arc_mode);

	/**
	 * @brief Check if camera is in arc-ball mode
	 *
	 * @return true if in arc-ball mode, false if in free-flight mode
	 */
	bool is_arc_camera_mode() const { return orbit_camera_mode_; }

	/**
	 * @brief Check if mouse cursor should be captured for camera control
	 *
	 * @return true if mouse should be captured (hidden and centered)
	 */
	bool should_capture_mouse() const;

	/**
	 * @brief Get mutable reference to camera for direct manipulation
	 *
	 * @return Camera& Reference to internal camera object
	 */
	Camera& get_camera() { return camera_; }

	// Callback for camera mode changes - Modern C++20 with perfect forwarding
	/**
	 * @brief Set callback for camera mode change notifications
	 *
	 * Registers a callback function that will be invoked whenever the camera
	 * mode switches between arc-ball and free-flight modes.
	 *
	 * @tparam Callable Function object type that accepts a boolean parameter
	 * @param callback Function to call when camera mode changes (true = arc-ball, false = free-flight)
	 */
	template<typename Callable>
		requires std::invocable<Callable, bool>
	void set_camera_mode_change_callback(Callable&& callback) {
		camera_mode_change_callback_ = std::forward<Callable>(callback);
	}

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
	 * @brief Get current rendering and visualization settings
	 *
	 * @return const Settings& Reference to current settings configuration
	 */
	const Settings& get_settings() const { return settings_; }

	/**
	 * @brief Finalize voxel instance data preparation for rendering
	 *
	 * Processes collected voxel data, applies color mapping based on energy
	 * values, and prepares instanced rendering data for GPU upload.
	 *
	 * @param mode Voxel coloring scheme (Absorption, Emittance, or Layers)
	 */
	void end_voxel_instances(VoxelMode mode = VoxelMode::Absorption);

	/**
	 * @brief Render all prepared voxel instances to screen
	 *
	 * Executes GPU instanced rendering of voxelized geometry using
	 * previously prepared instance data and appropriate shaders.
	 */
	void draw_voxel_instances();

	/**
	 * @brief Add triangle geometry instance for medium boundary visualization
	 *
	 * Queues a triangle for instanced rendering with specified vertices,
	 * color, and surface normal for proper lighting calculations.
	 *
	 * @param v0 First triangle vertex position
	 * @param v1 Second triangle vertex position
	 * @param v2 Third triangle vertex position
	 * @param color Triangle color and opacity (default: opaque white)
	 * @param normal Surface normal vector for lighting (default: +Z axis)
	 */
	void add_triangle_instance(const glm::vec3& v0,
							   const glm::vec3& v1,
							   const glm::vec3& v2,
							   const glm::vec4& color = glm::vec4(1.0f),
							   const glm::vec3& normal = glm::vec3(0.0f, 0.0f, 1.0f));

	/**
	 * @brief Render all queued triangle instances for medium geometry
	 *
	 * Executes instanced rendering of all added triangles for
	 * medium boundary and geometry visualization.
	 */
	void draw_triangle_instances();

	/**
	 * @brief Add line instance for medium wireframe visualization
	 *
	 * Queues a line segment for instanced rendering with gradient
	 * color support for enhanced visual clarity.
	 *
	 * @param start Line start position in world coordinates
	 * @param end Line end position in world coordinates
	 * @param color Line color (applied to both endpoints for solid color)
	 */
	void add_medium_line_instance(const glm::vec3& start,
								  const glm::vec3& end,
								  const glm::vec4& color = glm::vec4(1.0f));

	/**
	 * @brief Render all queued medium line instances
	 *
	 * Executes instanced rendering of all added line segments
	 * for medium wireframe and boundary visualization.
	 */
	void draw_medium_line_instances();

	/**
	 * @brief Render energy value labels as screen-space billboards
	 *
	 * Displays floating text labels showing energy values at key interaction
	 * points, with automatic screen positioning and visibility management.
	 *
	 * @param settings Current rendering settings for label visibility control
	 */
	void draw_labels(const Settings& settings);

	/**
	 * @brief Pre-compute energy labels from simulation data
	 *
	 * Analyzes emitter data and creates cached text labels with proper
	 * classification and positioning for efficient per-frame rendering.
	 */
	void cache_energy_labels();

	/**
	 * @brief Invalidate cached energy labels for regeneration
	 *
	 * Forces recalculation of energy labels when simulation data changes
	 * or display settings are modified.
	 */
	void invalidate_energy_label_cache();

	/**
	 * @brief Update screen positions for 3D world text labels
	 *
	 * Projects 3D world coordinates to screen space for proper billboard
	 * positioning. Called once per frame for optimal performance.
	 */
	void update_energy_label_screen_positions();

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
	 * @brief Invalidate cached photon path rendering instances
	 *
	 * Forces regeneration of line and point instance data for
	 * photon path visualization when path data changes.
	 */
	void invalidate_path_instances_cache();

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
	}

	/**
	 * @brief Generate adaptive color mapping for energy values
	 *
	 * Creates perceptually uniform color representation of energy values
	 * using adaptive scaling and logarithmic mapping for optimal contrast.
	 *
	 * @param energy Raw energy value to color-code
	 * @param min_energy Minimum energy value for color scale
	 * @param max_energy Maximum energy value for color scale
	 * @return glm::vec4 RGBA color value for energy visualization
	 */
	glm::vec4 get_adaptive_energy_color(float energy, float min_energy, float max_energy);

	/**
	 * @brief Generate layer-specific color mapping for multi-material visualization
	 *
	 * Creates consistent color coding for energy values within specific material
	 * layers, maintaining visual distinction between different materials.
	 *
	 * @param energy Energy value within the specified layer
	 * @param min_energy Minimum energy value for layer color scale
	 * @param max_energy Maximum energy value for layer color scale
	 * @param layer_id Unique identifier for material layer
	 * @return glm::vec4 RGBA color value with layer-specific hue coding
	 */
	glm::vec4 layer_energy_color(float energy, float min_energy, float max_energy, uint8_t layer_id);

private:
	void setup_opengl();
	void update_camera();
	void update_camera_target(const Simulator& simulator);

	// OpenGL state management helpers
	void use_shader_program(GLuint program_id) const;
	void enable_blending() const;
	void disable_blending() const;

	// Shader-based drawing functions
	void draw_volume(const Simulator& simulator);
	void draw_voxels(const Settings& settings);
	void draw_paths(const Settings& settings);

	// Utility methods
	bool is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const;

	// Shader management methods
	bool setup_voxel_instanced_rendering();
	bool setup_line_instanced_rendering();
	bool setup_point_instanced_rendering();
	bool setup_triangle_instanced_rendering();
	bool setup_medium_line_vao();

	// Cache energy range calculation
	void update_cached_energy_range(const Settings& settings) const;

	std::string load_shader_source(const std::string& file_path);
	GLuint create_shader_program(const std::string& vertex_source, const std::string& fragment_source);

private:
	// Vertex structures for fallback point rendering
	struct PointVertex
	{
		glm::vec3 position {};
		glm::vec4 color {1.0f};
	};

	// Instance data structure for voxel rendering
	struct VoxelInstance
	{
		glm::vec3 position {};
		glm::vec4 color {1.0f};
		float scale {1.0f};
		float depth {0.0f}; // For depth sorting
	};

	// Instance data structure for line rendering with gradient support
	struct LineInstance
	{
		glm::vec3 start {};
		glm::vec3 end {};
		glm::vec4 start_color {1.0f};
		glm::vec4 end_color {1.0f};
	};

	// Instance data structure for point rendering
	struct PointInstance
	{
		glm::vec3 position {};
		glm::vec4 color {1.0f};
		float size {1.0f};
	};

	// Instance data structure for medium triangle rendering
	struct TriangleInstance
	{
		glm::vec3 v0 {};
		glm::vec3 v1 {};
		glm::vec3 v2 {};
		glm::vec4 color {1.0f};
		glm::vec3 normal {}; // For lighting calculations
	};

	// Instance data structure for medium line rendering
	struct MediumLineInstance
	{
		glm::vec3 start {};
		glm::vec3 end {};
		glm::vec4 start_color {1.0f};
		glm::vec4 end_color {1.0f};
	};

	// Energy label structure for billboard text rendering
	struct EnergyLabel
	{
		glm::vec3 world_position;
		std::string text;
		glm::vec4 color;
		float scale {1.0f};
		glm::vec2 screen_position;          // Cached screen position
		bool screen_position_valid {false}; // Whether screen position is current
	};

	// Key state tracking for smooth movement
	struct KeyState
	{
		bool w_pressed {false};
		bool a_pressed {false};
		bool s_pressed {false};
		bool d_pressed {false};
		bool space_pressed {false};
		bool shift_pressed {false};
	};

private:
	// ========================================
	// OpenGL Rendering Resources
	// ========================================

	GLuint voxel_vao_ {0};                                ///< OpenGL vertex array object for voxel geometry
	GLuint voxel_vbo_ {0};                                ///< OpenGL vertex buffer for voxel vertex data
	GLuint voxel_instance_vbo_ {0};                       ///< OpenGL instance buffer for voxel data
	GLuint voxel_shader_ {0};                             ///< OpenGL shader program for voxel rendering
	std::vector<VoxelInstance> voxel_instances_;          ///< Container for voxel instance data

	GLuint lines_vao_ {0};                                ///< OpenGL vertex array object for line geometry
	GLuint lines_vbo_ {0};                                ///< OpenGL vertex buffer for line vertex data
	GLuint lines_instance_vbo_ {0};                       ///< OpenGL instance buffer for line data
	GLuint lines_shader_ {0};                             ///< OpenGL shader program for line rendering
	std::vector<LineInstance> line_instances_;            ///< Container for line instance data

	GLuint wireframe_vao_ {0};                            ///< OpenGL vertex array object for wireframe geometry
	GLuint wireframe_instance_vbo_ {0};                   ///< OpenGL instance buffer for wireframe data
	std::vector<MediumLineInstance> wireframe_instances_; ///< Container for medium wireframe instance data

	GLuint points_vao_ {0};                               ///< OpenGL vertex array object for point geometry
	GLuint points_vbo_ {0};                               ///< OpenGL vertex buffer for point vertex data
	GLuint points_instance_vbo_ {0};                      ///< OpenGL instance buffer for point data
	GLuint points_shader_ {0};                            ///< OpenGL shader program for point rendering
	std::vector<PointInstance> point_instances_;          ///< Container for point instance data

	GLuint triangles_vao_ {0};                            ///< OpenGL vertex array object for triangle geometry
	GLuint triangles_vbo_ {0};                            ///< OpenGL vertex buffer for triangle vertex data
	GLuint triangles_instance_vbo_ {0};                   ///< OpenGL instance buffer for triangle data
	GLuint triangles_shader_ {0};                         ///< OpenGL shader program for triangle rendering
	std::vector<TriangleInstance> triangle_instances_;    ///< Container for triangle instance data

	// ========================================
	// Core Renderer State
	// ========================================

	Camera camera_;                 ///< 3D camera system with orbital and free-flight modes
	bool orbit_camera_mode_ {true}; ///< Camera mode flag: true for Orbit, false for Free

	Settings settings_;             ///< Current rendering and visualization settings
	int viewport_width_ {800};      ///< Current viewport width in pixels
	int viewport_height_ {600};     ///< Current viewport height in pixels

	std::function<void(const std::string&, float, float, const glm::vec4&)>
		text_render_callback_;      ///< Text rendering callback for ImGui integration

	// ========================================
	// Performance Caching System
	// ========================================

	std::vector<EnergyLabel> cached_energy_labels_;           ///< Cached energy labels for billboard rendering
	bool energy_labels_cached_ {false};                       ///< Flag indicating if energy labels are current

	mutable float cached_min_energy_ {1.0f};                  ///< Cached minimum energy value for color mapping
	mutable float cached_max_energy_ {0.0f};                  ///< Cached maximum energy value for color mapping
	mutable bool energy_range_cached_ {false};                ///< Flag indicating if energy range cache is valid
	mutable size_t last_path_count_ {0};                      ///< Last known photon path count for cache invalidation
	mutable size_t last_voxel_data_version_ {0};              ///< Last known voxel data version for cache invalidation

	mutable float cached_surface_y_ {0.1f};                   ///< Cached surface Y-coordinate for geometry calculations
	mutable bool surface_cached_ {false};                     ///< Flag indicating if surface geometry cache is valid

	mutable std::vector<LineInstance> cached_line_instances_; ///< Cached photon path line instances
	mutable std::vector<PointInstance> cached_point_instances_; ///< Cached photon scattering point instances
	mutable bool path_instances_cached_ {false};                ///< Flag indicating if photon path instances are cached
	mutable size_t cached_photon_count_ {0};                    ///< Number of photons already processed in cache
	mutable uint64_t last_simulation_version_ {0};              ///< Last simulation version for cache invalidation

	mutable std::vector<TriangleInstance> cached_triangle_instances_;      ///< Cached medium geometry triangles
	mutable std::vector<MediumLineInstance> cached_medium_line_instances_; ///< Cached medium geometry wireframes
	mutable bool medium_geometry_cached_ {false};     ///< Flag indicating if medium geometry cache is valid
	mutable size_t last_medium_geometry_version_ {0}; ///< Last medium geometry version for cache invalidation

	// ========================================
	// GPU Buffer State Tracking
	// ========================================

	mutable bool line_buffer_uploaded_ {false};           ///< Flag indicating if line buffer is uploaded to GPU
	mutable bool point_buffer_uploaded_ {false};          ///< Flag indicating if point buffer is uploaded to GPU
	mutable bool voxel_buffer_uploaded_ {false};          ///< Flag indicating if voxel buffer is uploaded to GPU
	mutable bool voxel_instances_dirty_ {true};           ///< Flag to force voxel data recalculation

	mutable bool point_geometry_buffer_uploaded_ {false}; ///< Flag indicating if dynamic point geometry is uploaded
	mutable bool triangle_instances_uploaded_ {false};    ///< Flag indicating if triangle instances are uploaded
	mutable bool medium_line_instances_uploaded_ {false}; ///< Flag indicating if medium line instances are uploaded

	// ========================================
	// OpenGL State Management
	// ========================================

	mutable GLuint current_shader_program_ {0};       ///< Currently bound shader program ID (for state caching)
	mutable bool blend_enabled_ {false};              ///< Current OpenGL blend state (for state caching)

	mutable VoxelMode current_voxel_mode_ {
		VoxelMode::Absorption};                       ///< Current voxel rendering mode for transparency sorting

	mutable glm::vec3 cached_camera_position_ {0.0f}; ///< Cached camera position for change detection
	mutable glm::vec3 cached_camera_target_ {0.0f, 0.0f, -1.0f}; ///< Cached camera target for change detection

	// ========================================
	// Asynchronous Processing
	// ========================================

	mutable std::vector<VoxelInstance> background_sorted_voxels_;   ///< Background-sorted voxels for smooth rendering
	mutable std::future<void> sorting_future_;                      ///< Future for asynchronous voxel sorting
	mutable std::atomic<bool> background_sort_ready_ {false};       ///< Flag indicating if background sort is complete
	mutable std::atomic<bool> background_sort_in_progress_ {false}; ///< Flag indicating if background sort is running

	// ========================================
	// Shader Uniform Locations (Cached)
	// ========================================

	mutable GLint point_mvp_uniform_location_ {-1};  ///< Cached MVP uniform location for fallback point rendering
	mutable GLint voxel_mvp_uniform_location_ {-1};  ///< Cached MVP uniform location for instanced voxels
	mutable GLint point_size_uniform_location_ {-1}; ///< Cached point size uniform location
	mutable GLint line_instanced_mvp_uniform_location_ {-1};  ///< Cached MVP uniform location for instanced lines
	mutable GLint point_instanced_mvp_uniform_location_ {-1}; ///< Cached MVP uniform location for instanced points

	// ========================================
	// Input and Interaction State
	// ========================================

	bool camera_state_changed_ {true};                      ///< Flag indicating if camera state has changed

	std::function<void(bool)> camera_mode_change_callback_; ///< Callback for camera mode changes

	KeyState key_state_;                                    ///< Current keyboard input state for smooth movement

	Simulator* simulator_ {nullptr};                        ///< Reference to current simulator instance
};
