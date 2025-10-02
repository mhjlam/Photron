/**
 * @file renderer.cpp
 * @brief Implementation of OpenGL-based 3D visualization system
 *
 * Contains the complete implementation of the Renderer class for real-time
 * visualization of Monte Carlo photon transport results, including voxel
 * rendering, photon path display, and interactive camera controls.
 */

#include "renderer.hpp"

#include <algorithm>
#include <array>
#include <format>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "common/error_types.hpp"
#include "common/result.hpp"
#include "math/range.hpp"
#include "renderer/camera.hpp"
#include "renderer/settings.hpp"
#include "renderer/shader.hpp"
#include "simulator/config.hpp"
#include "simulator/layer.hpp"
#include "simulator/medium.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"
#include "simulator/voxel.hpp"

bool Renderer::initialize() {
	// Log initialization progress when logging is enabled
	if (Config::is_initialized() && Config::get().log()) {
		std::cout << "Initializing OpenGL renderer..." << std::endl;
	}

	// Configure OpenGL state
	setup_opengl();

	// Initialize VoxelRenderer (handles all voxel rendering)
	if (!voxel_renderer.initialize()) {
		std::cerr << "Failed to initialize VoxelRenderer" << std::endl;
		return false;
	}

	// Initialize PathRenderer (handles all photon path rendering)
	if (!path_renderer.initialize()) {
		std::cerr << "Failed to initialize PathRenderer" << std::endl;
		return false;
	}

	// Initialize GeometryRenderer (handles all medium geometry rendering)
	if (!geometry_renderer.initialize()) {
		std::cerr << "Failed to initialize GeometryRenderer" << std::endl;
		return false;
	}

	// Initialize LabelRenderer (handles all energy label rendering)
	if (!label_renderer.initialize()) {
		std::cerr << "Failed to initialize LabelRenderer" << std::endl;
		return false;
	}

	// Initialize camera with default orbital mode
	camera.set_fps_mode(!orbit_camera_mode);

	// Setup initial camera
	update_camera();

	if (Config::is_initialized() && Config::get().log()) {
		std::cout << "Renderer initialized successfully" << std::endl;
	}
	return true;
}

void Renderer::setup_opengl() {
	// Query OpenGL version
	const GLubyte* version = glGetString(GL_VERSION);
	const GLubyte* renderer = glGetString(GL_RENDERER);

	if (Config::is_initialized() && Config::get().log()) {
		std::cout << "OpenGL Version:  " << version << std::endl;
		std::cout << "Graphics Device: " << renderer << std::endl;
	}

	// Enable depth testing with polygon offset to reduce Z-fighting
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// Enable polygon offset to reduce Z-fighting between adjacent voxels
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 1.0f);

	// Enable blending for transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Enable point size control from shaders
	glEnable(GL_PROGRAM_POINT_SIZE);

	// Enable line width control - use thicker lines for better visibility like original
	glLineWidth(3.0f); // Much thicker lines to match backup

	// Enable line anti-aliasing for smooth rendering
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	// Set clear color to dark charcoal gray - easy on the eyes and neutral
	// Avoids conflicts with layer colors (red, blue, green, magenta, yellow, cyan)
	glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
}

void Renderer::update() {
	update_camera();
}

void Renderer::render(Simulator& simulator) {
	// Store simulator reference for rendering methods
	simulator_ = std::ref(simulator);

	// Invalidate caches only when simulation changes
	uint64_t current_sim_version = simulator.get_simulation_version();
	if (current_sim_version != last_simulation_version_) {
		// Simulation data has changed - invalidate rendering caches
		// Includes both energy visualization and photon path instances
		invalidate_all_caches();
		last_simulation_version_ = current_sim_version;
	}

	// Initialize camera targeting on first frame
	static bool first_frame = true;
	if (first_frame) {
		update_camera_target(simulator);
		first_frame = false;
	}

	// Clear frame buffers for new frame
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Update camera matrices and transformations
	update_camera();

	// Render geometry boundaries and wireframes
	geometry_renderer.render(simulator, settings_, camera);

	// Clear depth buffer for proper layering of volume elements
	glClear(GL_DEPTH_BUFFER_BIT);

	// Render voxelized absorption data using instanced rendering
	if (simulator_) {
		voxel_renderer.render_voxels(simulator_->get(), settings_, camera);
	}

	// Render photon transport paths with proper depth testing
	if (settings_.draw_paths && simulator_) {
		path_renderer.render_paths(settings_, simulator_->get(), camera);
	}

	// Render energy labels as screen-space billboards
	label_renderer.render(simulator, settings_, camera);
}

void Renderer::update_camera() {
	// Update camera matrices for viewport
	float aspect = static_cast<float>(viewport_width_) / static_cast<float>(viewport_height_);
	camera.set_aspect_ratio(aspect);
}

void Renderer::update_camera_target(const Simulator& simulator) {
	Range3 bounds = simulator.get_combined_bounds();

	// Set camera target to center of bounds
	glm::vec3 center = to_float(bounds.center());
	camera.set_target(center);

	// Set material elevation bounds for scroll wheel constraint (Y-axis)
	camera.set_elevation_bounds(static_cast<float>(bounds.min_bounds.y), static_cast<float>(bounds.max_bounds.y));
}

void Renderer::set_viewport(int width, int height) {
	viewport_width_ = width;
	viewport_height_ = height;
	glViewport(0, 0, width, height);

	// Update VoxelRenderer viewport
	voxel_renderer.set_viewport(width, height);

	// Update LabelRenderer viewport
	label_renderer.set_viewport(width, height);

	// Invalidate energy label cache when viewport changes to ensure proper coordinate conversion
	label_renderer.invalidate_cache();
}

void Renderer::auto_manage_energy_labels(Settings& settings) {
	if (!simulator_) {
		return;
	}
	label_renderer.auto_manage_labels(settings, simulator_->get().photons.size());
}

glm::vec2 Renderer::world_to_screen(const glm::vec3& world_pos, int screen_width, int screen_height) const {
	// Use the actual camera matrices that are being used for rendering
	glm::mat4 mvp = camera.get_mvp_matrix();

	// Transform world position to clip space
	glm::vec4 clip_pos = mvp * glm::vec4(world_pos, 1.0f);

	// Check if point is behind the camera or at infinity
	if (clip_pos.w <= 0.0f) {
		return glm::vec2(-1.0f, -1.0f); // Invalid point (behind camera)
	}

	// Perform perspective divide to get normalized device coordinates
	glm::vec3 ndc = glm::vec3(clip_pos) / clip_pos.w;

	// Check if point is outside the view frustum
	if (ndc.x < -1.0f || ndc.x > 1.0f || ndc.y < -1.0f || ndc.y > 1.0f || ndc.z < -1.0f || ndc.z > 1.0f) {
		return glm::vec2(-1.0f, -1.0f); // Invalid point (outside frustum)
	}

	// Convert from NDC [-1, 1] to screen coordinates [0, screen_size]
	float screen_x = (ndc.x + 1.0f) * 0.5f * screen_width;
	float screen_y = (1.0f - ndc.y) * 0.5f * screen_height; // Flip Y for screen coordinates

	return glm::vec2(screen_x, screen_y);
}

// Comprehensive cache invalidation
void Renderer::invalidate_all_caches() {
	// Invalidate label renderer cache
	label_renderer.invalidate_cache();

	// Invalidate path instance caches
	path_renderer.invalidate_cache();

	// Invalidate specialized renderer caches
	voxel_renderer.invalidate_cache();
	geometry_renderer.invalidate_cache();
}

// Only invalidate path cache when path-related settings change
void Renderer::set_settings(const Settings& new_settings) {
	// Check if path-related settings have changed
	bool path_settings_changed = (settings_.draw_paths != new_settings.draw_paths);

	// Update settings
	settings_ = new_settings;

	// Only invalidate path cache if path visualization settings changed
	if (path_settings_changed) {
		path_renderer.invalidate_cache();
	}
}
