/**
 * @file voxel_renderer.hpp
 * @brief Complete voxel rendering system for energy visualization
 *
 * Handles ALL aspects of voxel rendering including state management, OpenGL setup,
 * instanced rendering, energy-based coloring, transparency sorting, caching,
 * and performance optimizations. Supports multiple visualization modes including
 * absorption, emittance, and layer-based coloring.
 */

#pragma once

#include <atomic>
#include <future>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "renderer/settings.hpp"

// Forward declarations
class Simulator;
class Camera;
class Voxel;
enum class VoxelMode;

/**
 * @brief Instance data structure for voxel rendering
 */
struct VoxelInstance
{
	glm::vec3 position {};
	glm::vec4 color {1.0f};
	float scale {1.0f};
	float depth {0.0f}; // For depth sorting
};

/**
 * @brief Complete voxel rendering system with full OpenGL state management
 *
 * VoxelRenderer is a complete, self-contained rendering system that handles
 * ALL voxel-related rendering operations. It manages its own OpenGL resources,
 * caching, background processing, and GPU operations.
 */
class VoxelRenderer
{
public:
	/**
	 * @brief Construct VoxelRenderer and initialize OpenGL resources
	 */
	VoxelRenderer();

	/**
	 * @brief Destroy VoxelRenderer and clean up OpenGL resources
	 */
	~VoxelRenderer();

	/**
	 * @brief Initialize OpenGL resources (shaders, VAO, VBOs)
	 * @return true if initialization succeeded
	 */
	bool initialize();

	/**
	 * @brief Complete voxel rendering pipeline
	 * @param simulator Simulator containing voxel data
	 * @param settings Current rendering settings
	 * @param camera Camera for MVP matrix and depth calculations
	 */
	void render_voxels(const Simulator& simulator, const Settings& settings, const Camera& camera);

	/**
	 * @brief Invalidate cached data when simulation changes
	 */
	void invalidate_cache();

	/**
	 * @brief Update viewport dimensions for window resize events
	 */
	void set_viewport(int width, int height);

private:
	// OpenGL Resources
	GLuint voxel_shader_ = 0;
	GLuint voxel_vao_ = 0;
	GLuint voxel_vbo_ = 0;
	GLuint voxel_instance_vbo_ = 0;
	GLuint voxel_mvp_uniform_location_ = 0;

	// Rendering State
	std::vector<VoxelInstance> voxel_instances_;
	mutable std::vector<VoxelInstance> background_sorted_voxels_;
	mutable std::atomic<bool> background_sort_in_progress_ {false};
	mutable std::atomic<bool> background_sort_ready_ {false};
	mutable std::future<void> sorting_future_;

	// Caching and Performance
	bool voxel_instances_dirty_ = true;
	bool voxel_buffer_uploaded_ = false;
	VoxelMode current_voxel_mode_ = VoxelMode::Layers;
	glm::vec3 cached_camera_position_ {0.0f};
	glm::vec3 cached_camera_target_ {0.0f};

	// Energy Range Caching
	mutable bool energy_range_cached_ = false;
	mutable VoxelMode cached_range_mode_ = VoxelMode::Layers;
	mutable float cached_min_energy_ = 0.0f;
	mutable float cached_max_energy_ = 1.0f;

	// Global Energy Total Caching
	mutable bool global_energy_cached_ = false;
	mutable VoxelMode cached_energy_mode_ = VoxelMode::Layers;
	mutable double cached_global_absorption_total_ = 0.0;
	mutable double cached_global_emittance_total_ = 0.0;
	mutable double cached_global_combined_total_ = 0.0;

	// Viewport
	int viewport_width_ = 800;
	int viewport_height_ = 600;

	// Internal Methods
	void collect_voxel_instances(const Simulator& simulator, const Settings& settings, const Camera& camera);
	void update_cached_energy_range(const Simulator& simulator, const Settings& settings) const;
	void update_cached_global_energy(const Simulator& simulator, VoxelMode mode) const;
	void end_voxel_instances(VoxelMode mode, const Camera& camera);
	void draw_voxel_instances(const Camera& camera);
	bool setup_voxel_rendering();
	bool setup_voxel_geometry();
	glm::vec4 layer_energy_color(float energy, float min_energy, float max_energy, uint8_t layer_id) const;

	// Calculate voxel energy as percentage of total medium energy for proper visualization
	float calculate_voxel_energy_percentage(const Voxel* voxel, VoxelMode mode, const Simulator& simulator) const;
};
