/**
 * @file voxel_renderer.cpp
 * @brief Implementation of complete voxel rendering system
 */

#include "voxel_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ranges>

#include "common/config.hpp"
#include "math/math.hpp"
#include "renderer/camera.hpp"
#include "renderer/shader.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"
#include "simulator/voxel.hpp"

VoxelRenderer::VoxelRenderer() {
	// Constructor - OpenGL resources will be initialized in initialize()
}

VoxelRenderer::~VoxelRenderer() {
	// Clean up OpenGL resources
	if (voxel_shader_) {
		glDeleteProgram(voxel_shader_);
	}
	if (voxel_vao_) {
		glDeleteVertexArrays(1, &voxel_vao_);
	}
	if (voxel_vbo_) {
		glDeleteBuffers(1, &voxel_vbo_);
	}
	if (voxel_instance_vbo_) {
		glDeleteBuffers(1, &voxel_instance_vbo_);
	}

	// Wait for any background sorting to complete
	if (background_sort_in_progress_.load() && sorting_future_.valid()) {
		sorting_future_.wait();
	}
}

bool VoxelRenderer::initialize() {
	return setup_voxel_rendering();
}

void VoxelRenderer::render_voxels(const Simulator& simulator, const Settings& settings, const Camera& camera) {
	// Only draw voxels if voxel rendering is enabled
	if (!settings.draw_voxels) {
		return;
	}

	// Check if we need to rebuild voxel instances
	static VoxelMode last_voxel_mode = VoxelMode::Layers;
	bool need_rebuild = voxel_instances_dirty_ || settings.voxel_mode != last_voxel_mode;

	if (!need_rebuild && voxel_buffer_uploaded_ && !voxel_instances_.empty()) {
		// Use cached voxel instances, but still need to update depth sorting for camera changes
		end_voxel_instances(settings.voxel_mode, camera);
		draw_voxel_instances(camera);
		return;
	}

	// Mark instances as dirty if we're rebuilding due to mode change
	if (settings.voxel_mode != last_voxel_mode) {
		voxel_instances_dirty_ = true;
		voxel_buffer_uploaded_ = false;
		// Invalidate global energy cache when mode changes
		global_energy_cached_ = false;
	}

	last_voxel_mode = settings.voxel_mode;
	current_voxel_mode_ = settings.voxel_mode;

	// Update cached energy range
	update_cached_energy_range(simulator, settings);
	
	// Update cached global energy totals for efficiency
	update_cached_global_energy(simulator, settings.voxel_mode);

	// Set up instanced rendering
	if (voxel_instances_dirty_) {
		voxel_instances_.clear();
	}

	// Collect voxel instances
	collect_voxel_instances(simulator, settings, camera);

	// Process and render all voxels (OpenGL state handled in draw_voxel_instances)
	end_voxel_instances(settings.voxel_mode, camera);
	draw_voxel_instances(camera);
}

void VoxelRenderer::invalidate_cache() {
	voxel_instances_dirty_ = true;
	voxel_buffer_uploaded_ = false;
	energy_range_cached_ = false;
	global_energy_cached_ = false;

	// Wait for any background processing to complete
	if (background_sort_in_progress_.load() && sorting_future_.valid()) {
		sorting_future_.wait();
		background_sort_in_progress_.store(false);
	}
}

void VoxelRenderer::set_viewport(int width, int height) {
	viewport_width_ = width;
	viewport_height_ = height;
}

void VoxelRenderer::collect_voxel_instances(const Simulator& simulator,
											const Settings& settings,
											const Camera& camera) {
	// Get MCML grid parameters
	const auto& config = Config::get();
	const uint32_t nx = static_cast<uint32_t>(config.nx());
	const uint32_t ny = static_cast<uint32_t>(config.ny());
	const uint32_t nz = static_cast<uint32_t>(config.nz());
	const double vox_size = config.vox_size();
	const float voxel_scale = static_cast<float>(vox_size * 0.95f);

	// Pre-calculate camera data for culling
	glm::vec3 camera_pos = camera.get_position();
	glm::vec3 camera_target = camera.get_target();
	glm::vec3 camera_front = glm::normalize(camera_target - camera_pos);
	const float max_render_distance = 50.0f;

	// Get bounds for positioning
	auto bounds = simulator.get_combined_bounds();
	glm::vec3 bounds_min = to_float(bounds.min_bounds);

	// Collect voxel instances
	const uint32_t total_voxels = nx * ny * nz;
	for (uint32_t linear_idx = 0; linear_idx < total_voxels; ++linear_idx) {
		// Convert linear index to 3D coordinates
		const uint32_t iz = linear_idx / (nx * ny);
		const uint32_t iy = (linear_idx % (nx * ny)) / nx;
		const uint32_t ix = linear_idx % nx;

		Voxel* voxel = simulator.voxel_grid(ix, iy, iz);
		if (!voxel || !voxel->material) {
			continue;
		}

		// Calculate world position
		glm::vec3 voxel_pos = bounds_min
							  + glm::vec3(static_cast<float>(vox_size) * ix + static_cast<float>(vox_size) * 0.5f,
										  static_cast<float>(vox_size) * iy + static_cast<float>(vox_size) * 0.5f,
										  static_cast<float>(vox_size) * iz + static_cast<float>(vox_size) * 0.5f);

		// Distance culling
		glm::vec3 camera_offset = voxel_pos - camera_pos;
		float distance_squared = glm::dot(camera_offset, camera_offset);
		if (distance_squared > max_render_distance * max_render_distance) {
			continue;
		}

		// Calculate energy based on mode
		float total_energy;

		// Calculate voxel energy as percentage of total medium energy for proper visualization
		if (settings.voxel_mode == VoxelMode::Layers) {
			total_energy = 1.0f;
		} else {
			total_energy = calculate_voxel_energy_percentage(voxel, settings.voxel_mode, simulator);
		}

		glm::vec4 color(0.0f);

		if (settings.voxel_mode == VoxelMode::Layers) {
			// Layers mode: show all material voxels regardless of energy
			color = layer_energy_color(
				MathConstants::ENERGY_THRESHOLD, cached_min_energy_, cached_max_energy_, voxel->layer_id);
			color.a = 0.05f;
		}
		else {
			// Energy-based rendering using direct percentage mapping
			// total_energy is already the percentage (0.0 to 1.0) of this voxel relative to global total
			
			float dist_from_origin = glm::length(voxel_pos);
			static float cached_max_dist = glm::length(to_float(bounds.max_bounds));
			float normalized_dist = glm::clamp(dist_from_origin / cached_max_dist, 0.0f, 1.0f);

			bool use_emittance = (settings.voxel_mode == VoxelMode::Emittance);

			if (total_energy > MathConstants::ENERGY_THRESHOLD) {
				// Only enhance voxels that actually have energy
				// Apply modest gamma correction and scaling for better visibility
				float gamma_corrected = std::pow(total_energy, 0.4f);
				float enhanced_energy = gamma_corrected * 2.0f;  // More conservative scaling
				enhanced_energy = std::clamp(enhanced_energy, 0.0f, 1.0f);
				
				// Use original color mapping approach with conservative energy range
				float min_energy = 0.0f;
				float max_energy = 1.0f;
				
				float min_alpha = use_emittance ? 0.1f : 0.2f;
				float max_alpha = use_emittance ? 0.7f : 0.9f;
				float alpha = min_alpha + (max_alpha - min_alpha) * enhanced_energy * (1.0f - 0.2f * normalized_dist);
				
				color = layer_energy_color(enhanced_energy, min_energy, max_energy, voxel->layer_id);
				color.a = alpha;
			}
			else if (voxel->material != nullptr) {
				// Keep voxels with no energy very dim, same as layer mode
				color = layer_energy_color(MathConstants::ENERGY_THRESHOLD, 0.0f, 1.0f, voxel->layer_id);
				color.a = 0.05f;
			}
		}

		// Only add voxels that should be visible
		if (color.a > 0.0f) {
			float depth = glm::dot(camera_offset, camera_front);
			voxel_instances_.push_back({voxel_pos, color, voxel_scale, depth});
		}
	}

	voxel_instances_dirty_ = false;
}

void VoxelRenderer::update_cached_energy_range(const Simulator& simulator, const Settings& settings) const {
	if (energy_range_cached_ && cached_range_mode_ == settings.voxel_mode) {
		return;
	}

	std::vector<float> all_energies;

	// Collect energies from photon paths
	for (const Photon& photon : simulator.photons) {
		if (photon.path_head) {
			auto current = photon.path_head;
			while (current) {
				float energy = static_cast<float>(current->value);
				if (energy > 0.0f) {
					all_energies.push_back(energy);
				}
				current = current->next;
			}
		}
	}

	// Collect energies from voxels
	const auto& config = Config::get();
	const uint32_t nx = static_cast<uint32_t>(config.nx());
	const uint32_t ny = static_cast<uint32_t>(config.ny());
	const uint32_t nz = static_cast<uint32_t>(config.nz());

	if (nx > 0 && ny > 0 && nz > 0) {
		const uint32_t total_voxels = nx * ny * nz;

		all_energies.reserve(all_energies.size() + total_voxels / 10);

		for (uint32_t linear_idx = 0; linear_idx < total_voxels; ++linear_idx) {
			const uint32_t iz = linear_idx / (nx * ny);
			const uint32_t iy = (linear_idx % (nx * ny)) / nx;
			const uint32_t ix = linear_idx % nx;

			Voxel* voxel = simulator.voxel_grid(ix, iy, iz);
			if (voxel) {
				float total_energy;

				if (settings.voxel_mode == VoxelMode::Layers) {
					total_energy = (voxel->material != nullptr) ? 1.0f : 0.0f;
				} else {
					// Use the same simplified calculation for range
					total_energy = calculate_voxel_energy_percentage(voxel, settings.voxel_mode, simulator);
				}

				if (total_energy > 0.0000001f) {
					all_energies.push_back(total_energy);
				}
			}
		}
	}

	// Calculate range using percentile-based method
	cached_min_energy_ = MathConstants::ENERGY_THRESHOLD;
	cached_max_energy_ = 0.01f;

	if (!all_energies.empty()) {
		std::ranges::sort(all_energies);
		const size_t count = all_energies.size();
		const float p5 = all_energies[static_cast<size_t>(count * 0.05)];
		const float p95 = all_energies[static_cast<size_t>(count * 0.95)];

		cached_min_energy_ = std::max(p5, MathConstants::ENERGY_THRESHOLD);
		cached_max_energy_ = std::max(p95, cached_min_energy_ * 10.0f);
	}

	cached_range_mode_ = settings.voxel_mode;
	energy_range_cached_ = true;
}

void VoxelRenderer::update_cached_global_energy(const Simulator& simulator, VoxelMode mode) const {
	if (global_energy_cached_ && cached_energy_mode_ == mode) {
		return;
	}

	// Calculate total energy across ALL mediums for each energy type
	cached_global_absorption_total_ = 0.0;
	cached_global_emittance_total_ = 0.0;
	cached_global_combined_total_ = 0.0;

	const auto& mediums = simulator.get_mediums();
	for (const auto& medium : mediums) {
		const auto& volume = medium.get_volume();
		for (const auto& voxel_ptr : volume) {
			if (voxel_ptr && voxel_ptr->material) {
				cached_global_absorption_total_ += voxel_ptr->absorption;
				cached_global_emittance_total_ += voxel_ptr->total_emittance();
				cached_global_combined_total_ += voxel_ptr->absorption + voxel_ptr->total_emittance();
			}
		}
	}

	cached_energy_mode_ = mode;
	global_energy_cached_ = true;
}

void VoxelRenderer::end_voxel_instances(VoxelMode mode, const Camera& camera) {
	if (voxel_instances_.empty())
		return;

	// Only skip processing if truly nothing changed
	if (!voxel_instances_dirty_ && voxel_buffer_uploaded_) {
		// Check if we have a completed background sort ready
		if (background_sort_ready_.load() && !background_sort_in_progress_.load()) {
			// Background sort is complete - upload to GPU
			glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);
			glBufferData(GL_ARRAY_BUFFER,
						 background_sorted_voxels_.size() * sizeof(VoxelInstance),
						 background_sorted_voxels_.data(),
						 GL_STATIC_DRAW);

			// Swap the sorted data back to main instances
			voxel_instances_ = std::move(background_sorted_voxels_);
			background_sort_ready_.store(false);

			// Update cached camera position
			cached_camera_position_ = camera.get_position();
			cached_camera_target_ = camera.get_target();
		}

		// Check if we need to start a new background sort
		if (!background_sort_in_progress_.load()) {
			glm::vec3 current_camera_pos = camera.get_position();
			glm::vec3 current_camera_target = camera.get_target();

			// Only start background sort if camera moved significantly
			const float SIGNIFICANT_MOVEMENT = 0.05f;
			bool significant_camera_movement =
				glm::length(current_camera_pos - cached_camera_position_) > SIGNIFICANT_MOVEMENT
				|| glm::length(current_camera_target - cached_camera_target_) > SIGNIFICANT_MOVEMENT;

			// Limit background sorting to reasonable sizes
			if (significant_camera_movement && voxel_instances_.size() < 100000) {
				// Start background sorting
				background_sort_in_progress_.store(true);
				background_sorted_voxels_ = voxel_instances_; // Copy current data

				// Launch async sorting task
				sorting_future_ = std::async(std::launch::async, [this, current_camera_pos, current_camera_target]() {
					// Calculate depths on background thread
					glm::vec3 camera_front = glm::normalize(current_camera_target - current_camera_pos);

					for (auto& instance : background_sorted_voxels_) {
						glm::vec3 view_dir = instance.position - current_camera_pos;
						instance.depth = glm::dot(view_dir, camera_front);
					}

					// Sort on background thread
					std::sort(background_sorted_voxels_.begin(),
							  background_sorted_voxels_.end(),
							  [](const VoxelInstance& a, const VoxelInstance& b) {
								  return a.depth > b.depth; // Back to front
							  });

					// Mark as ready for GPU upload
					background_sort_ready_.store(true);
					background_sort_in_progress_.store(false);
				});
			}
		}

		// Current frame continues with existing data
		return;
	}

	// Store mode for consistency
	current_voxel_mode_ = mode;

	// Balance performance vs quality - sort small datasets
	if (voxel_instances_.size() < 50000) {
		// Sort by depth (back to front) for optimal transparency
		std::sort(voxel_instances_.begin(), voxel_instances_.end(), [](const VoxelInstance& a, const VoxelInstance& b) {
			return a.depth > b.depth; // Back to front
		});
	}

	// Upload instance data to GPU
	glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);
	glBufferData(
		GL_ARRAY_BUFFER, voxel_instances_.size() * sizeof(VoxelInstance), voxel_instances_.data(), GL_STATIC_DRAW);

	voxel_buffer_uploaded_ = true;
	voxel_instances_dirty_ = false;

	// Cache current camera position for change detection
	cached_camera_position_ = camera.get_position();
	cached_camera_target_ = camera.get_target();
}

void VoxelRenderer::draw_voxel_instances(const Camera& camera) {
	if (voxel_instances_.empty() || !voxel_shader_) {
		return;
	}

	// Enable transparency - KEEP ORIGINAL APPEARANCE (match original Renderer behavior)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Disable depth writing for transparent objects but keep depth testing
	glDepthMask(GL_FALSE);
	glEnable(GL_DEPTH_TEST);

	glUseProgram(voxel_shader_);

	// Set MVP matrix
	glm::mat4 mvp = camera.get_mvp_matrix();
	glUniformMatrix4fv(voxel_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

	// Bind VAO and draw instanced
	glBindVertexArray(voxel_vao_);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(voxel_instances_.size()));

	// Restore OpenGL state (match original Renderer behavior)
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glBindVertexArray(0);
}

bool VoxelRenderer::setup_voxel_rendering() {
	// Load voxel shaders
	std::string vertex_source = Shader::load_shader_source("shaders/voxels.vert");
	std::string fragment_source = Shader::load_shader_source("shaders/voxels.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load voxel shader sources" << std::endl;
		return false;
	}

	// Create shader program
	GLuint vertex_shader = Shader::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = Shader::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader)
			glDeleteShader(vertex_shader);
		if (fragment_shader)
			glDeleteShader(fragment_shader);
		return false;
	}

	voxel_shader_ = glCreateProgram();
	glAttachShader(voxel_shader_, vertex_shader);
	glAttachShader(voxel_shader_, fragment_shader);
	glLinkProgram(voxel_shader_);

	GLint success;
	glGetProgramiv(voxel_shader_, GL_LINK_STATUS, &success);
	if (!success) {
		char info_log[512];
		glGetProgramInfoLog(voxel_shader_, 512, nullptr, info_log);
		std::cerr << "Voxel shader program linking failed: " << info_log << std::endl;
		glDeleteProgram(voxel_shader_);
		voxel_shader_ = 0;
		glDeleteShader(vertex_shader);
		glDeleteShader(fragment_shader);
		return false;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	// Cache uniform location
	voxel_mvp_uniform_location_ = glGetUniformLocation(voxel_shader_, "uMVP");

	// Create VAO and VBOs
	glGenVertexArrays(1, &voxel_vao_);
	glGenBuffers(1, &voxel_vbo_);
	glGenBuffers(1, &voxel_instance_vbo_);

	glBindVertexArray(voxel_vao_);

	// Setup base voxel geometry
	if (!setup_voxel_geometry()) {
		std::cerr << "Failed to setup voxel geometry" << std::endl;
		return false;
	}

	// Setup instance data buffer
	glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);

	// Instance position (location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, position));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// Instance color (location 3)
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// Instance scale (location 4)
	glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, scale));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	glBindVertexArray(0);
	return true;
}

bool VoxelRenderer::setup_voxel_geometry() {
	// Create unit cube geometry (centered at origin)
	float vertices[] = {// Front face
						-0.5f,
						-0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,
						0.5f,
						-0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,
						0.5f,
						0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,
						0.5f,
						0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,
						-0.5f,
						0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,
						-0.5f,
						-0.5f,
						0.5f,
						0.0f,
						0.0f,
						1.0f,

						// Back face
						-0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,
						0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,
						0.5f,
						0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,
						0.5f,
						0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,
						-0.5f,
						0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,
						-0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						0.0f,
						-1.0f,

						// Left face
						-0.5f,
						0.5f,
						0.5f,
						-1.0f,
						0.0f,
						0.0f,
						-0.5f,
						0.5f,
						-0.5f,
						-1.0f,
						0.0f,
						0.0f,
						-0.5f,
						-0.5f,
						-0.5f,
						-1.0f,
						0.0f,
						0.0f,
						-0.5f,
						-0.5f,
						-0.5f,
						-1.0f,
						0.0f,
						0.0f,
						-0.5f,
						-0.5f,
						0.5f,
						-1.0f,
						0.0f,
						0.0f,
						-0.5f,
						0.5f,
						0.5f,
						-1.0f,
						0.0f,
						0.0f,

						// Right face
						0.5f,
						0.5f,
						0.5f,
						1.0f,
						0.0f,
						0.0f,
						0.5f,
						0.5f,
						-0.5f,
						1.0f,
						0.0f,
						0.0f,
						0.5f,
						-0.5f,
						-0.5f,
						1.0f,
						0.0f,
						0.0f,
						0.5f,
						-0.5f,
						-0.5f,
						1.0f,
						0.0f,
						0.0f,
						0.5f,
						-0.5f,
						0.5f,
						1.0f,
						0.0f,
						0.0f,
						0.5f,
						0.5f,
						0.5f,
						1.0f,
						0.0f,
						0.0f,

						// Bottom face
						-0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						-1.0f,
						0.0f,
						0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						-1.0f,
						0.0f,
						0.5f,
						-0.5f,
						0.5f,
						0.0f,
						-1.0f,
						0.0f,
						0.5f,
						-0.5f,
						0.5f,
						0.0f,
						-1.0f,
						0.0f,
						-0.5f,
						-0.5f,
						0.5f,
						0.0f,
						-1.0f,
						0.0f,
						-0.5f,
						-0.5f,
						-0.5f,
						0.0f,
						-1.0f,
						0.0f,

						// Top face
						-0.5f,
						0.5f,
						-0.5f,
						0.0f,
						1.0f,
						0.0f,
						0.5f,
						0.5f,
						-0.5f,
						0.0f,
						1.0f,
						0.0f,
						0.5f,
						0.5f,
						0.5f,
						0.0f,
						1.0f,
						0.0f,
						0.5f,
						0.5f,
						0.5f,
						0.0f,
						1.0f,
						0.0f,
						-0.5f,
						0.5f,
						0.5f,
						0.0f,
						1.0f,
						0.0f,
						-0.5f,
						0.5f,
						-0.5f,
						0.0f,
						1.0f,
						0.0f};

	glBindBuffer(GL_ARRAY_BUFFER, voxel_vbo_);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Normal attribute (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	return true;
}

glm::vec4 VoxelRenderer::layer_energy_color(float energy, float min_energy, float max_energy, uint8_t layer_id) const {
	// Clamp and normalize energy like the original function
	energy = std::clamp(energy, min_energy, max_energy);

	float linear_normalized = (energy - min_energy) / (max_energy - min_energy);
	linear_normalized = std::clamp(linear_normalized, 0.0f, 1.0f);

	float normalized = std::pow(linear_normalized, 0.3f);
	normalized = std::clamp(normalized, 0.0f, 1.0f);

	// Define base colors for different material types
	glm::vec3 base_colors[6] = {
		glm::vec3(1.0f, 0.4f, 0.4f), // Red theme for layer 0
		glm::vec3(0.4f, 0.4f, 1.0f), // Blue theme for layer 1
		glm::vec3(0.4f, 1.0f, 0.4f), // Green theme for layer 2
		glm::vec3(1.0f, 0.4f, 1.0f), // Magenta theme for layer 3
		glm::vec3(1.0f, 1.0f, 0.4f), // Yellow theme for layer 4
		glm::vec3(0.4f, 1.0f, 1.0f), // Cyan theme for layer 5
	};

	glm::vec3 base_color = base_colors[layer_id % 6];

	// Create energy gradient within the material's color theme
	if (normalized > 0.85f) {
		// Very high energy: brighten to near white
		return glm::vec4(glm::mix(base_color, glm::vec3(1.0f), 0.6f), 1.0f);
	}
	else if (normalized > 0.70f) {
		// High energy: bright version of base color
		float t = (normalized - 0.70f) / 0.15f;
		glm::vec3 bright_color = base_color * 1.4f;
		bright_color = glm::min(bright_color, glm::vec3(1.0f));
		return glm::vec4(glm::mix(base_color, bright_color, t), 1.0f);
	}
	else if (normalized > 0.40f) {
		// Medium energy: full base color
		return glm::vec4(base_color, 1.0f);
	}
	else if (normalized > 0.15f) {
		// Low-medium energy: darken the base color
		float t = (normalized - 0.15f) / 0.25f;
		glm::vec3 dark_color = base_color * 0.7f;
		return glm::vec4(glm::mix(dark_color, base_color, t), 1.0f);
	}
	else if (normalized > 0.05f) {
		// Low energy: darker version
		float t = (normalized - 0.05f) / 0.10f;
		glm::vec3 darker_color = base_color * 0.5f;
		glm::vec3 dark_color = base_color * 0.7f;
		return glm::vec4(glm::mix(darker_color, dark_color, t), 1.0f);
	}
	else {
		// Very low energy: very dark but still visible
		return glm::vec4(base_color * 0.3f, 1.0f);
	}
}

float VoxelRenderer::calculate_voxel_energy_percentage(const Voxel* voxel, VoxelMode mode, const Simulator& simulator) const {
	if (!voxel || !voxel->material) {
		return 0.0f;
	}

	// Ensure global energy totals are cached
	update_cached_global_energy(simulator, mode);

	// Get the voxel's absolute energy value
	double voxel_energy = 0.0;
	if (mode == VoxelMode::Absorption) {
		voxel_energy = voxel->absorption;
	} else if (mode == VoxelMode::Emittance) {
		voxel_energy = voxel->total_emittance();
	} else {
		voxel_energy = voxel->absorption + voxel->total_emittance();
	}

	// If there's no energy in this voxel, return 0
	if (voxel_energy <= 0.0) {
		return 0.0f;
	}

	// Use cached global total energy for this mode
	double total_energy = 0.0;
	if (mode == VoxelMode::Absorption) {
		total_energy = cached_global_absorption_total_;
	} else if (mode == VoxelMode::Emittance) {
		total_energy = cached_global_emittance_total_;
	} else {
		total_energy = cached_global_combined_total_;
	}

	// Return percentage of this voxel relative to global total (0.0 to 1.0)
	// This will be used by the percentile-based discretization system
	if (total_energy > 0.0) {
		return static_cast<float>(voxel_energy / total_energy);
	}
	return 0.0f;
}
