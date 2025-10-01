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
#include "simulator/simulator.hpp"
#include "simulator/config.hpp"
#include "simulator/layer.hpp"
#include "simulator/medium.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"
#include "simulator/voxel.hpp"

Renderer::Renderer() = default;

Renderer::~Renderer() {
	// Ensure background sorting completes before destroying resources
	if (background_sort_in_progress_.load() && sorting_future_.valid()) {
		sorting_future_.wait(); // Wait for background sort to complete
	}

	// Clean up advanced instanced rendering systems
	if (voxel_vao_) {
		glDeleteVertexArrays(1, &voxel_vao_);
	}
	if (voxel_vbo_) {
		glDeleteBuffers(1, &voxel_vbo_);
	}
	if (voxel_instance_vbo_) {
		glDeleteBuffers(1, &voxel_instance_vbo_);
	}
	if (lines_vao_) {
		glDeleteVertexArrays(1, &lines_vao_);
	}
	if (lines_vbo_) {
		glDeleteBuffers(1, &lines_vbo_);
	}
	if (lines_instance_vbo_) {
		glDeleteBuffers(1, &lines_instance_vbo_);
	}
	if (points_vao_) {
		glDeleteVertexArrays(1, &points_vao_);
	}
	if (points_vbo_) {
		glDeleteBuffers(1, &points_vbo_);
	}
	if (points_instance_vbo_) {
		glDeleteBuffers(1, &points_instance_vbo_);
	}
	if (triangles_vao_) {
		glDeleteVertexArrays(1, &triangles_vao_);
	}
	if (triangles_vbo_) {
		glDeleteBuffers(1, &triangles_vbo_);
	}
	if (triangles_instance_vbo_) {
		glDeleteBuffers(1, &triangles_instance_vbo_);
	}
	if (wireframe_instance_vbo_) {
		glDeleteBuffers(1, &wireframe_instance_vbo_);
	}

	// Clean up medium line VAO and VBO
	if (wireframe_vao_) {
		glDeleteVertexArrays(1, &wireframe_vao_);
	}
	if (wireframe_instance_vbo_) {
		glDeleteBuffers(1, &wireframe_instance_vbo_);
	}

	// Clean up all compiled shader programs
	if (voxel_shader_) {
		glDeleteProgram(voxel_shader_);
	}
	if (lines_shader_) {
		glDeleteProgram(lines_shader_);
	}
	if (points_shader_) {
		glDeleteProgram(points_shader_);
	}
	if (triangles_shader_) {
		glDeleteProgram(triangles_shader_);
	}
}

bool Renderer::initialize() {
	// Log initialization if debugging is enabled
	if (Config::is_initialized() && Config::get().log()) {
		std::cout << "Initializing modern OpenGL 4.5 renderer..." << std::endl;
	}

	// Configure OpenGL state
	setup_opengl();

	// Instanced voxel rendering for volume visualization
	if (!setup_voxel_instanced_rendering()) {
		std::cerr << "Failed to setup instanced voxel rendering" << std::endl;
		return false;
	}

	// Instanced line rendering for photon paths
	if (!setup_line_instanced_rendering()) {
		std::cerr << "Failed to setup instanced line rendering" << std::endl;
		return false;
	}

	if (!setup_point_instanced_rendering()) {
		std::cerr << "Failed to setup instanced point rendering" << std::endl;
		return false;
	}

	// Instanced medium geometry rendering
	if (!setup_triangle_instanced_rendering()) {
		std::cerr << "Failed to setup instanced triangle rendering" << std::endl;
		return false;
	}

	// Setup medium line VAO (uses same shader program but separate VAO/VBO)
	if (!setup_medium_line_vao()) {
		std::cerr << "Failed to setup medium line VAO" << std::endl;
		return false;
	}

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

	// Set clear color to dark charcoal gray - easy on the eyes and neutral
	// Avoids conflicts with layer colors (red, blue, green, magenta, yellow, cyan)
	glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
}

void Renderer::update() {
	update_camera();

	// Apply smooth movement for FPS mode
	if (!orbit_camera_mode_) {
		const float move_speed = 0.02f; // Reduced speed for smoother movement
		bool camera_moved = false;

		if (key_state_.w_pressed) {
			camera_.move_forward(move_speed);
			camera_moved = true;
		}
		if (key_state_.s_pressed) {
			camera_.move_backward(move_speed);
			camera_moved = true;
		}
		if (key_state_.a_pressed) {
			camera_.move_left(move_speed);
			camera_moved = true;
		}
		if (key_state_.d_pressed) {
			camera_.move_right(move_speed);
			camera_moved = true;
		}
		if (key_state_.space_pressed) {
			camera_.move_up(move_speed);
			camera_moved = true;
		}
		if (key_state_.shift_pressed) {
			camera_.move_down(move_speed);
			camera_moved = true;
		}

		// Mark camera state as changed if any movement occurred
		if (camera_moved) {
			camera_state_changed_ = true;
		}
	}
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

	// Update GUI element screen positions
	update_energy_label_screen_positions();

	// Render geometry boundaries and wireframes
	if (settings_.draw_volume) {
		draw_volume(simulator); // Medium geometry visualization
	}

	// Clear depth buffer for proper layering of volume elements
	glClear(GL_DEPTH_BUFFER_BIT);

	// Render voxelized absorption data using instanced rendering
	draw_voxels(settings_);

	// Render photon transport paths with proper depth testing
	if (settings_.draw_paths) {
		draw_paths(settings_); // Instanced path visualization
	}

	// Render energy labels as screen-space billboards
	draw_labels(settings_);
}

void Renderer::draw_volume(const Simulator& simulator) {
	// Combine wireframe (lines) and faces (triangles) for a complete boundary visualization
	// Use instanced rendering for better performance with complex meshes
	triangle_instances_.clear();
	wireframe_instances_.clear();

	// Face culling is already disabled by default, but ensure blending is set up properly
	enable_blending();
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glm::vec4 wireframe_color(0.7f, 0.7f, 0.7f, 0.8f); // Bright wireframe
	glm::vec4 face_color(0.2f, 0.2f, 0.2f, 0.1f);      // Subtle transparent faces

	// Process each layer's geometry
	const auto& layers = simulator.get_all_layers();
	for (const auto& layer : layers) {
		// For each triangle in the layer's mesh, we'll determine if it's part of a planar face
		// and group triangles that share the same plane into faces

		std::map<std::string, std::vector<glm::vec3>> planar_faces;

		for (const auto& triangle : layer.mesh) {
			glm::vec3 v0 = to_float(triangle.v0());
			glm::vec3 v1 = to_float(triangle.v1());
			glm::vec3 v2 = to_float(triangle.v2());

			// Calculate the triangle normal
			glm::vec3 normal = normalize(cross(v1 - v0, v2 - v0));

			// Calculate plane equation: ax + by + cz = d
			float d = dot(normal, v0);

			// Create a key for this plane (quantized to handle floating point precision)
			// clang-format off
			std::string plane_key = std::format("{},{},{},{}", 
												int(normal.x * 1000), 
												int(normal.y * 1000), 
												int(normal.z * 1000), 
												int(d * 1000));
			// clang-format on

			// Add vertices to this plane's face
			auto& face_vertices = planar_faces[plane_key];
			face_vertices.push_back(v0);
			face_vertices.push_back(v1);
			face_vertices.push_back(v2);
		}

		// Now process each planar face
		for (const auto& face_pair : planar_faces) {
			const std::vector<glm::vec3>& vertices = face_pair.second;

			if (vertices.size() == 3) {
				// Single triangle - render as triangle
				glm::vec3 normal = normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
				add_triangle_instance(vertices[0], vertices[1], vertices[2], face_color, normal);

				// Add wireframe edges
				add_medium_line_instance(vertices[0], vertices[1], wireframe_color);
				add_medium_line_instance(vertices[1], vertices[2], wireframe_color);
				add_medium_line_instance(vertices[2], vertices[0], wireframe_color);
			}
			else if (vertices.size() >= 6 && vertices.size() % 3 == 0) {
				// Multiple triangles forming a face - check if they form a quad

				// For now, let's find unique vertices and try to form a quadrilateral
				std::vector<glm::vec3> unique_vertices;
				const float epsilon = MathConstants::UV_MATCH_EPSILON;

				// Collect unique vertices
				for (const auto& v : vertices) {
					bool is_unique = true;
					for (const auto& uv : unique_vertices) {
						if (length(v - uv) < epsilon) {
							is_unique = false;
							break;
						}
					}
					if (is_unique) {
						unique_vertices.push_back(v);
					}
				}

				if (unique_vertices.size() == 4) {
					// We have a quadrilateral! Order the vertices properly
					// Find the centroid
					glm::vec3 center(0.0f);
					for (const auto& v : unique_vertices) {
						center += v;
					}
					center /= static_cast<float>(unique_vertices.size());

					// Calculate the face normal from first triangle
					glm::vec3 normal = normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));

					// Create a reference vector perpendicular to normal
					glm::vec3 ref_vec =
						abs(normal.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);
					glm::vec3 tangent = normalize(cross(normal, ref_vec));
					glm::vec3 bitangent = normalize(cross(normal, tangent));

					// Sort vertices by angle around the center
					// clang-format off
					std::ranges::sort(unique_vertices,
						[&center, &tangent, &bitangent](const glm::vec3& a, const glm::vec3& b) noexcept {
							const glm::vec3 dir_a = a - center;
							const glm::vec3 dir_b = b - center;

							const float angle_a = atan2f(dot(dir_a, bitangent), dot(dir_a, tangent));
							const float angle_b = atan2f(dot(dir_b, bitangent), dot(dir_b, tangent));

							return angle_a < angle_b;
					});
					// clang-format on

					// Render as two triangles forming a quad
					add_triangle_instance(
						unique_vertices[0], unique_vertices[1], unique_vertices[2], face_color, normal);
					add_triangle_instance(
						unique_vertices[0], unique_vertices[2], unique_vertices[3], face_color, normal);

					// Add wireframe edges for the quad
					add_medium_line_instance(unique_vertices[0], unique_vertices[1], wireframe_color);
					add_medium_line_instance(unique_vertices[1], unique_vertices[2], wireframe_color);
					add_medium_line_instance(unique_vertices[2], unique_vertices[3], wireframe_color);
					add_medium_line_instance(unique_vertices[3], unique_vertices[0], wireframe_color);
				}
				else {
					// Fallback: render all triangles individually
					for (size_t i = 0; i < vertices.size(); i += 3) {
						glm::vec3 normal =
							normalize(cross(vertices[i + 1] - vertices[i], vertices[i + 2] - vertices[i]));
						add_triangle_instance(vertices[i], vertices[i + 1], vertices[i + 2], face_color, normal);

						// Add wireframe edges
						add_medium_line_instance(vertices[i], vertices[i + 1], wireframe_color);
						add_medium_line_instance(vertices[i + 1], vertices[i + 2], wireframe_color);
						add_medium_line_instance(vertices[i + 2], vertices[i], wireframe_color);
					}
				}
			}
		}
	}

	// End collecting instances and draw
	triangle_instances_uploaded_ = false;
	medium_line_instances_uploaded_ = false;

	// Render triangles first, then lines to maintain exact rendering order
	draw_triangle_instances();
	draw_medium_line_instances();
}

void Renderer::update_camera() {
	// Update camera matrices for viewport
	float aspect = static_cast<float>(viewport_width_) / static_cast<float>(viewport_height_);
	camera_.set_aspect_ratio(aspect);
}

void Renderer::update_camera_target(const Simulator& simulator) {
	Range3 bounds = simulator.get_combined_bounds();

	// Set camera target to center of bounds
	glm::vec3 center = to_float(bounds.center());
	camera_.set_target(center);

	// Set material elevation bounds for scroll wheel constraint (Y-axis)
	camera_.set_elevation_bounds(static_cast<float>(bounds.min_bounds.y), static_cast<float>(bounds.max_bounds.y));
}

void Renderer::set_viewport(int width, int height) {
	viewport_width_ = width;
	viewport_height_ = height;
	glViewport(0, 0, width, height);

	// Invalidate energy label cache when viewport changes to ensure proper coordinate conversion
	invalidate_energy_label_cache();
}

// Input handlers (simplified for now)
void Renderer::handle_key_input(int key, int /*scancode*/, int action, int /*mods*/) {
	// Check if WASD is pressed while in Orbit mode - switch to Free mode automatically
	if (orbit_camera_mode_ && action == GLFW_PRESS) {
		if (key == GLFW_KEY_W || key == GLFW_KEY_A || key == GLFW_KEY_S || key == GLFW_KEY_D) {
			// Switch to Free camera mode
			set_camera_mode(false); // false = Free mode

			// Notify the UI about the mode change
			if (camera_mode_change_callback_) {
				camera_mode_change_callback_(false); // false = Free mode
			}
		}
	}

	// Track key states for smooth movement
	if (key == GLFW_KEY_W) {
		key_state_.w_pressed = (action != GLFW_RELEASE);
	}
	if (key == GLFW_KEY_A) {
		key_state_.a_pressed = (action != GLFW_RELEASE);
	}
	if (key == GLFW_KEY_S) {
		key_state_.s_pressed = (action != GLFW_RELEASE);
	}
	if (key == GLFW_KEY_D) {
		key_state_.d_pressed = (action != GLFW_RELEASE);
	}
	if (key == GLFW_KEY_SPACE) {
		key_state_.space_pressed = (action != GLFW_RELEASE);
	}
	if (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT) {
		key_state_.shift_pressed = (action != GLFW_RELEASE);
	}
}

void Renderer::handle_mouse_move(float xpos, float ypos) {
	camera_.handle_mouse_move(xpos, ypos);

	// Mark camera state as changed to update label positions
	camera_state_changed_ = true;
}

void Renderer::handle_mouse_button(int button, int action, int /*mods*/) {
	camera_.handle_mouse_button(button, action);
}

void Renderer::handle_mouse_scroll(float /*xoffset*/, float yoffset) {
	camera_.handle_mouse_scroll(yoffset);

	// Mark camera state as changed to update label positions
	camera_state_changed_ = true;
}

void Renderer::reset_camera() {
	// Reset camera to initial position and angles
	camera_.reset_to_initial();
}

void Renderer::set_camera_mode(bool is_arc_mode) {
	orbit_camera_mode_ = is_arc_mode;
	camera_.set_fps_mode(!is_arc_mode); // FPS mode when not Orbit mode
}

bool Renderer::should_capture_mouse() const {
	return camera_.should_capture_mouse();
}

// Instanced voxel rendering
void Renderer::draw_voxels(const Settings& settings) {
	if (!simulator_) {
		return;
	}

	// Only draw voxels if voxel rendering is enabled
	if (!settings.draw_voxels) {
		return;
	}

	// Check if we need to rebuild voxel instances
	// Only rebuild when simulation data changes or voxel mode changes
	static VoxelMode last_voxel_mode = VoxelMode::Layers;

	bool need_rebuild = voxel_instances_dirty_ || settings.voxel_mode != last_voxel_mode;

	if (!need_rebuild && voxel_buffer_uploaded_ && !voxel_instances_.empty()) {
		// Use cached voxel instances, but still need to update depth sorting for camera changes
		end_voxel_instances(settings.voxel_mode);
		draw_voxel_instances();
		return;
	}

	// Mark instances as dirty if we're rebuilding due to mode change
	if (settings.voxel_mode != last_voxel_mode) {
		voxel_instances_dirty_ = true;
		voxel_buffer_uploaded_ = false;
	}

	last_voxel_mode = settings.voxel_mode;

	// Setup OpenGL state for transparent voxel rendering (match original)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE); // Disable depth writing for transparency
	glDisable(GL_CULL_FACE);

	// Get MCML grid parameters
	const auto& config = Config::get();
	const uint32_t nx = static_cast<uint32_t>(config.nx());
	const uint32_t ny = static_cast<uint32_t>(config.ny());
	const uint32_t nz = static_cast<uint32_t>(config.nz());
	double vox_size = config.vox_size();
	double half_vox_size = vox_size * 0.5;
	float voxel_scale = static_cast<float>(vox_size * 0.95); // Slightly smaller to show gaps

	// Use cached energy range instead of expensive per-frame analysis
	update_cached_energy_range(settings);

	// Set up instanced rendering
	if (voxel_instances_dirty_) {
		voxel_instances_.clear();
	}

	// Pre-calculate camera data and bounds once for all voxels
	glm::vec3 camera_pos = camera_.get_position();
	glm::vec3 camera_target = camera_.get_target();
	glm::vec3 camera_front = glm::normalize(camera_target - camera_pos);
	float max_render_distance = 50.0f; // Don't render voxels too far away

	// Cache bounds calculation to avoid repeated calls
	auto bounds = simulator_->get().get_combined_bounds();

	// Render voxels using cached energy scaling with optimized iteration
	const uint32_t total_voxels = nx * ny * nz;

	// Pre-calculate constants outside the loop
	const float vox_size_f = static_cast<float>(vox_size);
	const float half_vox_size_f = static_cast<float>(half_vox_size);
	const glm::vec3 bounds_min = to_float(bounds.min_bounds);

	// Pre-determine energy calculation method
	const bool use_absorption = (settings.voxel_mode == VoxelMode::Absorption);
	const bool use_emittance = (settings.voxel_mode == VoxelMode::Emittance);
	const bool use_layers = (settings.voxel_mode == VoxelMode::Layers);

	// Cache distance calculation constants
	static float cached_max_dist = 0.0f;
	static bool max_dist_cached = false;
	if (!max_dist_cached) {
		cached_max_dist = glm::length(to_float(bounds.max_bounds));
		max_dist_cached = true;
	}

	// Pre-calculate alpha constants for each mode
	float min_alpha = 0.2f;
	float max_alpha = 0.9f;
	if (use_emittance) {
		min_alpha = 0.1f;
		max_alpha = 0.7f;
	}

	// Single linear loop with better cache locality (z-major order)
	for (uint32_t linear_idx = 0; linear_idx < total_voxels; ++linear_idx) {
		// Convert linear index to 3D coordinates (z-major for better cache performance)
		const uint32_t iz = linear_idx / (nx * ny);
		const uint32_t iy = (linear_idx % (nx * ny)) / nx;
		const uint32_t ix = linear_idx % nx;

		Voxel* voxel = simulator_->get().voxel_grid(ix, iy, iz);
		if (!voxel || !voxel->material) {
			continue; // Skip voxels with no material
		}

		// Calculate world position using pre-calculated constants
		// clang-format off
		const glm::vec3 voxel_pos = bounds_min + glm::vec3(vox_size_f * ix + half_vox_size_f,
														   vox_size_f * iy + half_vox_size_f,
														   vox_size_f * iz + half_vox_size_f);
		// clang-format on

		// Early distance culling with squared distance to avoid sqrt
		const glm::vec3 camera_offset = voxel_pos - camera_pos;
		const float distance_squared = glm::dot(camera_offset, camera_offset);
		const float max_distance_squared = max_render_distance * max_render_distance;

		if (distance_squared > max_distance_squared) {
			continue; // Skip voxels too far from camera
		}

		// Calculate energy based on pre-determined mode
		float total_energy;
		if (use_absorption) {
			total_energy = static_cast<float>(voxel->absorption);
		}
		else if (use_emittance) {
			total_energy = static_cast<float>(voxel->total_emittance());
		}
		else if (use_layers) {
			total_energy = 1.0f; // We already know voxel->material != nullptr
		}
		else {
			// Combined mode: calculate both values
			total_energy = static_cast<float>(voxel->absorption + voxel->total_emittance());
		}

		glm::vec4 color(0.0f);

		if (use_layers) {
			// Layers mode: show all material voxels regardless of energy
			color = layer_energy_color(
				MathConstants::ENERGY_THRESHOLD, cached_min_energy_, cached_max_energy_, voxel->layer_id);
			color.a = 0.05f; // Same very faint alpha as no-absorption voxels
		}
		else {
			// Original energy-based rendering logic for other modes

			// Use gamma correction for better mid-range contrast
			float normalized_energy = (total_energy - cached_min_energy_) / (cached_max_energy_ - cached_min_energy_);
			normalized_energy = std::clamp(normalized_energy, 0.0f, 1.0f);

			// Apply gamma correction for better contrast (gamma = 0.5 brightens mid-tones)
			float gamma_corrected = std::pow(normalized_energy, 0.5f);

			// Calculate distance from origin for depth-based alpha (using pre-calculated max_dist)
			float dist_from_origin = glm::length(voxel_pos);
			float normalized_dist = glm::clamp(dist_from_origin / cached_max_dist, 0.0f, 1.0f);

			// More visible alpha scaling based on gamma-corrected energy
			float alpha = min_alpha + (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.2f * normalized_dist);

			// Use layer-specific energy color mapping with dynamic range
			if (total_energy > MathConstants::ENERGY_THRESHOLD) {
				color = layer_energy_color(total_energy, cached_min_energy_, cached_max_energy_, voxel->layer_id);
				color.a = alpha;
			}
			else if (voxel->material != nullptr) {
				// Voxel in medium but no recorded energy: very faint material-colored hint
				color = layer_energy_color(
					MathConstants::ENERGY_THRESHOLD, cached_min_energy_, cached_max_energy_, voxel->layer_id);
				color.a = 0.05f; // Slightly more visible
			}
		}

		// Only add voxels that should be visible
		if (color.a > 0.0f) {
			// Calculate depth for sorting using pre-calculated camera offset
			float depth = glm::dot(camera_offset, camera_front);
			voxel_instances_.push_back({voxel_pos, color, voxel_scale, depth});
		}
	}

	// Render all voxels in a single instanced draw call
	end_voxel_instances(settings.voxel_mode);
	draw_voxel_instances();

	// Restore OpenGL state
	glDepthMask(GL_TRUE);    // Re-enable depth writing
	glEnable(GL_DEPTH_TEST); // Keep depth testing enabled
	glDisable(GL_BLEND);     // Disable blending
}

void Renderer::draw_paths(const Settings& settings) {
	if (!simulator_) {
		return;
	}

	// Use cached energy range for better performance
	update_cached_energy_range(settings);

	// Create adaptive logarithmic mapping function using cached values
	auto adaptive_log_color = [this](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, cached_min_energy_, cached_max_energy_);
	};

	// Use incremental caching instead of rebuilding every frame
	size_t current_photon_count = simulator_->get().photons.size();

	// Check if we need to rebuild cache completely or just add new photons
	if (!path_instances_cached_ || current_photon_count < cached_photon_count_) {
		// Complete rebuild needed (first time or photons were removed)
		cached_line_instances_.clear();
		cached_point_instances_.clear();
		cached_photon_count_ = 0;

		// Reset buffer upload flags when cache is invalidated
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;

		// Also reset dynamic geometry buffer flags
		point_geometry_buffer_uploaded_ = false;
	}
	else if (current_photon_count > cached_photon_count_) {
		// Incremental update - only process new photons
		// Keep existing cache, just mark buffers for re-upload
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;

		// Also reset dynamic geometry buffer flags
		point_geometry_buffer_uploaded_ = false;
	}

	// Only process photons if we have new ones to add
	if (current_photon_count > cached_photon_count_) {
		// Use cached surface calculation
		const auto& layers = simulator_->get().get_all_layers();
		if (!surface_cached_ && !layers.empty()) {
			cached_surface_y_ = -1000.0f;
			for (const auto& layer : layers) {
				for (const auto& triangle : layer.mesh) {
					cached_surface_y_ = std::max(cached_surface_y_, static_cast<float>(triangle.v0().y));
					cached_surface_y_ = std::max(cached_surface_y_, static_cast<float>(triangle.v1().y));
					cached_surface_y_ = std::max(cached_surface_y_, static_cast<float>(triangle.v2().y));
				}
			}
			surface_cached_ = true;
		}

		// Process only new photons (incremental caching)
		const auto paths = simulator_->get().photons;
		for (size_t i = cached_photon_count_; i < paths.size(); ++i) {
			const Photon& photon = paths[i];
			if (photon.path_head) {
				// Generate connected line segments with energy gradient
				auto current = photon.path_head;
				auto next = current ? current->next : nullptr;

				// Generate incident ray from source to material surface (cached surface)
				if (current && !simulator_->get().sources.empty()) {
					glm::vec3 first_interaction(static_cast<float>(current->position.x),
												static_cast<float>(current->position.y),
												static_cast<float>(current->position.z));

					const Source& source = simulator_->get().sources[0];
					glm::vec3 source_pos(static_cast<float>(source.origin.x),
										 static_cast<float>(source.origin.y),
										 static_cast<float>(source.origin.z));
					glm::vec3 source_dir(static_cast<float>(source.direction.x),
										 static_cast<float>(source.direction.y),
										 static_cast<float>(source.direction.z));

					// Use cached surface calculation
					if (source_dir.y != 0.0f) {
						float t = (cached_surface_y_ - source_pos.y) / source_dir.y;
						glm::vec3 surface_entry = source_pos + t * source_dir;

						// Add incident ray
						glm::vec4 incident_color(1.0f, 1.0f, 1.0f, 1.0f);
						cached_line_instances_.push_back({source_pos, surface_entry, incident_color, incident_color});

						// Add refracted ray if needed
						if (first_interaction.y < cached_surface_y_ - 0.001f) {
							glm::vec4 refracted_color(0.9f, 0.9f, 1.0f, 0.8f);
							cached_line_instances_.push_back(
								{surface_entry, first_interaction, refracted_color, refracted_color});
						}
					}
				}

				while (current && next) {
					// Use gradient colors for better energy visualization
					float energy1 = static_cast<float>(current->value);
					float energy2 = static_cast<float>(next->value);

					glm::vec4 start_color = adaptive_log_color(energy1);
					glm::vec4 end_color = adaptive_log_color(energy2);
					start_color.a = 1.0f;
					end_color.a = 1.0f;

					glm::vec3 start(static_cast<float>(current->position.x),
									static_cast<float>(current->position.y),
									static_cast<float>(current->position.z));

					glm::vec3 end(static_cast<float>(next->position.x),
								  static_cast<float>(next->position.y),
								  static_cast<float>(next->position.z));

					// Line segment with gradient colors
					cached_line_instances_.push_back({start, end, start_color, end_color});

					// Check for emitter connections at EVERY node during path traversal
					// Check if current node has an emit connection (external vertex)
					if (current->emit && current->emit->emitter) {
						const auto emitter = current->emit->emitter;

						glm::vec3 scatter_pos(static_cast<float>(current->position.x),
											  static_cast<float>(current->position.y),
											  static_cast<float>(current->position.z));

						glm::vec3 exit_point(static_cast<float>(emitter->position.x),
											 static_cast<float>(emitter->position.y),
											 static_cast<float>(emitter->position.z));

						// Use current node energy for coloring the exit segment
						float exit_energy = static_cast<float>(current->value);
						glm::vec4 exit_line_color = adaptive_log_color(exit_energy);
						exit_line_color.a = 0.8f; // Slightly transparent to distinguish from main path

						cached_line_instances_.push_back({scatter_pos, exit_point, exit_line_color, exit_line_color});
					}

					// Also check direct emitter connection (fallback)
					if (current->emitter) {
						const auto emitter = current->emitter;

						glm::vec3 scatter_pos(static_cast<float>(current->position.x),
											  static_cast<float>(current->position.y),
											  static_cast<float>(current->position.z));

						glm::vec3 exit_point(static_cast<float>(emitter->position.x),
											 static_cast<float>(emitter->position.y),
											 static_cast<float>(emitter->position.z));

						// Use current node energy for coloring the exit segment
						float exit_energy = static_cast<float>(current->value);
						glm::vec4 exit_line_color = adaptive_log_color(exit_energy);
						exit_line_color.a = 0.8f; // Slightly transparent to distinguish from main path

						cached_line_instances_.push_back({scatter_pos, exit_point, exit_line_color, exit_line_color});
					}

					// Move to next segment
					current = next;
					next = current->next;
				}

				// Handle the final node (which doesn't have a 'next' but might have emitters)
				if (current) {
					// Check for emitter connections at the final node
					if (current->emit && current->emit->emitter) {
						const auto emitter = current->emit->emitter;

						glm::vec3 last_scatter(static_cast<float>(current->position.x),
											   static_cast<float>(current->position.y),
											   static_cast<float>(current->position.z));

						glm::vec3 exit_point(static_cast<float>(emitter->position.x),
											 static_cast<float>(emitter->position.y),
											 static_cast<float>(emitter->position.z));

						// Use last node energy for coloring the exit segment
						float exit_energy = static_cast<float>(current->value);
						glm::vec4 exit_line_color = adaptive_log_color(exit_energy);
						exit_line_color.a = 0.8f; // Slightly transparent to distinguish from main path

						cached_line_instances_.push_back({last_scatter, exit_point, exit_line_color, exit_line_color});
					}

					// Also check direct emitter connection (fallback) for final node
					if (current->emitter) {
						const auto emitter = current->emitter;

						glm::vec3 last_scatter(static_cast<float>(current->position.x),
											   static_cast<float>(current->position.y),
											   static_cast<float>(current->position.z));

						glm::vec3 exit_point(static_cast<float>(emitter->position.x),
											 static_cast<float>(emitter->position.y),
											 static_cast<float>(emitter->position.z));

						// Use last node energy for coloring the exit segment
						float exit_energy = static_cast<float>(current->value);
						glm::vec4 exit_line_color = adaptive_log_color(exit_energy);
						exit_line_color.a = 0.8f; // Slightly transparent to distinguish from main path

						cached_line_instances_.push_back({last_scatter, exit_point, exit_line_color, exit_line_color});
					}
				}
			}
		}

		// Add all emitter direction vectors to instanced rendering cache with deduplication
		if (!simulator_->get().emitters.empty()) {
			// Group emitters by origin position and direction (with margin for similar directions)
			std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::shared_ptr<Emitter>>> direction_groups;
			const double POSITION_PRECISION = 100.0; // Group positions within 0.01 units
			const double DIRECTION_PRECISION = 50.0; // Group directions within ~0.02 radians (~1.1 degrees)

			for (const auto& emitter : simulator_->get().emitters) {
				// Only process significant emitters
				if (emitter->weight > 0.001) {
					// Quantize position and direction for grouping
					int pos_x = static_cast<int>(std::round(emitter->position.x * POSITION_PRECISION));
					int pos_y = static_cast<int>(std::round(emitter->position.y * POSITION_PRECISION));
					int pos_z = static_cast<int>(std::round(emitter->position.z * POSITION_PRECISION));
					int dir_x = static_cast<int>(std::round(emitter->direction.x * DIRECTION_PRECISION));
					int dir_y = static_cast<int>(std::round(emitter->direction.y * DIRECTION_PRECISION));
					int dir_z = static_cast<int>(std::round(emitter->direction.z * DIRECTION_PRECISION));

					auto group_key = std::make_tuple(pos_x, pos_y, pos_z, dir_x, dir_y, dir_z);
					direction_groups[group_key].push_back(emitter);
				}
			}

			// Render one emitter vector per group with averaged energy
			for (const auto& [group_key, grouped_emitters] : direction_groups) {
				if (grouped_emitters.empty()) {
					continue;
				}

				// Use the first emitter for position and direction
				const auto& representative = grouped_emitters[0];
				glm::vec3 start_pos(static_cast<float>(representative->position.x),
									static_cast<float>(representative->position.y),
									static_cast<float>(representative->position.z));

				glm::vec3 end_pos(static_cast<float>(representative->position.x + representative->direction.x * 0.05),
								  static_cast<float>(representative->position.y + representative->direction.y * 0.05),
								  static_cast<float>(representative->position.z + representative->direction.z * 0.05));

				// Average the energy of all emitters in this group
				double total_weight = 0.0;
				for (const auto& emitter : grouped_emitters) {
					total_weight += emitter->weight;
				}

				// Use percentage-based coloring with averaged energy
				float surface_refraction = static_cast<float>(simulator_->get().get_combined_surface_refraction());
				float energy_percentage = static_cast<float>(total_weight / surface_refraction);
				glm::vec4 direction_color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
				direction_color.a = 0.8f; // Slightly transparent for distinction

				cached_line_instances_.push_back({start_pos, end_pos, direction_color, direction_color});
			}
		}

		// Add specular reflection emitter for incident photon's surface reflection
		// This is integrated into the direction grouping above to avoid duplication
		double specular_reflection = simulator_->get().get_combined_specular_reflection();
		if (specular_reflection > 0.0 && !simulator_->get().sources.empty()) {
			const Source& source = simulator_->get().sources[0];

			// Check if there are any regular emitters at the same position with similar direction
			bool found_similar_emitter = false;
			glm::dvec3 specular_pos = source.intersect;
			glm::dvec3 specular_dir = glm::normalize(source.specular_direction);

			const double POSITION_TOLERANCE = 0.01; // Same as POSITION_PRECISION above
			const double ANGLE_TOLERANCE = 0.02;    // Same as DIRECTION_PRECISION above

			for (const auto& emitter : simulator_->get().emitters) {
				if (emitter->weight <= 0.001) {
					continue;
				}

				double pos_distance = glm::length(emitter->position - specular_pos);
				glm::dvec3 emitter_dir = glm::normalize(emitter->direction);
				double angle_diff = glm::length(emitter_dir - specular_dir);

				if (pos_distance <= POSITION_TOLERANCE && angle_diff <= ANGLE_TOLERANCE) {
					found_similar_emitter = true;
					break;
				}
			}

			// Only add specular reflection vector if no similar emitter exists
			if (!found_similar_emitter) {
				glm::vec3 reflection_start(static_cast<float>(source.intersect.x),
										   static_cast<float>(source.intersect.y),
										   static_cast<float>(source.intersect.z));

				// Make specular reflection vector twice as long as regular emitters (0.1 vs 0.05)
				glm::vec3 reflection_end(static_cast<float>(source.intersect.x + source.specular_direction.x * 0.1),
										 static_cast<float>(source.intersect.y + source.specular_direction.y * 0.1),
										 static_cast<float>(source.intersect.z + source.specular_direction.z * 0.1));

				// Use energy-based coloring - map specular reflection energy to color
				// Use the actual specular reflection value as energy for better color mapping
				float specular_energy = static_cast<float>(specular_reflection);
				glm::vec4 specular_color = get_adaptive_energy_color(specular_energy, 0.0f, 1.0f);
				specular_color.a = 1.0f; // Full opacity for better visibility

				cached_line_instances_.push_back({reflection_start, reflection_end, specular_color, specular_color});
			}
		}

		// Collect scatter points inside cache block to eliminate per-frame processing
		// First, collect scatter points from photon paths
		for (const Photon& photon : simulator_->get().photons) {
			if (photon.path_head) {
				auto current = photon.path_head;
				int vertex_count = 0;

				// Count total vertices to identify key points properly
				auto temp = current;
				while (temp) {
					vertex_count++;
					temp = temp->next;
				}

				// Only add markers at specific key points: incident, scatter, exit
				auto path_current = photon.path_head;
				std::shared_ptr<PhotonNode> prev = nullptr;
				std::shared_ptr<PhotonNode> next = nullptr;
				int current_index = 0;

				while (path_current) {
					glm::vec3 pos(static_cast<float>(path_current->position.x),
								  static_cast<float>(path_current->position.y),
								  static_cast<float>(path_current->position.z));

					bool should_mark = false;
					glm::vec4 marker_color {};

					if (current_index == 0) {
						// First vertex - incident point (bright white)
						should_mark = true;
						marker_color = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
					}
					else if (current_index == vertex_count - 1) {
						// Last vertex - exit point with adaptive energy-based coloring
						should_mark = true;
						float energy = static_cast<float>(path_current->value);
						marker_color = adaptive_log_color(energy);
					}
					else if (current_index > 0 && current_index < vertex_count - 1) {
						// Check for medium boundary crossings and path splits
						next = path_current->next;
						if (prev && next) {
							// Get positions
							glm::vec3 prev_pos(static_cast<float>(prev->position.x),
											   static_cast<float>(prev->position.y),
											   static_cast<float>(prev->position.z));

							glm::vec3 next_pos(static_cast<float>(next->position.x),
											   static_cast<float>(next->position.y),
											   static_cast<float>(next->position.z));

							// Check if this point represents a medium boundary crossing
							bool is_medium_boundary = false;

							// Surface entry/exit detection (z-coordinate near 0)
							if (std::abs(pos.z) < 0.001f && std::abs(prev_pos.z) > 0.001f) {
								is_medium_boundary = true; // Entry into medium
							}
							else if (std::abs(pos.z) < 0.001f && std::abs(next_pos.z) > 0.001f) {
								is_medium_boundary = true; // Exit from medium
							}

							// Check for path splits (if this vertex has emitted paths)
							bool has_emit = (path_current->emit != nullptr);

							if (is_medium_boundary || has_emit) {
								should_mark = true;
								// Use adaptive energy-based coloring for all boundary/split points
								float energy = static_cast<float>(path_current->value);
								marker_color = adaptive_log_color(energy);
							}
						}
					}

					if (should_mark) {
						// Add scatter points to point instances
						PointInstance point_instance;
						point_instance.position = pos;
						point_instance.color = marker_color;
						point_instance.size = 6.0f; // Smaller size for scatter points
						cached_point_instances_.push_back(point_instance);
					}

					prev = path_current;
					path_current = path_current->next;
					current_index++;
				}
			}
		}

		// Collect emitter points inside cache block
		// Second, add emitter exit points (these are the accurate surface boundary points)
		if (!simulator_->get().emitters.empty()) {
			// Add emitter points to the point instances vector
			for (const auto& emitter : simulator_->get().emitters) {
				// Use emitter->position which contains the corrected surface intersection coordinates
				glm::vec3 exit_pos(static_cast<float>(emitter->position.x),
								   static_cast<float>(emitter->position.y),
								   static_cast<float>(emitter->position.z));

				// Use percentage-based coloring instead of absolute weights for consistency
				// Convert absolute weight to percentage of total energy budget
				float surface_refraction = static_cast<float>(simulator_->get().get_combined_surface_refraction());
				float energy_percentage = static_cast<float>(emitter->weight / surface_refraction);
				glm::vec4 exit_color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
				exit_color.a = 1.0f; // Full opacity for emitter points

				// Add emitter points to point instances
				PointInstance point_instance;
				point_instance.position = exit_pos;
				point_instance.color = exit_color;
				point_instance.size = 6.0f; // Same size as scatter points
				cached_point_instances_.push_back(point_instance);
			}

			// Add surface specular reflection as an emitter point
			double surface_specular_reflection = simulator_->get().get_combined_specular_reflection();
			if (surface_specular_reflection > 0.0 && !simulator_->get().sources.empty()) {
				const Source& source = simulator_->get().sources[0];

				// Use actual source intersection point
				glm::vec3 surface_entry(static_cast<float>(source.intersect.x),
										static_cast<float>(source.intersect.y),
										static_cast<float>(source.intersect.z));

				// Use percentage-based coloring for surface reflection
				double surface_refraction = simulator_->get().get_combined_surface_refraction();
				float energy_percentage = static_cast<float>(surface_specular_reflection / surface_refraction);
				glm::vec4 surface_point_color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
				surface_point_color.a = 1.0f; // Full opacity for surface reflection point

				PointInstance surface_point;
				surface_point.position = surface_entry;
				surface_point.color = surface_point_color;
				surface_point.size = 8.0f;    // Slightly larger for surface reflection point
				cached_point_instances_.push_back(surface_point);
			}
		}

		// Update cached photon count and mark cache as up to date
		cached_photon_count_ = current_photon_count;
		path_instances_cached_ = true;
	}

	// Render cached line instances
	if (!cached_line_instances_.empty() && lines_shader_) {
		glUseProgram(lines_shader_);

		// Use cached uniform location, compute MVP directly
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		glUniformMatrix4fv(line_instanced_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

		// Only upload buffer when data has changed, not every frame
		if (!line_buffer_uploaded_) {
			glBindBuffer(GL_ARRAY_BUFFER, lines_instance_vbo_);
			glBufferData(GL_ARRAY_BUFFER,
						 cached_line_instances_.size() * sizeof(LineInstance),
						 cached_line_instances_.data(),
						 GL_STATIC_DRAW);
			line_buffer_uploaded_ = true;
		}

		// Render using uploaded buffer
		glBindVertexArray(lines_vao_);
		glDrawArraysInstanced(GL_LINES, 0, 2, static_cast<GLsizei>(cached_line_instances_.size()));
		glBindVertexArray(0);

		glUseProgram(0);
	}

	// Render cached scatter points and emitter points using instanced rendering
	if (!cached_point_instances_.empty() && points_shader_) {
		// Set up all required OpenGL state for point rendering
		// Ensures points render correctly independent of other geometry
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_DEPTH_TEST); // Ensure depth testing is enabled
		glDisable(GL_CULL_FACE); // Ensure face culling is disabled for points
		enable_blending();       // Use state management helper
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// Ensure shader program is properly bound (fix for volume geometry dependency)
		glUseProgram(points_shader_);
		current_shader_program_ = points_shader_;

		// Use cached uniform location, compute MVP directly
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		glUniformMatrix4fv(point_instanced_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

		// Only upload buffer when data has changed, not every frame
		if (!point_buffer_uploaded_) {
			glBindBuffer(GL_ARRAY_BUFFER, points_instance_vbo_);
			glBufferData(GL_ARRAY_BUFFER,
						 cached_point_instances_.size() * sizeof(PointInstance),
						 cached_point_instances_.data(),
						 GL_STATIC_DRAW);
			point_buffer_uploaded_ = true;
		}

		// Bind VAO and render using uploaded buffer
		glBindVertexArray(points_vao_);
		glDrawArraysInstanced(GL_POINTS, 0, 1, static_cast<GLsizei>(cached_point_instances_.size()));
		glBindVertexArray(0);

		glUseProgram(0);

		// Restore OpenGL state once after rendering
		disable_blending(); // Use state management helper
		glDisable(GL_PROGRAM_POINT_SIZE);
	}
}

// ========================================
// GEOMETRY SETUP HELPER FUNCTIONS FOR INSTANCED RENDERING
// ========================================

static bool setup_point_geometry(GLuint vbo) {
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	float point_vertex[3] = {0.0f, 0.0f, 0.0f}; // Single point at origin
	glBufferData(GL_ARRAY_BUFFER, sizeof(point_vertex), point_vertex, GL_STATIC_DRAW);

	// Vertex position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	return true;
}

static bool setup_line_geometry(GLuint vbo) {
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	// clang-format off
	float line_vertices[] = {
		0.0f, 0.0f, 0.0f, // Start point
		1.0f, 0.0f, 0.0f  // End point
	};
	// clang-format on

	glBufferData(GL_ARRAY_BUFFER, sizeof(line_vertices), line_vertices, GL_STATIC_DRAW);

	// Vertex position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	return true;
}

static bool setup_voxel_geometry(GLuint vbo) {
	// Create unit cube geometry (centered at origin)
	// Unit cube vertices: position (3) + normal (3) per vertex

	// clang-format off
	float vertices[] = {
		// Front face
		-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
		 0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
		 0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
		 0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
		-0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,
		-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,

		// Back face
		-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
		 0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
		 0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
		 0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
		-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
		-0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

		// Left face
		-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
		-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
		-0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
		-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
		-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
		-0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,

		// Right face
		 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
		 0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
		 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
		 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
		 0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
		 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,

		// Bottom face
		-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
		 0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
		 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
		 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
		-0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
		-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,

		// Top face
		-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
		 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
		 0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
		 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
		-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
		-0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f
	};
	// clang-format on

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Normal attribute (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	return true;
}

// ========================================
// CONSOLIDATED SHADER-BASED RENDERING METHODS
// ========================================

GLuint Renderer::create_shader_program(const std::string& vertex_source, const std::string& fragment_source) {
	GLuint vertex_shader = Shader::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = Shader::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader) {
			glDeleteShader(vertex_shader);
		}
		if (fragment_shader) {
			glDeleteShader(fragment_shader);
		}
		return 0;
	}

	GLuint program = glCreateProgram();
	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	GLint success;
	glGetProgramiv(program, GL_LINK_STATUS, &success);

	if (!success) {
		std::array<char, 512> info_log {};
		glGetProgramInfoLog(program, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Program linking failed: " << info_log.data() << std::endl;
		glDeleteProgram(program);
		program = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	return program;
}

std::string Renderer::load_shader_source(const std::string& file_path) {
	return Shader::load_shader_source(file_path);
}

void Renderer::auto_manage_energy_labels(Settings& settings) {
	if (!simulator_) {
		return;
	}

	static bool auto_disabled_labels = false;
	bool many_photons = (simulator_->get().photons.size() > 10);

	// Auto-disable when crossing the 10 photon threshold
	if (many_photons && !auto_disabled_labels) {
		settings.draw_labels = false;
		auto_disabled_labels = true;
	}
	// Reset auto-disable flag when photon count drops back down
	else if (!many_photons && auto_disabled_labels) {
		auto_disabled_labels = false;
	}
}

void Renderer::cache_energy_labels() {
	cached_energy_labels_.clear();
	energy_labels_cached_ = false;

	if (!simulator_) {
		return;
	}

	// Use emitter data for accurate energy labels with proper classification
	const auto& emitters = simulator_->get().emitters;

	if (emitters.empty()) {
		energy_labels_cached_ = true;
		return;
	}

	// ADAPTIVE DENSITY-BASED GROUPING
	// Calculate total photon count for percentage-based thresholds
	int total_photons = static_cast<int>(emitters.size());
	if (total_photons == 0) {
		energy_labels_cached_ = true;
		return;
	}

	// First pass: Calculate local density (photons per unit area) for each emitter
	std::map<std::shared_ptr<Emitter>, double> emitter_density;
	const double DENSITY_RADIUS = 0.1;                                               // Radius for density calculation
	const double DENSITY_AREA = MathConstants::PI * DENSITY_RADIUS * DENSITY_RADIUS; // Circle area

	for (const auto& emitter : emitters) {
		// Skip very low energy exits
		if (emitter->weight < 0.001) {
			continue;
		}

		int nearby_count = 0;
		for (const auto& other : emitters) {
			if (other->weight < 0.001) {
				continue;
			}

			double distance = glm::length(emitter->position - other->position);
			if (distance <= DENSITY_RADIUS) {
				nearby_count++;
			}
		}

		// Calculate density as photons per unit area
		emitter_density[emitter] = static_cast<double>(nearby_count) / DENSITY_AREA;
	}

	// Calculate density percentiles for adaptive thresholds
	std::vector<double> all_densities;
	for (const auto& [emitter, density] : emitter_density) {
		all_densities.push_back(density);
	}
	std::sort(all_densities.begin(), all_densities.end());

	double density_90th = all_densities[static_cast<size_t>(all_densities.size() * 0.9)];
	double density_75th = all_densities[static_cast<size_t>(all_densities.size() * 0.75)];
	double density_50th = all_densities[static_cast<size_t>(all_densities.size() * 0.5)];

	// Second pass: Group emitters by position using adaptive grouping radius
	std::map<std::tuple<int, int, int>, std::vector<std::shared_ptr<Emitter>>> position_groups;

	for (const auto& emitter : emitters) {
		if (emitter->weight < 0.001)
			continue; // Skip very low energy exits

		// Determine grouping multiplier based on density percentile
		double density = emitter_density[emitter];
		double grouping_multiplier;

		if (density >= density_90th) {
			grouping_multiplier = 5.0;  // Large groups for top 10% density (0.2 units)
		}
		else if (density >= density_75th) {
			grouping_multiplier = 10.0; // Medium groups for top 25% density (0.1 units)
		}
		else if (density >= density_50th) {
			grouping_multiplier = 20.0; // Small groups for above median density (0.05 units)
		}
		else {
			grouping_multiplier = 50.0; // Fine granularity for below median density (0.02 units)
		}

		// Group by position using adaptive precision
		int pos_x = static_cast<int>(std::round(emitter->position.x * grouping_multiplier));
		int pos_y = static_cast<int>(std::round(emitter->position.y * grouping_multiplier));
		int pos_z = static_cast<int>(std::round(emitter->position.z * grouping_multiplier));
		auto pos_key = std::make_tuple(pos_x, pos_y, pos_z);

		position_groups[pos_key].push_back(emitter);
	}

	// Create labels for each position group
	for (const auto& [pos_key, emitters_at_pos] : position_groups) {
		if (emitters_at_pos.empty())
			continue;

		// Use the first emitter's position as the label position
		const auto& representative = emitters_at_pos[0];
		glm::vec3 label_pos(static_cast<float>(representative->position.x),
							static_cast<float>(representative->position.y),
							static_cast<float>(representative->position.z));

		// Sum up energy and classify by exit type
		double total_reflected_energy = 0.0;
		double total_transmitted_energy = 0.0;
		double total_unclassified_energy = 0.0;

		for (const auto& emitter : emitters_at_pos) {
			switch (emitter->exit_type) {
				case Emitter::ExitType::REFLECTED: total_reflected_energy += emitter->weight; break;
				case Emitter::ExitType::TRANSMITTED: total_transmitted_energy += emitter->weight; break;
				default: total_unclassified_energy += emitter->weight; break;
			}
		}

		// Create label with proper classification and NORMALIZED percentages
		std::string label_text;
		glm::vec4 label_color(1.0f, 1.0f, 1.0f, 1.0f);

		double total_energy = total_reflected_energy + total_transmitted_energy + total_unclassified_energy;

		// Normalize by total number of photons, not surface refraction
		// Each photon starts with weight 1.0, so total initial energy = num_photons
		double total_initial_energy = static_cast<double>(simulator_->get().photons.size());

		// Calculate normalized percentage (matches console energy conservation calculation)
		double energy_percent = (total_energy / total_initial_energy) * 100.0;

		// Format percentage, showing "<1%" instead of "0%" for very small values
		int rounded_percent = static_cast<int>(energy_percent);
		std::string percent_text =
			(rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);

		if (total_reflected_energy > total_transmitted_energy && total_reflected_energy > total_unclassified_energy) {
			// Predominantly reflected
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.8f, 0.2f, 1.0f); // Bright green for reflection (clearly visible)
		}
		else if (total_transmitted_energy > total_reflected_energy
				 && total_transmitted_energy > total_unclassified_energy) {
			// Predominantly transmitted
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.6f, 1.0f, 1.0f); // Bright blue for transmission (clearly visible)
		}
		else {
			// Mixed or unclassified
			label_text = percent_text;
			label_color = glm::vec4(0.8f, 0.8f, 0.8f, 1.0f); // Gray for mixed/unclassified
		}

		cached_energy_labels_.push_back({.world_position = label_pos,
										 .text = label_text,
										 .color = label_color,
										 .scale = 1.0f,
										 .screen_position = glm::vec2(0.0f),
										 .screen_position_valid = false});
	}

	// Add surface specular reflection label - position at tip of reflection vector
	double specular_reflection = simulator_->get().get_combined_specular_reflection();
	if (specular_reflection > 0.0 && !simulator_->get().sources.empty()) {
		const Source& source = simulator_->get().sources[0];

		// Position label at the tip of the specular reflection vector (twice as long as regular emitters)
		glm::vec3 surface_label_pos(static_cast<float>(source.intersect.x + source.specular_direction.x * 0.1),
									static_cast<float>(source.intersect.y + source.specular_direction.y * 0.1),
									static_cast<float>(source.intersect.z + source.specular_direction.z * 0.1));

		// Calculate normalized percentage for surface reflection
		double surface_refraction = simulator_->get().get_combined_surface_refraction();
		double energy_percent = (specular_reflection / surface_refraction) * 100.0;

		// Format percentage, showing "<1%" instead of "0%" for very small values
		int rounded_percent = static_cast<int>(energy_percent);
		std::string surface_percent_text =
			(rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);

		// Use bright purple color for surface specular reflection to match saturation of green and blue
		glm::vec4 surface_label_color = glm::vec4(0.8f, 0.3f, 1.0f, 1.0f);

		cached_energy_labels_.push_back({.world_position = surface_label_pos,
										 .text = surface_percent_text,
										 .color = surface_label_color,
										 .scale = 1.0f,
										 .screen_position = glm::vec2(0.0f),
										 .screen_position_valid = false});
	}

	energy_labels_cached_ = true;
}

void Renderer::invalidate_energy_label_cache() {
	energy_labels_cached_ = false;
	cached_energy_labels_.clear();
	camera_state_changed_ = true; // Force screen position recalculation
}

void Renderer::update_energy_label_screen_positions() {
	if (!energy_labels_cached_ || cached_energy_labels_.empty()) {
		return;
	}

	// Check if camera has actually changed
	glm::vec3 current_position = camera_.get_position();
	// Update screen positions only when camera changes (simplified without detailed tracking)
	if (camera_state_changed_) {
		camera_state_changed_ = false;

		// Update all cached screen positions
		for (auto& label : cached_energy_labels_) {
			label.screen_position = world_to_screen(label.world_position, viewport_width_, viewport_height_);
			label.screen_position_valid =
				(label.screen_position.x >= 0 && label.screen_position.x < viewport_width_
				 && label.screen_position.y >= 0 && label.screen_position.y < viewport_height_);
		}
	}
}

void Renderer::draw_labels(const Settings& settings) {
	if (!simulator_ || !settings.draw_labels)
		return;

	// Cache energy labels if not already cached
	if (!energy_labels_cached_) {
		cache_energy_labels();
	}

	// Render cached text labels using pre-calculated screen positions
	if (text_render_callback_) {
		for (const auto& label : cached_energy_labels_) {
			if (label.screen_position_valid) {
				text_render_callback_(label.text, label.screen_position.x, label.screen_position.y, label.color);
			}
		}
	}
}

glm::vec4 Renderer::get_adaptive_energy_color(float energy, float min_energy, float max_energy) {
	// Enhanced non-linear normalization with wider color range and better high-energy distinction
	// Most voxels have very low energy, so we need to expand that range visually
	// while providing clear distinction at high energies (80%-100%)

	// Clamp energy to valid range
	energy = std::clamp(energy, min_energy, max_energy);

	// First do linear normalization
	float linear_normalized = (energy - min_energy) / (max_energy - min_energy);
	linear_normalized = std::clamp(linear_normalized, 0.0f, 1.0f);

	// Apply power function to expand low-energy visualization
	// Using power of 0.25 for even better low-energy expansion
	float normalized = std::pow(linear_normalized, 0.25f);
	normalized = std::clamp(normalized, 0.0f, 1.0f);

	// Enhanced color zones with better high-energy distinction and clearer visual hierarchy
	// Color progression: brilliant blue-white > white > bright yellow > yellow > orange > red > dark red

	if (normalized > 0.95f) {
		// Very high energy (95-100%): brilliant white with slight blue tint
		float t = (normalized - 0.95f) / 0.05f;
		return glm::vec4(1.0f, 1.0f, 1.0f + t * 0.2f, 1.0f); // Slightly blue-white
	}
	else if (normalized > 0.90f) {
		// High energy (90-95%): pure white
		return glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	}
	else if (normalized > 0.85f) {
		// High-medium energy (85-90%): white to very bright yellow
		float t = (normalized - 0.85f) / 0.05f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 1.0f - t * 0.4f; // From white (1.0) to very bright yellow (0.6)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.80f) {
		// 80-85% energy: very bright yellow to bright yellow
		float t = (normalized - 0.80f) / 0.05f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.6f - t * 0.2f; // From very bright yellow (0.6) to bright yellow (0.4)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.70f) {
		// 70-80% energy: bright yellow to yellow
		float t = (normalized - 0.70f) / 0.10f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.4f - t * 0.2f; // From bright yellow (0.4) to yellow (0.2)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.55f) {
		// Medium-high energy: yellow to warmer yellow
		float t = (normalized - 0.55f) / 0.15f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.2f - t * 0.1f; // From yellow (0.2) to warmer yellow (0.1)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.40f) {
		// Medium energy: warmer yellow to orange
		float t = (normalized - 0.40f) / 0.15f;
		float r = 1.0f;
		float g = 1.0f - t * 0.3f; // From 1.0 (yellow) to 0.7 (orange)
		float b = 0.1f - t * 0.1f; // From warmer yellow (0.1) to orange (0.0)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.25f) {
		// Medium-low energy: orange to red-orange
		float t = (normalized - 0.25f) / 0.15f;
		float r = 1.0f;
		float g = 0.7f - t * 0.4f; // From 0.7 (orange) to 0.3 (red-orange)
		float b = 0.0f;
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.15f) {
		// Low energy: red-orange to red
		float t = (normalized - 0.15f) / 0.10f;
		float r = 1.0f;
		float g = 0.3f - t * 0.3f; // From 0.3 (red-orange) to 0.0 (red)
		float b = 0.0f;
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.08f) {
		// Low-medium red: bright red to medium red
		float t = (normalized - 0.08f) / 0.07f;
		float r = 1.0f - t * 0.15f; // From 1.0 (bright red) to 0.85 (medium red)
		float g = 0.0f;
		float b = 0.0f;
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.04f) {
		// Medium red: medium red to darker red
		float t = (normalized - 0.04f) / 0.04f;
		float r = 0.85f - t * 0.15f; // From 0.85 (medium red) to 0.7 (darker red)
		float g = 0.0f;
		float b = 0.0f;
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.02f) {
		// Dark red: darker red to dark red
		float t = (normalized - 0.02f) / 0.02f;
		float r = 0.7f - t * 0.15f; // From 0.7 (darker red) to 0.55 (dark red)
		float g = 0.0f;
		float b = 0.0f;
		return glm::vec4(r, g, b, 1.0f);
	}
	else {
		// Very dark red: darkest visible red
		return glm::vec4(0.55f, 0.0f, 0.0f, 1.0f);
	}
}

glm::vec4 Renderer::layer_energy_color(float energy, float min_energy, float max_energy, uint8_t layer_id) {
	// Clamp and normalize energy like the original function
	energy = std::clamp(energy, min_energy, max_energy);

	float linear_normalized = (energy - min_energy) / (max_energy - min_energy);
	linear_normalized = std::clamp(linear_normalized, 0.0f, 1.0f);

	float normalized = std::pow(linear_normalized, 0.3f);
	normalized = std::clamp(normalized, 0.0f, 1.0f);

	// Define base colors for different material types
	glm::vec3 base_colors[6] = {
		glm::vec3(1.0f, 0.4f, 0.4f),                  // Red theme for layer 0
		glm::vec3(0.4f, 0.4f, 1.0f),                  // Blue theme for layer 1
		glm::vec3(0.4f, 1.0f, 0.4f),                  // Green theme for layer 2
		glm::vec3(1.0f, 0.4f, 1.0f),                  // Magenta theme for layer 3
		glm::vec3(1.0f, 1.0f, 0.4f),                  // Yellow theme for layer 4
		glm::vec3(0.4f, 1.0f, 1.0f),                  // Cyan theme for layer 5
	};

	glm::vec3 base_color = base_colors[layer_id % 6]; // Create energy gradient within the material's color theme
	if (normalized > 0.85f) {
		// Very high energy: brighten to near white
		return glm::vec4(glm::mix(base_color, glm::vec3(1.0f), 0.6f), 1.0f);
	}
	else if (normalized > 0.70f) {
		// High energy: bright version of base color
		float t = (normalized - 0.70f) / 0.15f;
		glm::vec3 bright_color = base_color * 1.4f;             // Brighten
		bright_color = glm::min(bright_color, glm::vec3(1.0f)); // Clamp to valid range
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

glm::vec2 Renderer::world_to_screen(const glm::vec3& world_pos, int screen_width, int screen_height) const {
	// Use the actual camera matrices that are being used for rendering
	glm::mat4 mvp = camera_.get_mvp_matrix();

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

// ========================================
// Instanced Voxel Rendering Implementation
// ========================================

void Renderer::end_voxel_instances(VoxelMode mode) {
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
			cached_camera_position_ = camera_.get_position();
			cached_camera_target_ = camera_.get_target();
		}

		// Check if we need to start a new background sort
		if (!background_sort_in_progress_.load()) {
			glm::vec3 current_camera_pos = camera_.get_position();
			glm::vec3 current_camera_target = camera_.get_target();

			// Only start background sort if camera moved significantly
			const float SIGNIFICANT_MOVEMENT = 0.05f; // Threshold for starting background sort
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

	// Balance performance vs quality
	// Use consistent sorting behavior for all modes
	bool should_sort_for_quality = false;

	if (voxel_instances_.size() < 50000) {
		should_sort_for_quality = true; // Sort small datasets
	}

	if (should_sort_for_quality) {
		// Sort by depth (back to front) for optimal transparency
		std::sort(voxel_instances_.begin(), voxel_instances_.end(), [](const VoxelInstance& a, const VoxelInstance& b) {
			return a.depth > b.depth; // Back to front
		});
	}

	// Upload instance data to GPU
	glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);
	glBufferData(GL_ARRAY_BUFFER,
				 voxel_instances_.size() * sizeof(VoxelInstance),
				 voxel_instances_.data(),
				 GL_STATIC_DRAW); // Use STATIC_DRAW for better performance

	voxel_buffer_uploaded_ = true;
	voxel_instances_dirty_ = false;

	// Cache current camera position for change detection
	cached_camera_position_ = camera_.get_position();
	cached_camera_target_ = camera_.get_target();
}
void Renderer::draw_voxel_instances() {
	if (voxel_instances_.empty() || !voxel_shader_) {
		return;
	}

	// Enable transparency - KEEP ORIGINAL APPEARANCE
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Disable depth writing for transparent objects but keep depth testing
	glDepthMask(GL_FALSE);
	glEnable(GL_DEPTH_TEST);

	glUseProgram(voxel_shader_);

	// Set MVP matrix using cached location
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(voxel_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

	// Bind VAO and draw instanced
	glBindVertexArray(voxel_vao_);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(voxel_instances_.size()));

	// Restore OpenGL state
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glBindVertexArray(0);
}

bool Renderer::setup_voxel_instanced_rendering() {
	using namespace RenderSetup;

	ShaderConfig shader_config {"shaders/voxels.vert", "shaders/voxels.frag", "voxel instanced"};

	std::vector<UniformConfig> uniforms {{"uMVP", &voxel_mvp_uniform_location_}};

	std::vector<VertexAttributeConfig> instance_attributes {
		{2, 3, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, position), 1}, // Instance
																										// position
		{3, 4, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, color), 1}, // Instance color
		{4, 1, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, scale), 1}  // Instance scale
	};

	return setup_instanced_rendering<VoxelInstance>(shader_config,
													uniforms,
													setup_voxel_geometry,
													instance_attributes,
													voxel_shader_,
													voxel_vao_,
													voxel_vbo_,
													voxel_instance_vbo_);
}

bool Renderer::setup_point_instanced_rendering() {
	using namespace RenderSetup;

	ShaderConfig shader_config {"shaders/points.vert", "shaders/points.frag", "point instanced"};
	std::vector<UniformConfig> uniforms {{"uMVP", &point_instanced_mvp_uniform_location_}};

	std::vector<VertexAttributeConfig> instance_attributes {
		{1, 3, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, position), 1}, // Instance
																										// position
		{2, 4, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, color), 1}, // Instance color
		{3, 1, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, size), 1}   // Instance size
	};

	return setup_instanced_rendering<PointInstance>(shader_config,
													uniforms,
													setup_point_geometry,
													instance_attributes,
													points_shader_,
													points_vao_,
													points_vbo_,
													points_instance_vbo_);
}

// Cache energy range calculation to avoid expensive per-frame analysis
void Renderer::update_cached_energy_range(const Settings& settings) const {
	if (!simulator_) {
		energy_range_cached_ = false;
		return;
	}

	// Use proper cache invalidation instead of checking every frame
	if (energy_range_cached_) {
		return; // Use cached values - they're already invalidated by simulation version changes
	}

	// Recalculate energy range using optimized method
	std::vector<float> all_energies;

	// Collect energies from photon paths
	for (const Photon& photon : simulator_->get().photons) {
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

	// Collect energies from voxels using optimized single-loop method
	const auto& config = Config::get();
	const uint32_t nx = static_cast<uint32_t>(config.nx());
	const uint32_t ny = static_cast<uint32_t>(config.ny());
	const uint32_t nz = static_cast<uint32_t>(config.nz());

	if (nx > 0 && ny > 0 && nz > 0) {
		const uint32_t total_voxels = nx * ny * nz;

		// Pre-determine energy calculation method outside the loop for better performance
		const bool use_absorption = (settings.voxel_mode == VoxelMode::Absorption);
		const bool use_emittance = (settings.voxel_mode == VoxelMode::Emittance);
		const bool use_layers = (settings.voxel_mode == VoxelMode::Layers);

		// Reserve space to avoid repeated allocations (estimate ~10% of voxels have energy)
		all_energies.reserve(all_energies.size() + total_voxels / 10);

		// Single loop with linear indexing for better cache performance
		for (uint32_t linear_idx = 0; linear_idx < total_voxels; ++linear_idx) {
			// Convert linear index to 3D coordinates (z-major order for cache efficiency)
			const uint32_t iz = linear_idx / (nx * ny);
			const uint32_t iy = (linear_idx % (nx * ny)) / nx;
			const uint32_t ix = linear_idx % nx;

			Voxel* voxel = simulator_->get().voxel_grid(ix, iy, iz);
			if (voxel) {
				float total_energy;

				if (use_absorption) {
					total_energy = static_cast<float>(voxel->absorption);
				}
				else if (use_emittance) {
					total_energy = static_cast<float>(voxel->total_emittance());
				}
				else if (use_layers) {
					total_energy = (voxel->material != nullptr) ? 1.0f : 0.0f;
				}
				else {
					// Combined mode: calculate both values only when needed
					total_energy = static_cast<float>(voxel->absorption + voxel->total_emittance());
				}

				if (total_energy > 0.0000001f) {
					all_energies.push_back(total_energy);
				}
			}
		}
	}

	// Calculate range using percentile-based method for accurate colors
	cached_min_energy_ = MathConstants::ENERGY_THRESHOLD; // Much lower threshold for emittance detection
	cached_max_energy_ = 0.01f;                           // Default fallback

	if (!all_energies.empty()) {
		// Use percentile-based scaling for better contrast
		std::ranges::sort(all_energies);

		// Use percentile-based scaling for better contrast
		const size_t count = all_energies.size();
		const float p5 = all_energies[static_cast<size_t>(count * 0.05)];   // 5th percentile
		const float p95 = all_energies[static_cast<size_t>(count * 0.95)];  // 95th percentile

		cached_min_energy_ = std::max(p5, MathConstants::ENERGY_THRESHOLD); // Lower minimum
		cached_max_energy_ = std::max(p95, cached_min_energy_ * 10.0f);     // Ensure reasonable range
	}

	// Update cache state
	energy_range_cached_ = true;
}

// Comprehensive cache invalidation
void Renderer::invalidate_all_caches() {
	// Invalidate energy-related caches
	energy_range_cached_ = false;
	invalidate_energy_label_cache();

	// Invalidate path instance caches
	invalidate_path_instances_cache();

	// Invalidate voxel instance caches
	voxel_instances_dirty_ = true;
	voxel_buffer_uploaded_ = false;

	// Invalidate surface geometry cache
	surface_cached_ = false;
}

void Renderer::invalidate_path_instances_cache() {
	path_instances_cached_ = false;
	cached_line_instances_.clear();
	cached_point_instances_.clear();
	cached_photon_count_ = 0; // Reset photon count for complete rebuild

	// Reset buffer upload flags when cache is invalidated
	line_buffer_uploaded_ = false;
	point_buffer_uploaded_ = false;

	// Also reset dynamic geometry buffer flags
	point_geometry_buffer_uploaded_ = false;
}

// Only invalidate path cache when path-related settings change
void Renderer::set_settings(const Settings& new_settings) {
	// Check if path-related settings have changed
	bool path_settings_changed = (settings_.draw_paths != new_settings.draw_paths);

	// Update settings
	settings_ = new_settings;

	// Only invalidate path cache if path visualization settings changed
	if (path_settings_changed) {
		invalidate_path_instances_cache();
	}
}

bool Renderer::setup_line_instanced_rendering() {
	using namespace RenderSetup;

	ShaderConfig shader_config {"shaders/lines.vert", "shaders/lines.frag", "line instanced"};
	std::vector<UniformConfig> uniforms {{"uMVP", &line_instanced_mvp_uniform_location_}};

	std::vector<VertexAttributeConfig> instance_attributes {
		// Instance start position
		{1, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start), 1},

		// Instance end position
		{2, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end), 1},

		// Instance start color
		{3, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start_color), 1},

		// Instance end color
		{4, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end_color), 1}
	};

	bool success = setup_instanced_rendering<LineInstance>(shader_config,
														   uniforms,
														   setup_line_geometry,
														   instance_attributes,
														   lines_shader_,
														   lines_vao_,
														   lines_vbo_,
														   lines_instance_vbo_);

	if (success) {
		// Create separate VBO for medium line instances using the same VAO
		glBindVertexArray(lines_vao_);
		glGenBuffers(1, &wireframe_instance_vbo_);
		glBindVertexArray(0);
	}

	return success;
}

// Check if a point is inside the mesh geometry
bool Renderer::is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const {
	glm::dvec3 dpoint(point.x, point.y, point.z);
	return simulator.is_point_inside_geometry(dpoint);
}

// OpenGL state management helpers to minimize redundant state changes

void Renderer::use_shader_program(GLuint program_id) const {
	if (current_shader_program_ != program_id) {
		glUseProgram(program_id);
		current_shader_program_ = program_id;
	}
}

void Renderer::enable_blending() const {
	if (!blend_enabled_) {
		glEnable(GL_BLEND);
		blend_enabled_ = true;
	}
}

void Renderer::disable_blending() const {
	if (blend_enabled_) {
		glDisable(GL_BLEND);
		blend_enabled_ = false;
	}
}

/**
 * @brief Setup OpenGL resources for instanced triangle rendering
 * @return true if setup successful, false otherwise
 *
 * Creates VAO, VBO, and shader program for efficient batch triangle rendering.
 * Uses specialized triangle instancing shaders for medium geometry visualization.
 */
bool Renderer::setup_triangle_instanced_rendering() {
	// Load and compile triangle instance shaders
	std::string vertex_source = load_shader_source("shaders/triangles.vert");
	std::string fragment_source = load_shader_source("shaders/triangles.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load triangle instanced shader sources" << std::endl;
		return false;
	}

	triangles_shader_ = create_shader_program(vertex_source, fragment_source);
	if (!triangles_shader_) {
		std::cerr << "Failed to create triangle instanced shader program" << std::endl;
		return false;
	}

	// Create VAO for triangle instances
	glGenVertexArrays(1, &triangles_vao_);
	glBindVertexArray(triangles_vao_);

	// Create base triangle geometry (single triangle)
	glGenBuffers(1, &triangles_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo_);

	// Simple triangle vertices (will be transformed per instance)
	float triangle_vertices[] = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};

	glBufferData(GL_ARRAY_BUFFER, sizeof(triangle_vertices), triangle_vertices, GL_STATIC_DRAW);

	// Setup vertex positions (layout location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Create instance data buffer
	glGenBuffers(1, &triangles_instance_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_instance_vbo_);

	// Setup instance attributes (per-triangle data)
	// v0 (layout location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v0));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// v1 (layout location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v1));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// v2 (layout location 3)
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, v2));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// Color (layout location 4)
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	// Normal (layout location 5)
	glVertexAttribPointer(
		5, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleInstance), (void*)offsetof(TriangleInstance, normal));
	glEnableVertexAttribArray(5);
	glVertexAttribDivisor(5, 1);

	glBindVertexArray(0);
	return true;
}

/**
 * @brief Setup separate VAO for medium line rendering
 * @return true if setup successful, false otherwise
 *
 * Creates VAO and VBO for medium lines that shares the same shader program as regular lines
 * but has separate vertex buffer to avoid interference between photon paths and medium geometry.
 */
bool Renderer::setup_medium_line_vao() {
	// Create VAO for medium line instances
	glGenVertexArrays(1, &wireframe_vao_);
	glBindVertexArray(wireframe_vao_);

	// Bind the same base line geometry VBO as regular lines
	glBindBuffer(GL_ARRAY_BUFFER, lines_vbo_);

	// Setup vertex positions (layout location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Create separate instance data buffer for medium lines
	glGenBuffers(1, &wireframe_instance_vbo_);
	glBindBuffer(GL_ARRAY_BUFFER, wireframe_instance_vbo_);

	// Setup instance attributes (per-line data)
	// Start position (layout location 1)
	glVertexAttribPointer(
		1, 3, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, start));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// End position (layout location 2)
	glVertexAttribPointer(
		2, 3, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, end));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// Start color (layout location 3)
	glVertexAttribPointer(
		3, 4, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, start_color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// End color (layout location 4)
	glVertexAttribPointer(
		4, 4, GL_FLOAT, GL_FALSE, sizeof(MediumLineInstance), (void*)offsetof(MediumLineInstance, end_color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	glBindVertexArray(0);
	return true;
}

/**
 * @brief Add a triangle instance to the batch
 */
void Renderer::add_triangle_instance(const glm::vec3& v0,
									 const glm::vec3& v1,
									 const glm::vec3& v2,
									 const glm::vec4& color,
									 const glm::vec3& normal) {
	TriangleInstance instance;
	instance.v0 = v0;
	instance.v1 = v1;
	instance.v2 = v2;
	instance.color = color;
	instance.normal = normal;
	triangle_instances_.push_back(instance);
}

/**
 * @brief Render all triangle instances in a single batch draw call
 */
void Renderer::draw_triangle_instances() {
	if (triangle_instances_.empty() || !triangles_shader_) {
		return;
	}

	if (!triangle_instances_uploaded_) {
		glBindBuffer(GL_ARRAY_BUFFER, triangles_instance_vbo_);
		glBufferData(GL_ARRAY_BUFFER,
					 triangle_instances_.size() * sizeof(TriangleInstance),
					 triangle_instances_.data(),
					 GL_DYNAMIC_DRAW);
		triangle_instances_uploaded_ = true;
	}

	use_shader_program(triangles_shader_);

	glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
	GLint mvp_location = glGetUniformLocation(triangles_shader_, "mvp");
	if (mvp_location >= 0) {
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
	}

	glBindVertexArray(triangles_vao_);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 3, static_cast<GLsizei>(triangle_instances_.size()));
	glBindVertexArray(0);
}

/**
 * @brief Add a medium line instance to the batch
 */
void Renderer::add_medium_line_instance(const glm::vec3& start, const glm::vec3& end, const glm::vec4& color) {
	MediumLineInstance instance;
	instance.start = start;
	instance.end = end;
	instance.start_color = color;
	instance.end_color = color; // Same color for both ends (solid line)
	wireframe_instances_.push_back(instance);
}

/**
 * @brief Render all medium line instances in a single batch draw call
 */
void Renderer::draw_medium_line_instances() {
	if (wireframe_instances_.empty() || !lines_shader_) {
		return;
	}

	if (!medium_line_instances_uploaded_) {
		glBindBuffer(GL_ARRAY_BUFFER, wireframe_instance_vbo_);
		
		// Upload medium line instances to the separate medium line VBO
		glBufferData(GL_ARRAY_BUFFER,
					 wireframe_instances_.size() * sizeof(MediumLineInstance),
					 wireframe_instances_.data(),
					 GL_DYNAMIC_DRAW);
		medium_line_instances_uploaded_ = true;
	}

	use_shader_program(lines_shader_);

	glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
	glUniformMatrix4fv(line_instanced_mvp_uniform_location_, 1, GL_FALSE, &mvp[0][0]);

	// Use separate VAO for medium lines (no buffer rebinding needed)
	glBindVertexArray(wireframe_vao_);
	glDrawArraysInstanced(GL_LINES, 0, 2, static_cast<GLsizei>(wireframe_instances_.size()));
	glBindVertexArray(0);
}
