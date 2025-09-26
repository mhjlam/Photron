#include "renderer.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <format>
#include <iostream>
#include <sstream>
#include <set>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "math/range.hpp"
#include "renderer/camera.hpp"
#include "renderer/settings.hpp"
#include "simulator/config.hpp"
#include "simulator/config.hpp"
#include "simulator/layer.hpp"
#include "simulator/medium.hpp"
#include "simulator/simulator.hpp"

Renderer::Renderer() = default;

Renderer::~Renderer() {
	// Clean up background sorting thread if active
	if (background_sort_in_progress_.load() && sorting_future_.valid()) {
		sorting_future_.wait(); // Wait for background sort to complete
	}
	
	// Cleanup OpenGL resources
	if (line_vao_) glDeleteVertexArrays(1, &line_vao_);
	if (line_vbo_) glDeleteBuffers(1, &line_vbo_);
	if (point_vao_) glDeleteVertexArrays(1, &point_vao_);
	if (point_vbo_) glDeleteBuffers(1, &point_vbo_);
	if (triangle_vao_) glDeleteVertexArrays(1, &triangle_vao_);
	if (triangle_vbo_) glDeleteBuffers(1, &triangle_vbo_);
	
	// Cleanup instanced rendering resources
	if (voxel_cube_vao_) glDeleteVertexArrays(1, &voxel_cube_vao_);
	if (voxel_cube_vbo_) glDeleteBuffers(1, &voxel_cube_vbo_);
	if (voxel_instance_vbo_) glDeleteBuffers(1, &voxel_instance_vbo_);
	if (line_instanced_vao_) glDeleteVertexArrays(1, &line_instanced_vao_);
	if (line_instanced_vbo_) glDeleteBuffers(1, &line_instanced_vbo_);
	if (line_instance_vbo_) glDeleteBuffers(1, &line_instance_vbo_);
	
	// Cleanup shader programs
	if (line_shader_program_) glDeleteProgram(line_shader_program_);
	if (point_shader_program_) glDeleteProgram(point_shader_program_);
	if (triangle_shader_program_) glDeleteProgram(triangle_shader_program_);
	if (voxel_shader_program_) glDeleteProgram(voxel_shader_program_);
	if (line_instanced_shader_program_) glDeleteProgram(line_instanced_shader_program_);
}

bool Renderer::initialize() {
	if (Config::is_initialized() && Config::get().log()) {
		std::cout << "Initializing modern OpenGL 4.5 renderer..." << std::endl;
	}

	setup_opengl();

	// Initialize consolidated rendering systems
	if (!setup_line_rendering()) {
		std::cerr << "Failed to setup consolidated line rendering" << std::endl;
		return false;
	}

	if (!setup_point_rendering()) {
		std::cerr << "Failed to setup consolidated point rendering" << std::endl;
		return false;
	}

	if (!setup_triangle_rendering()) {
		std::cerr << "Failed to setup consolidated triangle rendering" << std::endl;
		return false;
	}

	if (!setup_voxel_instanced_rendering()) {
		std::cerr << "Failed to setup instanced voxel rendering" << std::endl;
		return false;
	}

	if (!setup_line_instanced_rendering()) {
		std::cerr << "Failed to setup instanced line rendering" << std::endl;
		return false;
	}

	if (!setup_point_instanced_rendering()) {
		std::cerr << "Failed to setup instanced point rendering" << std::endl;
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
	if (!is_arc_camera_mode_) {
		const float move_speed = 0.02f; // Reduced speed for smoother movement

		if (key_state_.w_pressed)
			camera_.move_forward(move_speed);
		if (key_state_.s_pressed)
			camera_.move_backward(move_speed);
		if (key_state_.a_pressed)
			camera_.move_left(move_speed);
		if (key_state_.d_pressed)
			camera_.move_right(move_speed);
		if (key_state_.q_pressed)
			camera_.move_down(move_speed);
		if (key_state_.e_pressed)
			camera_.move_up(move_speed);
	}
}

void Renderer::render(Simulator& simulator) {
	// Store simulator reference for use in drawing functions
	simulator_ = &simulator;
	
	// PERFORMANCE FIX: Only invalidate caches when simulation data actually changes
	uint64_t current_sim_version = simulator.get_simulation_version();
	if (current_sim_version != last_simulation_version_) {
		// When simulation data changes, invalidate all relevant caches
		// This includes both energy caches AND path instances since the path data may have changed
		invalidate_all_caches();
		
		last_simulation_version_ = current_sim_version;
	}

	// Update camera target based on simulation bounds
	static bool first_frame = true;
	if (first_frame) {
		update_camera_target(simulator);
		first_frame = false;
	}

	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Update camera and matrices
	update_camera();

	// Update energy label screen positions once per frame (fixes label positioning stability)
	update_energy_label_screen_positions();

	// Draw everything based on current settings (user request 4 - re-enable ImGui menu)
	if (settings_.draw_volume) {
		draw_volume(simulator);  // Geometry-aware bounds
	}

	// Clear depth buffer after drawing wireframes so voxels and lines can have proper depth ordering
	glClear(GL_DEPTH_BUFFER_BIT);

	// Render voxels using optimized instanced rendering
	draw_voxels_instanced(settings_);

	// Draw paths after voxels with proper depth testing
	if (settings_.draw_paths) {
		draw_paths_instanced(settings_); // High-performance instanced version
	}

	// Draw energy labels as billboards (user request 3)
	draw_labels(settings_);
}

void Renderer::draw_coordinate_axes() {
	// Draw coordinate axes with very muted colors to not interfere with physics visualization
	begin_lines();

	// X axis - Very muted red
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.5f, 0.0f, 0.0f), glm::vec4(0.4f, 0.1f, 0.1f, 0.3f));

	// Y axis - Very muted green (shortened to avoid confusion with physics rays)
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.5f, 0.0f), glm::vec4(0.1f, 0.4f, 0.1f, 0.3f));

	// Z axis - Very muted blue
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.5f), glm::vec4(0.1f, 0.1f, 0.4f, 0.3f));

	end_lines();
	draw_lines();
}

void Renderer::draw_test_geometry() {
	// Draw a bright magenta triangle for testing
	begin_lines();

	glm::vec4 magenta(1.0f, 0.0f, 1.0f, 1.0f);
	add_line(glm::vec3(-0.5f, -0.5f, 0.0f), glm::vec3(0.5f, -0.5f, 0.0f), magenta);
	add_line(glm::vec3(0.5f, -0.5f, 0.0f), glm::vec3(0.0f, 0.5f, 0.0f), magenta);
	add_line(glm::vec3(0.0f, 0.5f, 0.0f), glm::vec3(-0.5f, -0.5f, 0.0f), magenta);

	end_lines();
	draw_lines();

	// Draw some test points
	begin_points();
	add_point(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec4(1.0f, 1.0f, 0.0f, 1.0f));
	add_point(glm::vec3(-1.0f, 1.0f, -1.0f), glm::vec4(0.0f, 1.0f, 1.0f, 1.0f));
	end_points();
	draw_points();
}

void Renderer::draw_volume(const Simulator& simulator) {
	// Combine wireframe (lines) and faces (triangles) for a complete boundary visualization
	begin_lines();
	begin_triangles();

	// Face culling is already disabled by default, but ensure blending is set up properly
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glm::vec4 wireframe_color(0.7f, 0.7f, 0.7f, 0.8f);  // Bright wireframe
	glm::vec4 face_color(0.2f, 0.2f, 0.2f, 0.1f);       // Subtle transparent faces

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
			std::string plane_key = std::format("{},{},{},{}",
				int(normal.x * 1000),
				int(normal.y * 1000), 
				int(normal.z * 1000),
				int(d * 1000));
			
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
				add_triangle(vertices[0], vertices[1], vertices[2], face_color);
				
				// Add wireframe edges
				add_line(vertices[0], vertices[1], wireframe_color);
				add_line(vertices[1], vertices[2], wireframe_color);
				add_line(vertices[2], vertices[0], wireframe_color);
			}
			else if (vertices.size() >= 6 && vertices.size() % 3 == 0) {
				// Multiple triangles forming a face - check if they form a quad
				
				// For now, let's find unique vertices and try to form a quadrilateral
				std::vector<glm::vec3> unique_vertices;
				const float epsilon = 0.001f;
				
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
					glm::vec3 ref_vec = abs(normal.x) < 0.9f ? glm::vec3(1.0f, 0.0f, 0.0f) : glm::vec3(0.0f, 1.0f, 0.0f);
					glm::vec3 tangent = normalize(cross(normal, ref_vec));
					glm::vec3 bitangent = normalize(cross(normal, tangent));
					
					// Sort vertices by angle around the center - Modern C++20
					std::ranges::sort(unique_vertices, 
						[&center, &tangent, &bitangent](const glm::vec3& a, const glm::vec3& b) noexcept {
							const glm::vec3 dir_a = a - center;
							const glm::vec3 dir_b = b - center;
							
							const float angle_a = atan2f(dot(dir_a, bitangent), dot(dir_a, tangent));
							const float angle_b = atan2f(dot(dir_b, bitangent), dot(dir_b, tangent));
							
							return angle_a < angle_b;
						});
					
					// Render as two triangles forming a quad
					add_triangle(unique_vertices[0], unique_vertices[1], unique_vertices[2], face_color);
					add_triangle(unique_vertices[0], unique_vertices[2], unique_vertices[3], face_color);
					
					// Add wireframe edges for the quad
					add_line(unique_vertices[0], unique_vertices[1], wireframe_color);
					add_line(unique_vertices[1], unique_vertices[2], wireframe_color);
					add_line(unique_vertices[2], unique_vertices[3], wireframe_color);
					add_line(unique_vertices[3], unique_vertices[0], wireframe_color);
				}
				else {
					// Fallback: render all triangles individually
					for (size_t i = 0; i < vertices.size(); i += 3) {
						add_triangle(vertices[i], vertices[i+1], vertices[i+2], face_color);
						
						// Add wireframe edges
						add_line(vertices[i], vertices[i+1], wireframe_color);
						add_line(vertices[i+1], vertices[i+2], wireframe_color);
						add_line(vertices[i+2], vertices[i], wireframe_color);
					}
				}
			}
		}
	}

	end_lines();
	end_triangles();
	draw_lines();
	draw_triangles();

	// No need to change OpenGL state - keep defaults (face culling disabled, blending enabled)
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
	if (is_arc_camera_mode_ && action == GLFW_PRESS) {
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
	if (key == GLFW_KEY_W)
		key_state_.w_pressed = (action != GLFW_RELEASE);
	if (key == GLFW_KEY_A)
		key_state_.a_pressed = (action != GLFW_RELEASE);
	if (key == GLFW_KEY_S)
		key_state_.s_pressed = (action != GLFW_RELEASE);
	if (key == GLFW_KEY_D)
		key_state_.d_pressed = (action != GLFW_RELEASE);
	if (key == GLFW_KEY_Q)
		key_state_.q_pressed = (action != GLFW_RELEASE);
	if (key == GLFW_KEY_E)
		key_state_.e_pressed = (action != GLFW_RELEASE);
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
	is_arc_camera_mode_ = is_arc_mode;
	camera_.set_fps_mode(!is_arc_mode); // FPS mode when not Orbit mode
}

bool Renderer::should_capture_mouse() const {
	return camera_.should_capture_mouse();
}

void Renderer::draw_voxels(const Settings& settings) {
	if (!simulator_) {
		return;
	}

	// Only draw voxels if voxel rendering is enabled
	if (!settings.draw_voxels) {
		return;
	}

	// Simple and robust approach: use the simulator's voxelization data
	// No complex clipping planes needed - just trust the inside/outside classification

	// Setup OpenGL state for transparent voxel rendering
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE); // Disable depth writing for transparency
	glDisable(GL_CULL_FACE);

	// Enable global clipping planes (disabled - using per-voxel clipping instead)
	// for (size_t i = 0; i < std::min(clipping_planes.size(), size_t(6)); i++) {
	//     glEnable(GL_CLIP_DISTANCE0 + i);
	// }

	begin_triangles();

	// Get MCML grid parameters - EXACTLY like backup
	const auto& config = Config::get();
	int nx = config.nx();
	int ny = config.ny();
	int nz = config.nz();
	double voxsize = config.vox_size();
	double half_voxsize = voxsize * 0.5;

	// Collect voxels with distance for depth sorting (user request 4)
	struct VoxelRenderData
	{
		Voxel* voxel;
		glm::vec3 position;
		float distance_to_camera;
		glm::vec4 color;
	};

	std::vector<VoxelRenderData> voxels_to_render;

	// Get camera position for distance calculation
	glm::vec3 camera_pos = camera_.get_position();

	// DEBUG: Print when starting voxel rendering
	std::cout << "VOXEL RENDERING START: mode=" << static_cast<int>(settings.voxel_mode) 
	          << " (0=Absorption, 1=Emittance, 2=Layers, 3=Combined)" << std::endl;

	// First pass: analyze energy distribution to find dynamic range for adaptive scaling
	std::vector<float> all_energies;
	float total_accumulated_energy = 0.0f;

	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// Use Volume coordinate-based access instead of linear indexing
				if (ix >= 0 && iy >= 0 && iz >= 0 && 
					static_cast<uint32_t>(ix) < static_cast<uint32_t>(nx) && 
					static_cast<uint32_t>(iy) < static_cast<uint32_t>(ny) && 
					static_cast<uint32_t>(iz) < static_cast<uint32_t>(nz)) {
					
					Voxel* voxel = simulator_->voxel_grid(static_cast<uint32_t>(ix), static_cast<uint32_t>(iy), static_cast<uint32_t>(iz));

					if (voxel) {
						float absorption = static_cast<float>(voxel->absorption);
						float emittance = static_cast<float>(voxel->total_emittance()); // Combines all emittance types including specular reflection

						// DEBUG: Print ALL voxels that have emittance values
						if (emittance > 1e-10f) {
							std::cout << "RENDERER EMITTANCE: grid(" << ix << "," << iy << "," << iz 
									  << ") voxel(" << voxel->ix() << "," << voxel->iy() << "," << voxel->iz()
									  << ") emittance=" << emittance << " absorption=" << absorption << std::endl;
						}

						// DEBUG: Show when we're processing known emittance coordinates  
						if ((voxel->ix() == 12 && voxel->iy() == 2 && voxel->iz() == 12) || 
						    (voxel->ix() == 12 && voxel->iy() == 4 && voxel->iz() == 11) ||
						    (voxel->ix() == 12 && voxel->iy() == 3 && voxel->iz() == 12) ||
						    (voxel->ix() == 11 && voxel->iy() == 4 && voxel->iz() == 12) ||
						    (voxel->ix() == 12 && voxel->iy() == 2 && voxel->iz() == 13) ||
						    (voxel->ix() == 12 && voxel->iy() == 3 && voxel->iz() == 13)) {
							std::cout << "KNOWN EMITTANCE VOXEL: grid(" << ix << "," << iy << "," << iz 
									  << ") voxel(" << voxel->ix() << "," << voxel->iy() << "," << voxel->iz()
									  << ") emittance=" << emittance << " mode=" << static_cast<int>(settings.voxel_mode) << std::endl;
						}

						float total_energy;
						if (settings.voxel_mode == VoxelMode::Absorption) {
							total_energy = absorption;
						}
						else if (settings.voxel_mode == VoxelMode::Emittance) {
							total_energy = emittance;
						}
						else if (settings.voxel_mode == VoxelMode::Layers) {
							// For Layers mode, use a constant value for all voxels with material
							total_energy = (voxel->material != nullptr) ? 1.0f : 0.0f;
						}
						else {
							total_energy = absorption + emittance;
						}

						if (total_energy > 0.0000001f) {
							all_energies.push_back(total_energy);
							total_accumulated_energy += total_energy;
						}
					}
				}
			}
		}
	}

	// Calculate adaptive energy range for better contrast
	float min_energy = 1e-8f; // Much lower threshold for emittance detection
	float max_energy = 0.01f; // Default fallback

	if (!all_energies.empty()) {
		// Modern C++20: Use ranges::sort for cleaner syntax
		std::ranges::sort(all_energies);

		// Use percentile-based scaling for better contrast
		const size_t count = all_energies.size();
		const float p5 = all_energies[static_cast<size_t>(count * 0.05)];  // 5th percentile
		const float p95 = all_energies[static_cast<size_t>(count * 0.95)]; // 95th percentile

		min_energy = std::max(p5, 1e-8f);                            // Lower minimum
		max_energy = std::max(p95, min_energy * 10.0f);              // Ensure reasonable range
	}

	// Second pass: collect all voxels with their render data using adaptive scaling
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// Use Volume coordinate-based access
				if (ix >= 0 && iy >= 0 && iz >= 0 && 
					static_cast<uint32_t>(ix) < static_cast<uint32_t>(nx) && 
					static_cast<uint32_t>(iy) < static_cast<uint32_t>(ny) && 
					static_cast<uint32_t>(iz) < static_cast<uint32_t>(nz)) {
					
					Voxel* voxel = simulator_->voxel_grid(static_cast<uint32_t>(ix), static_cast<uint32_t>(iy), static_cast<uint32_t>(iz));

					if (voxel) {
						// Calculate center position using the exact same bounds used for voxelization
						auto bounds = simulator_->mediums[0].get_bounds();
						double x = bounds.min_bounds.x + (voxsize * ix) + half_voxsize;
						double y = bounds.min_bounds.y + (voxsize * iy) + half_voxsize;
						double z = bounds.min_bounds.z + (voxsize * iz) + half_voxsize;

						// Check if this voxel is actually inside the mesh geometry AND has material
						bool is_in_geometry = simulator_->is_point_inside_geometry({x,y,z});
						bool has_tissue = (voxel->material != nullptr);
						
						// For complex geometries like pyramids, we need to check if any part of the voxel
						// intersects with the geometry, not just the center point
						bool voxel_intersects_geometry = false;
						
						if (has_tissue) {
							// Check voxel center first (fastest check)
							if (is_in_geometry) {
								voxel_intersects_geometry = true;
							} else {
								// If center is outside, check the 8 corners to see if any are inside
								double half_voxel = voxsize * 0.5;
								std::vector<glm::dvec3> corners = {
									{x - half_voxel, y - half_voxel, z - half_voxel},
									{x + half_voxel, y - half_voxel, z - half_voxel},
									{x - half_voxel, y + half_voxel, z - half_voxel},
									{x + half_voxel, y + half_voxel, z - half_voxel},
									{x - half_voxel, y - half_voxel, z + half_voxel},
									{x + half_voxel, y - half_voxel, z + half_voxel},
									{x - half_voxel, y + half_voxel, z + half_voxel},
									{x + half_voxel, y + half_voxel, z + half_voxel}
								};
								
								// If any corner is inside the geometry, this voxel intersects
								for (const auto& corner : corners) {
									if (simulator_->is_point_inside_geometry(corner)) {
										voxel_intersects_geometry = true;
										break;
									}
								}
							}
						}
						
						// Only render voxels that have material AND intersect with the geometry
						if (!has_tissue || !voxel_intersects_geometry) {
							continue;
						}

						float fx = static_cast<float>(x);
						float fy = static_cast<float>(y);
						float fz = static_cast<float>(z);
						glm::vec3 voxel_pos(fx, fy, fz);

						// Use real emittance data from simulation - including specular reflection at entry
						float absorption = static_cast<float>(voxel->absorption);
						float emittance = static_cast<float>(voxel->total_emittance()); // Combines emittance_transmitted + emittance_reflected + emittance_diffuse

						float total_energy;
						if (settings.voxel_mode == VoxelMode::Absorption) {
							total_energy = absorption;
						}
						else if (settings.voxel_mode == VoxelMode::Emittance) {
							total_energy = emittance;
						}
						else if (settings.voxel_mode == VoxelMode::Layers) {
							// For Layers mode, we ignore energy and just show material types
							total_energy = 1.0f;                   // Constant energy for all voxels with material
						}
						else {
							total_energy = absorption + emittance; // Fallback to combined
						}

						// Show even the tiniest energy interactions
						glm::vec4 color(0.0f, 0.0f, 0.0f, 0.0f); // Default transparent

						// Check if this voxel is actually inside the mesh geometry
						glm::dvec3 voxel_center = glm::dvec3(x, y, z);
						bool is_inside_geometry = simulator_->is_point_inside_geometry(voxel_center);

						// In Layers mode, show all voxels with material regardless of energy
						bool should_render_voxel = false;
						if (settings.voxel_mode == VoxelMode::Layers) {
							should_render_voxel = (voxel->material != nullptr) && is_inside_geometry;
						}
						else {
							should_render_voxel = (total_energy > 1e-8f || voxel->material != nullptr) && is_inside_geometry;
						}

						if (should_render_voxel) {
							// For Layers mode, use simplified rendering logic
							if (settings.voxel_mode == VoxelMode::Layers) {
								if (voxel->material != nullptr) {
									// Match the appearance of Absorption mode when there's no absorption
									// Use the same energy and alpha values as "no energy" voxels in absorption mode
									color = get_layer_specific_energy_color(1e-8f, min_energy, max_energy,
																			voxel->material->id());
									color.a = 0.05f; // Same very faint alpha as no-absorption voxels
								}
							}
							else {
								// Original energy-based rendering logic for other modes

								// Use gamma correction for better mid-range contrast
								float normalized_energy = (total_energy - min_energy) / (max_energy - min_energy);
								normalized_energy = std::clamp(normalized_energy, 0.0f, 1.0f);

								// Apply gamma correction for better contrast (gamma = 0.5 brightens mid-tones)
								float gamma_corrected = std::pow(normalized_energy, 0.5f);

								// Calculate distance from origin for depth-based alpha
								float max_dist = glm::length(to_float(bounds.max_bounds));
								float dist_from_origin = glm::length(voxel_pos);
								float normalized_dist = glm::clamp(dist_from_origin / max_dist, 0.0f, 1.0f);

								// More visible alpha scaling based on gamma-corrected energy
								float min_alpha = 0.2f; // Higher minimum for emittance visibility
								float max_alpha = 0.9f; // Higher maximum visibility

								// Special case for emittance mode: boost visibility significantly
								if (settings.voxel_mode == VoxelMode::Emittance) {
									min_alpha = 0.6f; // Much more visible
									max_alpha = 1.0f; // Full opacity for high emittance
								}

								float alpha =
									min_alpha
									+ (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.2f * normalized_dist);

								// Use layer-specific energy color mapping with dynamic range
								if (total_energy > 1e-8f) { // Lower threshold
									color = get_layer_specific_energy_color(total_energy, min_energy, max_energy,
																			voxel->material->id());
									color.a = alpha;
								}
								else if (voxel->material != nullptr) {
									// Voxel in medium but no recorded energy: very faint material-colored hint
									color = get_layer_specific_energy_color(1e-8f, min_energy, max_energy,
																			voxel->material->id());
									color.a = 0.05f; // Slightly more visible
								}
							}

							// Calculate distance to camera for depth sorting
							float distance = glm::length(voxel_pos - camera_pos);

							// Add to render list (only if it has visible color)
							voxels_to_render.push_back({voxel, voxel_pos, distance, color});
							
							// DEBUG: Print when voxels with emittance are added to render list
							if (settings.voxel_mode == VoxelMode::Emittance && emittance > 1e-8f) {
								std::cout << "RENDER EMITTANCE: grid(" << ix << "," << iy << "," << iz 
										  << ") voxel(" << voxel->ix() << "," << voxel->iy() << "," << voxel->iz()
										  << ") emittance=" << emittance << " pos=(" << fx << "," << fy << "," << fz << ")"
										  << " color=(" << color.r << "," << color.g << "," << color.b << "," << color.a << ")" << std::endl;
							}
						}
					}
				}
			}
		}
	}

	// PERFORMANCE OPTIMIZATION: Skip expensive sorting for large datasets  
	// Sorting 981k+ voxels every frame is prohibitively expensive (O(n log n) = ~19M operations)
	// Only sort when voxel count is reasonable for perfect transparency
	if (voxels_to_render.size() < 50000) {
		// Sort voxels back-to-front for optimal transparency blending
		std::ranges::sort(voxels_to_render, [](const VoxelRenderData& a, const VoxelRenderData& b) noexcept {
			return a.distance_to_camera > b.distance_to_camera; // Furthest first
		});
	}
	// For large datasets: rely on depth buffer + alpha blending without sorting
	// For large datasets: rely on depth buffer + alpha blending without sorting

	// Now render all voxels in batches with depth offsets to prevent Z-fighting
	const size_t BATCH_SIZE = 300; // Render 300 voxels at a time (safe limit)
	size_t total_voxels = voxels_to_render.size();

	for (size_t batch_start = 0; batch_start < total_voxels; batch_start += BATCH_SIZE) {
		begin_triangles();         // Start new batch
		size_t batch_end = std::min(batch_start + BATCH_SIZE, total_voxels);

		// Add triangles for current batch
		for (size_t i = batch_start; i < batch_end; i++) {
			const auto& voxel_data = voxels_to_render[i];
			float half = static_cast<float>(half_voxsize * 0.95);
			glm::vec3 pos = voxel_data.position;

			// Add small depth offset based on voxel index to prevent Z-fighting
			float depth_offset = static_cast<float>(i) * 0.000001f;
			pos.z += depth_offset;

			// Create voxel corners
			glm::vec3 corners[8] = {
				glm::vec3(pos.x - half, pos.y - half, pos.z - half), // 0
				glm::vec3(pos.x + half, pos.y - half, pos.z - half), // 1
				glm::vec3(pos.x + half, pos.y + half, pos.z - half), // 2
				glm::vec3(pos.x - half, pos.y + half, pos.z - half), // 3
				glm::vec3(pos.x - half, pos.y - half, pos.z + half), // 4
				glm::vec3(pos.x + half, pos.y - half, pos.z + half), // 5
				glm::vec3(pos.x + half, pos.y + half, pos.z + half), // 6
				glm::vec3(pos.x - half, pos.y + half, pos.z + half)  // 7
			};

			// Render complete voxel - trust the simulator's voxelization
			// The simulator already determined this voxel should be rendered
			add_triangle(corners[4], corners[5], corners[6], voxel_data.color);
			add_triangle(corners[4], corners[6], corners[7], voxel_data.color);
			add_triangle(corners[1], corners[0], corners[3], voxel_data.color);
			add_triangle(corners[1], corners[3], corners[2], voxel_data.color);
			add_triangle(corners[0], corners[4], corners[7], voxel_data.color);
			add_triangle(corners[0], corners[7], corners[3], voxel_data.color);
			add_triangle(corners[1], corners[2], corners[6], voxel_data.color);
			add_triangle(corners[1], corners[6], corners[5], voxel_data.color);
			add_triangle(corners[0], corners[1], corners[5], voxel_data.color);
			add_triangle(corners[0], corners[5], corners[4], voxel_data.color);
			add_triangle(corners[3], corners[7], corners[6], voxel_data.color);
			add_triangle(corners[3], corners[6], corners[2], voxel_data.color);
		}
		end_triangles(); // End the current batch

		// Render current batch normally (no clipping planes needed)
		draw_triangles();
	}

	// Restore OpenGL state to defaults
	glDepthMask(GL_TRUE);      // Re-enable depth writing
	glEnable(GL_DEPTH_TEST);   // Re-enable depth testing
	glDisable(GL_CULL_FACE);   // Keep face culling disabled (default state)
	
	// Disable all clipping planes
	for (int i = 0; i < 6; i++) {
		glDisable(GL_CLIP_DISTANCE0 + i);
	}
}

// HIGH-PERFORMANCE instanced voxel rendering - replaces the slow triangle-based approach
void Renderer::draw_voxels_instanced(const Settings& settings) {
	if (!simulator_) {
		return;
	}

	// Only draw voxels if voxel rendering is enabled
	if (!settings.draw_voxels) {
		return;
	}

	// PERFORMANCE OPTIMIZATION: Check if we need to rebuild voxel instances
	// Only rebuild when simulation data changes or voxel mode changes
	static VoxelMode last_voxel_mode = VoxelMode::Layers;
	
	bool need_rebuild = voxel_instances_dirty_ || 
	                    settings.voxel_mode != last_voxel_mode;
						
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
	int nx = config.nx();
	int ny = config.ny();
	int nz = config.nz();
	double voxsize = config.vox_size();
	double half_voxsize = voxsize * 0.5;
	float voxel_scale = static_cast<float>(voxsize * 0.95); // Slightly smaller to show gaps

	// PERFORMANCE OPTIMIZATION: Use cached energy range instead of expensive per-frame analysis
	update_cached_energy_range(settings);

	// Set up instanced rendering for maximum performance
	begin_voxel_instances();
	
	// PERFORMANCE: Pre-calculate camera data and bounds once for all voxels
	glm::vec3 camera_pos = camera_.get_position();
	glm::vec3 camera_target = camera_.get_target();
	glm::vec3 camera_front = glm::normalize(camera_target - camera_pos);
	float max_render_distance = 50.0f; // Don't render voxels too far away
	
	// PERFORMANCE: Cache bounds calculation - don't call it 981k times!
	auto bounds = simulator_->get_combined_bounds();
	
	// Render voxels using cached energy scaling - PERFORMANCE OPTIMIZED
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// Bounds check
				if (ix >= 0 && iy >= 0 && iz >= 0 && 
					static_cast<uint32_t>(ix) < static_cast<uint32_t>(nx) && 
					static_cast<uint32_t>(iy) < static_cast<uint32_t>(ny) && 
					static_cast<uint32_t>(iz) < static_cast<uint32_t>(nz)) {
					
					Voxel* voxel = simulator_->voxel_grid(static_cast<uint32_t>(ix), static_cast<uint32_t>(iy), static_cast<uint32_t>(iz));
					if (!voxel || !voxel->material) continue; // Skip voxels with no material

					// Calculate energy based on current mode
					float absorption = static_cast<float>(voxel->absorption);
					float emittance = static_cast<float>(voxel->total_emittance()); // Includes specular reflection from entry
					
					float total_energy;
					
					if (settings.voxel_mode == VoxelMode::Absorption) {
						total_energy = absorption;
					}
					else if (settings.voxel_mode == VoxelMode::Emittance) {
						total_energy = emittance;
					}
					else if (settings.voxel_mode == VoxelMode::Layers) {
						total_energy = 1.0f; // We already know voxel->material != nullptr
					}
					else {
						total_energy = absorption + emittance;
					}

					// IMPORTANT: Don't skip voxels with no energy! 
					// They should still render as very faint material-colored hints

					// Calculate world position using cached bounds (PERFORMANCE FIX)
					double x = bounds.min_bounds.x + (voxsize * ix) + half_voxsize;
					double y = bounds.min_bounds.y + (voxsize * iy) + half_voxsize;
					double z = bounds.min_bounds.z + (voxsize * iz) + half_voxsize;



					glm::vec3 voxel_pos(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));

					// PERFORMANCE: Distance-based culling to reduce rendered voxels
					float distance_to_camera = glm::length(voxel_pos - camera_pos);
					if (distance_to_camera > max_render_distance) {
						continue; // Skip voxels too far from camera
					}

					// PERFORMANCE OPTIMIZATION: Skip expensive geometry tests
					// The voxel already has material assigned, so it's already been determined to be in geometry
					// during the voxelization process. No need to re-test every frame.

					// EXACT ORIGINAL COLOR CALCULATION - DO NOT CHANGE!
					glm::vec4 color(0.0f);

					if (settings.voxel_mode == VoxelMode::Layers) {
						// Layers mode: show all material voxels regardless of energy
						if (voxel->material != nullptr) {
							// Match the appearance of Absorption mode when there's no absorption
							// Use the same energy and alpha values as "no energy" voxels in absorption mode
							color = get_layer_specific_energy_color(1e-8f, cached_min_energy_, cached_max_energy_, voxel->material->id());
							color.a = 0.05f; // Same very faint alpha as no-absorption voxels
						}
					} else {
						// Original energy-based rendering logic for other modes

						// Use gamma correction for better mid-range contrast
						float normalized_energy = (total_energy - cached_min_energy_) / (cached_max_energy_ - cached_min_energy_);
						normalized_energy = std::clamp(normalized_energy, 0.0f, 1.0f);

						// Apply gamma correction for better contrast (gamma = 0.5 brightens mid-tones)
						float gamma_corrected = std::pow(normalized_energy, 0.5f);

						// Calculate distance from origin for depth-based alpha (cached calculation)
						static float cached_max_dist = 0.0f;
						static bool max_dist_cached = false;
						if (!max_dist_cached) {
							cached_max_dist = glm::length(to_float(bounds.max_bounds));
							max_dist_cached = true;
						}
						
						float dist_from_origin = glm::length(voxel_pos);
						float normalized_dist = glm::clamp(dist_from_origin / cached_max_dist, 0.0f, 1.0f);

						// More visible alpha scaling based on gamma-corrected energy
						float min_alpha = 0.2f; // Higher minimum for emittance visibility
						float max_alpha = 0.9f; // Higher maximum visibility

						// EMITTANCE MODE FIX: Use more reasonable alpha values
						// The issue was min_alpha = 0.6f made all emittance voxels too opaque
						if (settings.voxel_mode == VoxelMode::Emittance) {
							// Use energy-proportional alpha for better emittance visualization
							min_alpha = 0.1f;  // Lower minimum - let the energy determine visibility
							max_alpha = 0.7f;  // Reasonable maximum for emittance
						}

						float alpha = min_alpha + (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.2f * normalized_dist);

						// Use layer-specific energy color mapping with dynamic range
						if (total_energy > 1e-8f) { // Lower threshold
							color = get_layer_specific_energy_color(total_energy, cached_min_energy_, cached_max_energy_, voxel->material->id());
							color.a = alpha;
						}
						else if (voxel->material != nullptr) {
							// Voxel in medium but no recorded energy: very faint material-colored hint
							color = get_layer_specific_energy_color(1e-8f, cached_min_energy_, cached_max_energy_, voxel->material->id());
							color.a = 0.05f; // Slightly more visible
						}
					}
					// Only add voxels that should be visible
					if (color.a > 0.0f) {
						// Calculate depth for sorting (distance from camera) - use pre-calculated camera data
						glm::vec3 view_dir = voxel_pos - camera_pos;
						float depth = glm::dot(view_dir, camera_front);
						
						add_voxel_instance(voxel_pos, color, voxel_scale, depth);
					}

				}
			}
		}
	}

	// Render ALL voxels in a single high-performance instanced draw call
	end_voxel_instances(settings.voxel_mode);
	draw_voxel_instances();

	// Restore OpenGL state
	glDepthMask(GL_TRUE);      // Re-enable depth writing
	glEnable(GL_DEPTH_TEST);   // Keep depth testing enabled
	glDisable(GL_BLEND);       // Disable blending
}

void Renderer::draw_paths(const Settings& /*settings*/) {
	if (!simulator_)
		return;

	// First pass: analyze energy distribution for adaptive logarithmic mapping
	std::vector<float> all_energies;
	for (const Photon& photon : simulator_->get_paths()) {
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

	// Calculate adaptive logarithmic mapping parameters
	float min_energy = 1.0f, max_energy = 0.0f;
	if (!all_energies.empty()) {
		auto [min_it, max_it] = std::minmax_element(all_energies.begin(), all_energies.end());
		min_energy = *min_it;
		max_energy = *max_it;
	}

	// Create adaptive logarithmic mapping function using helper method
	auto adaptive_log_color = [this, min_energy, max_energy](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, min_energy, max_energy);
	};

	// Draw photon path histories with adaptive energy-based coloring
	begin_lines();

	for (const Photon& photon : simulator_->get_paths()) {
		if (photon.path_head) {
			// Draw connected line segments with energy gradient exactly like backup
			auto current = photon.path_head;
			auto next = current ? current->next : nullptr;

			// Draw physically correct incident ray from source to material surface
			if (current && !simulator_->sources.empty()) {
				glm::vec3 first_interaction(static_cast<float>(current->position.x),
											static_cast<float>(current->position.y),
											static_cast<float>(current->position.z));

				// Use actual source parameters from the configuration
				const Source& source = simulator_->sources[0]; // Use first source
				glm::vec3 source_pos(static_cast<float>(source.origin.x), static_cast<float>(source.origin.y),
									 static_cast<float>(source.origin.z));
				glm::vec3 source_dir(static_cast<float>(source.direction.x), static_cast<float>(source.direction.y),
									 static_cast<float>(source.direction.z));

				// Calculate surface entry point (find topmost Y coordinate of any layer)
				float surface_y = 0.1f; // Default top surface
				const auto& layers = simulator_->get_all_layers();
				if (!layers.empty()) {
					// Find the highest Y coordinate among all layers by examining triangle vertices
					surface_y = -1000.0f; // Start very low
					for (const auto& layer : layers) {
						for (const auto& triangle : layer.mesh) {
							float y0 = static_cast<float>(triangle.v0().y);
							float y1 = static_cast<float>(triangle.v1().y);
							float y2 = static_cast<float>(triangle.v2().y);
							if (y0 > surface_y)
								surface_y = y0;
							if (y1 > surface_y)
								surface_y = y1;
							if (y2 > surface_y)
								surface_y = y2;
						}
					}
				}

				// Calculate where ray hits the surface
				if (source_dir.y != 0.0f) {
					float t = (surface_y - source_pos.y) / source_dir.y;
					glm::vec3 surface_entry = source_pos + t * source_dir;

					// Draw incident ray from source to surface
					glm::vec4 incident_color(1.0f, 1.0f, 1.0f, 1.0f); // Bright white
					add_line(source_pos, surface_entry, incident_color);

					// If photon actually enters material, draw refracted ray from surface to first interaction
					if (first_interaction.y < surface_y - 0.001f) {
						glm::vec4 refracted_color(0.9f, 0.9f, 1.0f, 0.8f); // Light blue
						add_line(surface_entry, first_interaction, refracted_color);
					}
				}
			}

			while (current && next) {
				// Use adaptive logarithmic coloring based on actual energy distribution
				float energy1 = static_cast<float>(current->value);
				float energy2 = static_cast<float>(next->value);

				glm::vec4 start_color = adaptive_log_color(energy1);
				glm::vec4 end_color = adaptive_log_color(energy2);
				start_color.a = 1.0f;
				end_color.a = 1.0f;

				glm::vec3 start(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								static_cast<float>(current->position.z));
				glm::vec3 end(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
							  static_cast<float>(next->position.z));

				// Draw actual photon path segments with energy-based coloring
				// Create gradient by adding multiple line segments with interpolated colors
				const int gradient_segments = 10;
				for (int i = 0; i < gradient_segments; i++) {
					float t1 = static_cast<float>(i) / static_cast<float>(gradient_segments);
					float t2 = static_cast<float>(i + 1) / static_cast<float>(gradient_segments);

					glm::vec3 seg_start = start + t1 * (end - start);
					glm::vec3 seg_end = start + t2 * (end - start);

					// Always draw the segment - don't check geometry intersection as it causes gaps
					glm::vec4 seg_color = start_color * (1.0f - (t1 + t2) * 0.5f) + end_color * ((t1 + t2) * 0.5f);
					seg_color.a = 1.0f;

					add_line(seg_start, seg_end, seg_color);
				}

				// Move to next segment
				current = next;
				next = current->next;
			}

			// Also draw emitted paths if they exist with adaptive energy-based coloring
			current = photon.path_head;
			while (current) {
				if (current->emit) {
					// Use the ORIGINAL emittance direction as computed by the simulator
					// The simulator should already compute the correct scattered/reflected direction
					float emit_energy = static_cast<float>(current->value);
					glm::vec4 emit_color = adaptive_log_color(emit_energy);
					emit_color.a = 1.0f;

					glm::vec3 start(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
									static_cast<float>(current->position.z));
					glm::vec3 emit_end(static_cast<float>(current->emit->position.x),
									   static_cast<float>(current->emit->position.y),
									   static_cast<float>(current->emit->position.z));
					
					// Always draw emitted paths - don't check geometry intersection
					add_line(start, emit_end, emit_color);
				}
				current = current->next;
			}
		}
	}

	end_lines();
	draw_lines();

	// Draw scatter/interaction markers with better visualization (user feedback)
	begin_points();

	for (const Photon& photon : simulator_->get_paths()) {
		if (photon.path_head) {
			auto current = photon.path_head;
			int vertex_count = 0;

			// Count total vertices to identify key points properly
			auto temp = current;
			while (temp) {
				vertex_count++;
				temp = temp->next;
			}

			// Only add markers at specific key points: incident, scatter, exit (user request 2)
			auto path_current = photon.path_head;
			std::shared_ptr<PhotonNode> prev = nullptr;
			std::shared_ptr<PhotonNode> next = nullptr;
			int current_index = 0;

			while (path_current) {
				glm::vec3 pos(static_cast<float>(path_current->position.x),
							  static_cast<float>(path_current->position.y),
							  static_cast<float>(path_current->position.z));

				bool should_mark = false;
				glm::vec4 marker_color;

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
					if (prev && next && simulator_) {
						// Get positions
						glm::vec3 prev_pos(static_cast<float>(prev->position.x), static_cast<float>(prev->position.y),
										   static_cast<float>(prev->position.z));
						glm::vec3 next_pos(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
										   static_cast<float>(next->position.z));

						// Check if this point represents a medium boundary crossing
						// by checking if we're at the surface (z ≈ 0) or at material boundaries
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
					add_point(pos, marker_color);
				}

				prev = path_current;
				path_current = path_current->next;
				current_index++;
			}
		}
	}

	end_points();
	draw_points();
}

void Renderer::draw_paths_instanced(const Settings& settings) {
	if (!simulator_)
		return;

	// PERFORMANCE OPTIMIZATION: Use cached energy range instead of expensive per-frame analysis
	update_cached_energy_range(settings);

	// Create adaptive logarithmic mapping function using cached values
	auto adaptive_log_color = [this](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, cached_min_energy_, cached_max_energy_);
	};

	// PERFORMANCE OPTIMIZATION: Use incremental caching instead of rebuilding every frame
	size_t current_photon_count = simulator_->get_paths().size();
	
	// Check if we need to rebuild cache completely or just add new photons
	if (!path_instances_cached_ || current_photon_count < cached_photon_count_) {
		// Complete rebuild needed (first time or photons were removed)
		cached_line_instances_.clear();
		cached_point_instances_.clear();
		cached_photon_count_ = 0;
		// Reset buffer upload flags when cache is invalidated
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;
	} else if (current_photon_count > cached_photon_count_) {
		// Incremental update - only process new photons
		// Keep existing cache, just mark buffers for re-upload
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;
	}
	
	// Only process photons if we have new ones to add
	if (current_photon_count > cached_photon_count_) {

		// PERFORMANCE: Use cached surface calculation - moved to class member for proper caching
		const auto& layers = simulator_->get_all_layers();
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

		// PERFORMANCE: Process only NEW photons (incremental caching)
		const auto paths = simulator_->get_paths();
		for (size_t i = cached_photon_count_; i < paths.size(); ++i) {
			const Photon& photon = paths[i];
		if (photon.path_head) {
			// Generate connected line segments with energy gradient
			auto current = photon.path_head;
			auto next = current ? current->next : nullptr;

			// Generate incident ray from source to material surface (cached surface)
			if (current && !simulator_->sources.empty()) {
				glm::vec3 first_interaction(static_cast<float>(current->position.x),
											static_cast<float>(current->position.y),
											static_cast<float>(current->position.z));

				const Source& source = simulator_->sources[0];
				glm::vec3 source_pos(static_cast<float>(source.origin.x), static_cast<float>(source.origin.y),
									 static_cast<float>(source.origin.z));
				glm::vec3 source_dir(static_cast<float>(source.direction.x), static_cast<float>(source.direction.y),
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
						cached_line_instances_.push_back({surface_entry, first_interaction, refracted_color, refracted_color});
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

				glm::vec3 start(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								static_cast<float>(current->position.z));
				glm::vec3 end(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
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

		// PERFORMANCE: Add all emitter direction vectors to instanced rendering cache with deduplication
		if (!simulator_->emitters.empty()) {
			// Group emitters by origin position and direction (with margin for similar directions)
			std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::shared_ptr<Emitter>>> direction_groups;
			const double POSITION_PRECISION = 100.0; // Group positions within 0.01 units
			const double DIRECTION_PRECISION = 50.0;  // Group directions within ~0.02 radians (~1.1 degrees)
			
			for (const auto& emitter : simulator_->emitters) {
				if (emitter->weight > 0.001) { // Only process significant emitters
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
				if (grouped_emitters.empty()) continue;
				
				// Use the first emitter for position and direction
				const auto& representative = grouped_emitters[0];
				glm::vec3 start_pos(
					static_cast<float>(representative->position.x),
					static_cast<float>(representative->position.y),
					static_cast<float>(representative->position.z)
				);
				
				glm::vec3 end_pos(
					static_cast<float>(representative->position.x + representative->direction.x * 0.05),
					static_cast<float>(representative->position.y + representative->direction.y * 0.05),
					static_cast<float>(representative->position.z + representative->direction.z * 0.05)
				);
				
				// Average the energy of all emitters in this group
				double total_weight = 0.0;
				for (const auto& emitter : grouped_emitters) {
					total_weight += emitter->weight;
				}
				
				// Use percentage-based coloring with averaged energy
				float surface_refraction = static_cast<float>(simulator_->get_combined_surface_refraction());
				float energy_percentage = static_cast<float>(total_weight / surface_refraction);
				glm::vec4 direction_color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
				direction_color.a = 0.8f; // Slightly transparent for distinction
				
				cached_line_instances_.push_back({start_pos, end_pos, direction_color, direction_color});
			}
		}

		// Add specular reflection emitter for incident photon's surface reflection
		// This is integrated into the direction grouping above to avoid duplication
		double specular_reflection = simulator_->get_combined_specular_reflection();
		if (specular_reflection > 0.0 && !simulator_->sources.empty()) {
			const Source& source = simulator_->sources[0];
			
			// Check if there are any regular emitters at the same position with similar direction
			bool found_similar_emitter = false;
			glm::dvec3 specular_pos = source.intersect;
			glm::dvec3 specular_dir = glm::normalize(source.specular_direction);
			
			const double POSITION_TOLERANCE = 0.01; // Same as POSITION_PRECISION above
			const double ANGLE_TOLERANCE = 0.02;    // Same as DIRECTION_PRECISION above
			
			for (const auto& emitter : simulator_->emitters) {
				if (emitter->weight <= 0.001) continue;
				
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
				glm::vec3 reflection_end(
					static_cast<float>(source.intersect.x + source.specular_direction.x * 0.1),
					static_cast<float>(source.intersect.y + source.specular_direction.y * 0.1),
					static_cast<float>(source.intersect.z + source.specular_direction.z * 0.1)
				);
				
				// Use energy-based coloring - map specular reflection energy to color
				// Use the actual specular reflection value as energy for better color mapping
				float specular_energy = static_cast<float>(specular_reflection);
				glm::vec4 specular_color = get_adaptive_energy_color(specular_energy, 0.0f, 1.0f);
				specular_color.a = 1.0f; // Full opacity for better visibility
				
				cached_line_instances_.push_back({reflection_start, reflection_end, specular_color, specular_color});
			}
		}

		// PERFORMANCE: Collect scatter points INSIDE cache block to eliminate per-frame processing
		// First, collect scatter points from photon paths
		for (const Photon& photon : simulator_->get_paths()) {
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
					glm::vec4 marker_color{};

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
						if (prev && next && simulator_) {
							// Get positions
							glm::vec3 prev_pos(static_cast<float>(prev->position.x), static_cast<float>(prev->position.y),
											   static_cast<float>(prev->position.z));
							glm::vec3 next_pos(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
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

		// PERFORMANCE: Collect emitter points INSIDE cache block  
		// Second, add emitter exit points (these are the accurate surface boundary points)
		if (simulator_ && !simulator_->emitters.empty()) {
			// Add emitter points to the point instances vector
			for (const auto& emitter : simulator_->emitters) {
				// Use emitter->position which contains the corrected surface intersection coordinates
				glm::vec3 exit_pos(
					static_cast<float>(emitter->position.x),
					static_cast<float>(emitter->position.y), 
					static_cast<float>(emitter->position.z)
				);
				
				// Use percentage-based coloring instead of absolute weights for consistency
				// Convert absolute weight to percentage of total energy budget
				float surface_refraction = static_cast<float>(simulator_->get_combined_surface_refraction());
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
			double surface_specular_reflection = simulator_->get_combined_specular_reflection();
			if (surface_specular_reflection > 0.0 && !simulator_->sources.empty()) {
				const Source& source = simulator_->sources[0];
				
				// Use actual source intersection point
				glm::vec3 surface_entry(static_cast<float>(source.intersect.x),
										static_cast<float>(source.intersect.y),
										static_cast<float>(source.intersect.z));

				// Use percentage-based coloring for surface reflection
				double surface_refraction = simulator_->get_combined_surface_refraction();
				float energy_percentage = static_cast<float>(surface_specular_reflection / surface_refraction);
				glm::vec4 surface_point_color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
				surface_point_color.a = 1.0f; // Full opacity for surface reflection point

				PointInstance surface_point;
				surface_point.position = surface_entry;
				surface_point.color = surface_point_color;
				surface_point.size = 8.0f; // Slightly larger for surface reflection point
				cached_point_instances_.push_back(surface_point);
			}
		}

		// Update cached photon count and mark cache as up to date
		cached_photon_count_ = current_photon_count;
		path_instances_cached_ = true;
		
	}

	// Render cached line instances (MASSIVE performance improvement)
	if (!cached_line_instances_.empty() && line_instanced_shader_program_) {
		glUseProgram(line_instanced_shader_program_);
		
		// PERFORMANCE: Use cached uniform location, compute MVP directly
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		glUniformMatrix4fv(line_instanced_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));
		
		// PERFORMANCE FIX: Only upload buffer when data has changed, not every frame!
		if (!line_buffer_uploaded_) {
			glBindBuffer(GL_ARRAY_BUFFER, line_instance_vbo_);
			glBufferData(GL_ARRAY_BUFFER, cached_line_instances_.size() * sizeof(LineInstance), 
						 cached_line_instances_.data(), GL_STATIC_DRAW);
			line_buffer_uploaded_ = true;
		}
		
		// Render using uploaded buffer
		glBindVertexArray(line_instanced_vao_);
		glDrawArraysInstanced(GL_LINES, 0, 2, static_cast<GLsizei>(cached_line_instances_.size()));
		glBindVertexArray(0);
		
		glUseProgram(0);
	}

	// Render cached scatter points and emitter points using instanced rendering
	if (!cached_point_instances_.empty()) {
		// PERFORMANCE: Minimize OpenGL state changes during camera movement
		// Set up blending state once before rendering points
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
		glUseProgram(point_instanced_shader_program_);
		
		// PERFORMANCE: Use cached uniform location, compute MVP directly
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		glUniformMatrix4fv(point_instanced_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

		// PERFORMANCE FIX: Only upload buffer when data has changed, not every frame!
		if (!point_buffer_uploaded_) {
			glBindBuffer(GL_ARRAY_BUFFER, point_instance_vbo_);
			glBufferData(GL_ARRAY_BUFFER, cached_point_instances_.size() * sizeof(PointInstance), 
						 cached_point_instances_.data(), GL_STATIC_DRAW);
			point_buffer_uploaded_ = true;
		}

		// Bind VAO and render using uploaded buffer
		glBindVertexArray(point_instanced_vao_);
		glDrawArraysInstanced(GL_POINTS, 0, 1, static_cast<GLsizei>(cached_point_instances_.size()));
		glBindVertexArray(0);
		
		glUseProgram(0);
		
		// PERFORMANCE: Restore OpenGL state once after rendering
		glDisable(GL_BLEND);
		glDisable(GL_PROGRAM_POINT_SIZE);
	}
}



// ========================================
// CONSOLIDATED SHADER-BASED RENDERING METHODS
// ========================================

void Renderer::begin_lines() {
	line_vertices_.clear();
}

void Renderer::add_line(const glm::vec3& start, const glm::vec3& end, const glm::vec4& color) {
	line_vertices_.push_back({start, color});
	line_vertices_.push_back({end, color});
}

void Renderer::end_lines() {
	if (line_vertices_.empty())
		return;

	glBindBuffer(GL_ARRAY_BUFFER, line_vbo_);
	glBufferData(GL_ARRAY_BUFFER, line_vertices_.size() * sizeof(LineVertex), line_vertices_.data(), GL_DYNAMIC_DRAW);
}

void Renderer::draw_lines() {
	if (line_vertices_.empty() || !line_shader_program_)
		return;

	glUseProgram(line_shader_program_);

	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(line_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

	glBindVertexArray(line_vao_);
	glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(line_vertices_.size()));
	glBindVertexArray(0);
}

void Renderer::begin_points() {
	point_vertices_.clear();
}

void Renderer::add_point(const glm::vec3& position, const glm::vec4& color) {
	point_vertices_.push_back({position, color});
}

void Renderer::end_points() {
	if (point_vertices_.empty())
		return;

	glBindBuffer(GL_ARRAY_BUFFER, point_vbo_);
	glBufferData(GL_ARRAY_BUFFER, point_vertices_.size() * sizeof(PointVertex), point_vertices_.data(),
				 GL_DYNAMIC_DRAW);
}

void Renderer::draw_points() {
	if (point_vertices_.empty() || !point_shader_program_)
		return;

	glUseProgram(point_shader_program_);

	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(point_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));
	glUniform1f(point_size_uniform_location_, 8.0f); // Smaller spheres as requested

	glBindVertexArray(point_vao_);
	glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(point_vertices_.size()));
	glBindVertexArray(0);
}

void Renderer::begin_triangles() {
	triangle_vertices_.clear();
}

void Renderer::add_triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec4& color) {
	triangle_vertices_.push_back({v1, color});
	triangle_vertices_.push_back({v2, color});
	triangle_vertices_.push_back({v3, color});
}

void Renderer::end_triangles() {
	if (triangle_vertices_.empty())
		return;

	glBindBuffer(GL_ARRAY_BUFFER, triangle_vbo_);
	glBufferData(GL_ARRAY_BUFFER, triangle_vertices_.size() * sizeof(TriangleVertex), triangle_vertices_.data(),
				 GL_DYNAMIC_DRAW);
}

void Renderer::draw_triangles() {
	// Call with empty clipping planes for normal rendering
	std::vector<glm::vec4> empty_planes;
	draw_triangles_with_clipping(empty_planes);
}

void Renderer::draw_triangles_with_clipping(const std::vector<glm::vec4>& clipping_planes) {
	if (triangle_vertices_.empty() || !triangle_shader_program_)
		return;

	glUseProgram(triangle_shader_program_);

	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(triangle_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));

	// Set clipping plane uniforms using cached locations
	int num_planes = std::min(static_cast<int>(clipping_planes.size()), 6);
	glUniform1i(triangle_num_planes_uniform_location_, num_planes);
	
	// Enable OpenGL clipping planes
	for (int i = 0; i < num_planes && i < 6; i++) {
		glEnable(GL_CLIP_DISTANCE0 + i);
	}
	
	if (num_planes > 0 && triangle_clip_planes_uniform_location_ != -1) {
		glUniform4fv(triangle_clip_planes_uniform_location_, num_planes, glm::value_ptr(clipping_planes[0]));
	}

	glBindVertexArray(triangle_vao_);
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(triangle_vertices_.size()));
	glBindVertexArray(0);
	
	// Disable clipping planes after rendering
	for (int i = 0; i < num_planes && i < 6; i++) {
		glDisable(GL_CLIP_DISTANCE0 + i);
	}
}

bool Renderer::setup_line_rendering() {
	// Create shader program
	std::string vertex_source = load_shader_source("shaders/lines.vert");
	std::string fragment_source = load_shader_source("shaders/lines.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load line shaders" << std::endl;
		return false;
	}

	line_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!line_shader_program_) {
		return false;
	}

	// PERFORMANCE: Cache uniform location to avoid glGetUniformLocation every frame
	line_mvp_uniform_location_ = glGetUniformLocation(line_shader_program_, "uMVP");

	// Create VAO and VBO
	glGenVertexArrays(1, &line_vao_);
	glGenBuffers(1, &line_vbo_);

	glBindVertexArray(line_vao_);
	glBindBuffer(GL_ARRAY_BUFFER, line_vbo_);

	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), (void*)0);
	glEnableVertexAttribArray(0);

	// Color attribute (location 1)
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(LineVertex), (void*)offsetof(LineVertex, color));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
	return true;
}

bool Renderer::setup_point_rendering() {
	// Create shader program
	std::string vertex_source = load_shader_source("shaders/points.vert");
	std::string fragment_source = load_shader_source("shaders/points.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load point shaders" << std::endl;
		return false;
	}

	point_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!point_shader_program_) {
		return false;
	}

	// PERFORMANCE: Cache uniform locations to avoid glGetUniformLocation every frame
	point_mvp_uniform_location_ = glGetUniformLocation(point_shader_program_, "uMVP");
	point_size_uniform_location_ = glGetUniformLocation(point_shader_program_, "uPointSize");

	// Create VAO and VBO
	glGenVertexArrays(1, &point_vao_);
	glGenBuffers(1, &point_vbo_);

	glBindVertexArray(point_vao_);
	glBindBuffer(GL_ARRAY_BUFFER, point_vbo_);

	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)0);
	glEnableVertexAttribArray(0);

	// Color attribute (location 1)
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, color));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
	return true;
}

bool Renderer::setup_triangle_rendering() {
	// Create shader program
	std::string vertex_source = load_shader_source("shaders/triangles.vert");
	std::string fragment_source = load_shader_source("shaders/triangles.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load triangle shaders" << std::endl;
		return false;
	}

	triangle_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!triangle_shader_program_) {
		return false;
	}

	// PERFORMANCE: Cache uniform locations to avoid glGetUniformLocation every frame
	triangle_mvp_uniform_location_ = glGetUniformLocation(triangle_shader_program_, "uMVP");
	triangle_num_planes_uniform_location_ = glGetUniformLocation(triangle_shader_program_, "uNumClipPlanes");
	triangle_clip_planes_uniform_location_ = glGetUniformLocation(triangle_shader_program_, "uClipPlanes");

	// Create VAO and VBO
	glGenVertexArrays(1, &triangle_vao_);
	glGenBuffers(1, &triangle_vbo_);

	glBindVertexArray(triangle_vao_);
	glBindBuffer(GL_ARRAY_BUFFER, triangle_vbo_);

	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(TriangleVertex), (void*)0);
	glEnableVertexAttribArray(0);

	// Color attribute (location 1)
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(TriangleVertex), (void*)offsetof(TriangleVertex, color));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
	return true;
}

GLuint Renderer::create_shader_program(const std::string& vertex_source, const std::string& fragment_source) {
	GLuint vertex_shader = compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader)
			glDeleteShader(vertex_shader);
		if (fragment_shader)
			glDeleteShader(fragment_shader);
		return 0;
	}

	GLuint program = glCreateProgram();
	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	GLint success;
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success) {
		std::array<char, 512> info_log{};
		glGetProgramInfoLog(program, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Program linking failed: " << info_log.data() << std::endl;
		glDeleteProgram(program);
		program = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	return program;
}

GLuint Renderer::compile_shader(const std::string& source, GLenum shader_type) {
	GLuint shader = glCreateShader(shader_type);
	const char* source_cstr = source.c_str();
	glShaderSource(shader, 1, &source_cstr, nullptr);
	glCompileShader(shader);

	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		std::array<char, 512> info_log{};
		glGetShaderInfoLog(shader, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "Shader compilation failed: " << info_log.data() << std::endl;
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

std::string Renderer::load_shader_source(const std::string& file_path) {
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open shader file: " << file_path << std::endl;
		return "";
	}

	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

void Renderer::auto_manage_energy_labels(Settings& settings) {
	if (!simulator_) return;
	
	static bool auto_disabled_labels = false;
	bool many_photons = (simulator_->get_paths().size() > 10);
	
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

	if (!simulator_)
		return;

	// Collect energy labels from photon paths
	for (const Photon& photon : simulator_->get_paths()) {
		if (photon.path_head) {
			auto current = photon.path_head;
			int vertex_count = 0;

			// Count vertices
			auto temp = current;
			while (temp) {
				vertex_count++;
				temp = temp->next;
			}

			// Add labels for key points
			int current_index = 0;
			while (current) {
				bool should_label = false;
				std::string label_text;
				glm::vec4 label_color(1.0f, 1.0f, 1.0f, 1.0f);

				if (current_index == 0) {
					// Incidence point - check if at surface (reflected) or inside (transmitted into material)
					should_label = true;
					float energy_percent = static_cast<float>(current->value * 100.0);

					glm::vec3 pos(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								  static_cast<float>(current->position.z));

					// Check if at or near surface (incident interface)
					bool at_surface = (std::abs(pos.z) < 0.01f); // Within 0.01 units of surface

					if (at_surface) {
						// At incident interface - this is reflected light
						if (energy_percent < 1.0f) {
							label_text = "Reflected: <1%";
						}
						else {
							label_text = std::format("Reflected: {}%", static_cast<int>(energy_percent));
						}
					}
					else {
						// Inside material - this is just the incident energy
						if (energy_percent < 1.0f) {
							label_text = "<1%";
						}
						else {
							label_text = std::format("{}%", static_cast<int>(energy_percent));
						}
					}
					label_color = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f); // White for incident
				}
				else if (current_index == vertex_count - 1) {
					// End point - check if at bottom surface (transmitted) or internal (absorbed/scattered)
					should_label = true;
					float energy_percent = static_cast<float>(current->value * 100.0);

					glm::vec3 pos(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								  static_cast<float>(current->position.z));

					// Calculate material bottom boundary
					float bottom_y = -0.1f; // Default bottom
					const auto& layers = simulator_->get_all_layers();
					if (!layers.empty()) {
						bottom_y = 1000.0f; // Start very high
						for (const auto& layer : layers) {
							for (const auto& triangle : layer.mesh) {
								float y0 = static_cast<float>(triangle.v0().y);
								float y1 = static_cast<float>(triangle.v1().y);
								float y2 = static_cast<float>(triangle.v2().y);
								if (y0 < bottom_y)
									bottom_y = y0;
								if (y1 < bottom_y)
									bottom_y = y1;
								if (y2 < bottom_y)
									bottom_y = y2;
							}
						}
					}

					// Calculate material top boundary as well
					float top_y = 0.1f;   // Default top
					if (!layers.empty()) {
						top_y = -1000.0f; // Start very low
						for (const auto& layer : layers) {
							for (const auto& triangle : layer.mesh) {
								float y0 = static_cast<float>(triangle.v0().y);
								float y1 = static_cast<float>(triangle.v1().y);
								float y2 = static_cast<float>(triangle.v2().y);
								if (y0 > top_y)
									top_y = y0;
								if (y1 > top_y)
									top_y = y1;
								if (y2 > top_y)
									top_y = y2;
							}
						}
					}

					// Check if at interfaces
					bool at_top = (std::abs(pos.y - top_y) < 0.05f);       // Increased tolerance for surface detection
					bool at_bottom = (std::abs(pos.y - bottom_y) < 0.05f); // Increased tolerance for surface detection

					if (at_top) {
						// At incident interface - this is reflected light exiting back through the top
						if (energy_percent < 1.0f) {
							label_text = "Reflected: <1%";
						}
						else {
							label_text = std::format("Reflected: {}%", static_cast<int>(energy_percent));
						}
					}
					else if (at_bottom) {
						// At exit interface - this is transmitted light
						if (energy_percent < 1.0f) {
							label_text = "Transmitted: <1%";
						}
						else {
							label_text = std::format("Transmitted: {}%", static_cast<int>(energy_percent));
						}
					}
					else {
						// Inside material - just show percentage
						if (energy_percent < 1.0f) {
							label_text = "<1%";
						}
						else {
							label_text = std::format("{}%", static_cast<int>(energy_percent));
						}
					}

					// Use logarithmic color mapping like the spheres
					float energy = std::max(0.1f, static_cast<float>(current->value));
					float log_energy = std::log10(energy + 0.01f) + 2.0f;
					log_energy = std::clamp(log_energy / 3.0f, 0.0f, 1.0f);

					if (log_energy > 0.8f) {
						label_color = glm::vec4(1.0f, 1.0f, 0.9f, 1.0f);
					}
					else if (log_energy > 0.5f) {
						float t = (log_energy - 0.5f) / 0.3f;
						label_color = glm::vec4(1.0f, 0.6f + 0.4f * t, 0.3f * t, 1.0f);
					}
					else {
						float t = log_energy / 0.5f;
						label_color = glm::vec4(0.6f + 0.4f * t, 0.2f * t, 0.2f * t, 1.0f);
					}
				}
				else {
					// Check for scatter/boundary points
					if (current->emit || current_index > 0) {
						// Look for medium boundary crossings or significant energy drops
						glm::vec3 pos(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
									  static_cast<float>(current->position.z));

						// Check if at surface (z ≈ 0) or has emitted path
						bool is_boundary = (std::abs(pos.z) < 0.001f) || (current->emit != nullptr);

						if (is_boundary) {
							should_label = true;
							float energy_percent = static_cast<float>(current->value * 100.0);

							// Determine label type based on position
							// Calculate material boundaries
							float top_y = 0.1f, bottom_y = -0.1f;
							const auto& layers = simulator_->get_all_layers();
							if (!layers.empty()) {
								top_y = -1000.0f;
								bottom_y = 1000.0f;
								for (const auto& layer : layers) {
									for (const auto& triangle : layer.mesh) {
										float y0 = static_cast<float>(triangle.v0().y);
										float y1 = static_cast<float>(triangle.v1().y);
										float y2 = static_cast<float>(triangle.v2().y);
										if (y0 > top_y)
											top_y = y0;
										if (y1 > top_y)
											top_y = y1;
										if (y2 > top_y)
											top_y = y2;
										if (y0 < bottom_y)
											bottom_y = y0;
										if (y1 < bottom_y)
											bottom_y = y1;
										if (y2 < bottom_y)
											bottom_y = y2;
									}
								}
							}

							bool at_top = (std::abs(pos.y - top_y) < 0.05f);
							bool at_bottom = (std::abs(pos.y - bottom_y) < 0.05f);

							// Instead of using geometric position to guess reflection/transmission,
							// use the actual physics-based classification from the simulation.
							// The simulation now correctly classifies reflection vs transmission
							// based on proper directional analysis.
							
							if (at_top || at_bottom) {
								// For exit points, show the energy percentage without incorrect R/T labels
								// The actual R/T classification is shown in the overlay statistics
								if (energy_percent < 1.0f) {
									label_text = "Exit: <1%";
								}
								else {
									label_text = std::format("Exit: {}%", static_cast<int>(energy_percent));
								}
							}
							else {
								// Internal scattering - just percentage
								if (energy_percent < 1.0f) {
									label_text = "<1%";
								}
								else {
									label_text = std::format("{}%", static_cast<int>(energy_percent));
								}
							}

							// Color based on energy like the spheres
							float energy = std::max(0.1f, static_cast<float>(current->value));
							float log_energy = std::log10(energy + 0.01f) + 2.0f;
							log_energy = std::clamp(log_energy / 3.0f, 0.0f, 1.0f);

							if (log_energy > 0.8f) {
								label_color = glm::vec4(1.0f, 1.0f, 0.9f, 1.0f);
							}
							else if (log_energy > 0.5f) {
								float t = (log_energy - 0.5f) / 0.3f;
								label_color = glm::vec4(1.0f, 0.6f + 0.4f * t, 0.3f * t, 1.0f);
							}
							else {
								float t = log_energy / 0.5f;
								label_color = glm::vec4(0.6f + 0.4f * t, 0.2f * t, 0.2f * t, 1.0f);
							}
						}
					}
				}

				if (should_label) {
					EnergyLabel label;
					label.world_position =
						glm::vec3(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								  static_cast<float>(current->position.z));
					label.text = label_text;
					label.color = label_color;
					label.scale = 1.0f; // Base scale

					cached_energy_labels_.push_back(label);
				}

				current = current->next;
				current_index++;
			}
		}
	}

	// Add surface scattering label at incident point (unified with other energy labels)
	double specular_reflection = simulator_->get_combined_specular_reflection();
	if (specular_reflection > 0.0) {
		// Calculate incident surface position
		glm::vec3 source_pos(0.05f, 0.0f, 0.05f); // From the incident ray calculations
		glm::vec3 source_dir(0.0f, 1.0f, 0.0f);   // Upward direction

		// Find topmost surface point
		float surface_y = 0.1f;
		const auto& layers = simulator_->get_all_layers();
		if (!layers.empty()) {
			surface_y = -1000.0f;
			for (const auto& layer : layers) {
				for (const auto& triangle : layer.mesh) {
					float y0 = static_cast<float>(triangle.v0().y);
					float y1 = static_cast<float>(triangle.v1().y);
					float y2 = static_cast<float>(triangle.v2().y);
					if (y0 > surface_y)
						surface_y = y0;
					if (y1 > surface_y)
						surface_y = y1;
					if (y2 > surface_y)
						surface_y = y2;
				}
			}
		}

		// Calculate surface entry point
		if (source_dir.y != 0.0f) {
			float t = (surface_y - source_pos.y) / source_dir.y;
			glm::vec3 surface_entry = source_pos + t * source_dir;

			// Add surface scattering label
			float surface_scattering = static_cast<float>(simulator_->get_combined_specular_reflection());
			float scatter_percent = (surface_scattering / static_cast<float>(simulator_->get_combined_surface_refraction())) * 100.0f;
			
			// Use percentage for consistent coloring with other energy labels
			float energy_percentage = scatter_percent / 100.0f; // Convert back to 0-1 range

			std::string label_text;
			if (scatter_percent < 1.0f) {
				label_text = "Reflected: <1%";
			}
			else {
				label_text = std::format("Reflected: {}%", static_cast<int>(scatter_percent));
			}

			EnergyLabel surface_label;
			surface_label.world_position = surface_entry;
			surface_label.text = label_text;
			surface_label.color = get_adaptive_energy_color(energy_percentage, 0.0f, 1.0f);
			surface_label.scale = 1.0f;

			cached_energy_labels_.push_back(surface_label);
		}
	}

	energy_labels_cached_ = true;
}

void Renderer::cache_energy_labels_from_emitters() {
	cached_energy_labels_.clear();
	energy_labels_cached_ = false;

	if (!simulator_)
		return;

	// Use emitter data for accurate energy labels with proper classification
	const auto& emitters = simulator_->emitters;
	
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
	const double DENSITY_RADIUS = 0.1; // Radius for density calculation
	const double DENSITY_AREA = std::numbers::pi * DENSITY_RADIUS * DENSITY_RADIUS; // Circle area
	
	for (const auto& emitter : emitters) {
		if (emitter->weight < 0.001) continue; // Skip very low energy exits
		
		int nearby_count = 0;
		for (const auto& other : emitters) {
			if (other->weight < 0.001) continue;
			
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
		if (emitter->weight < 0.001) continue; // Skip very low energy exits
		
		// Determine grouping multiplier based on density percentile
		double density = emitter_density[emitter];
		double grouping_multiplier;
		
		if (density >= density_90th) {
			grouping_multiplier = 5.0;   // Large groups for top 10% density (0.2 units)
		} else if (density >= density_75th) {
			grouping_multiplier = 10.0;  // Medium groups for top 25% density (0.1 units)  
		} else if (density >= density_50th) {
			grouping_multiplier = 20.0;  // Small groups for above median density (0.05 units)
		} else {
			grouping_multiplier = 50.0;  // Fine granularity for below median density (0.02 units)
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
		if (emitters_at_pos.empty()) continue;
		
		// Use the first emitter's position as the label position
		const auto& representative = emitters_at_pos[0];
		glm::vec3 label_pos(
			static_cast<float>(representative->position.x),
			static_cast<float>(representative->position.y),
			static_cast<float>(representative->position.z)
		);
		
		// Sum up energy and classify by exit type
		double total_reflected_energy = 0.0;
		double total_transmitted_energy = 0.0;
		double total_unclassified_energy = 0.0;
		
		for (const auto& emitter : emitters_at_pos) {
			switch (emitter->exit_type) {
				case Emitter::ExitType::REFLECTED:
					total_reflected_energy += emitter->weight;
					break;
				case Emitter::ExitType::TRANSMITTED:
					total_transmitted_energy += emitter->weight;
					break;
				default:
					total_unclassified_energy += emitter->weight;
					break;
			}
		}
		
		// Create label with proper classification and NORMALIZED percentages
		std::string label_text;
		glm::vec4 label_color(1.0f, 1.0f, 1.0f, 1.0f);
		
		double total_energy = total_reflected_energy + total_transmitted_energy + total_unclassified_energy;
		
		// CRITICAL FIX: Normalize by total number of photons, not surface refraction
		// Each photon starts with weight 1.0, so total initial energy = num_photons
		double total_initial_energy = static_cast<double>(simulator_->photons.size());
		
		// Calculate normalized percentage (matches console energy conservation calculation)
		double energy_percent = (total_energy / total_initial_energy) * 100.0;
		
		// Format percentage, showing "<1%" instead of "0%" for very small values
		int rounded_percent = static_cast<int>(energy_percent);
		std::string percent_text = (rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);
		
		if (total_reflected_energy > total_transmitted_energy && total_reflected_energy > total_unclassified_energy) {
			// Predominantly reflected
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.8f, 0.2f, 1.0f); // Bright green for reflection (clearly visible)
		} else if (total_transmitted_energy > total_reflected_energy && total_transmitted_energy > total_unclassified_energy) {
			// Predominantly transmitted
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.6f, 1.0f, 1.0f); // Bright blue for transmission (clearly visible)
		} else {
			// Mixed or unclassified
			label_text = percent_text;
			label_color = glm::vec4(0.8f, 0.8f, 0.8f, 1.0f); // Gray for mixed/unclassified
		}
		
		cached_energy_labels_.push_back({
			.world_position = label_pos,
			.text = label_text,
			.color = label_color,
			.scale = 1.0f,
			.screen_position = glm::vec2(0.0f),
			.screen_position_valid = false
		});
	}

	// Add surface specular reflection label - position at tip of reflection vector
	double specular_reflection = simulator_->get_combined_specular_reflection();
	if (specular_reflection > 0.0 && !simulator_->sources.empty()) {
		const Source& source = simulator_->sources[0];
		
		// Position label at the tip of the specular reflection vector (twice as long as regular emitters)
		glm::vec3 surface_label_pos(
			static_cast<float>(source.intersect.x + source.specular_direction.x * 0.1),
			static_cast<float>(source.intersect.y + source.specular_direction.y * 0.1),
			static_cast<float>(source.intersect.z + source.specular_direction.z * 0.1)
		);
		
		// Calculate normalized percentage for surface reflection
		double surface_refraction = simulator_->get_combined_surface_refraction();
		double energy_percent = (specular_reflection / surface_refraction) * 100.0;
		
		// Format percentage, showing "<1%" instead of "0%" for very small values
		int rounded_percent = static_cast<int>(energy_percent);
		std::string surface_percent_text = (rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);
		
		// Use bright purple color for surface specular reflection to match saturation of green and blue
		glm::vec4 surface_label_color = glm::vec4(0.8f, 0.3f, 1.0f, 1.0f);
		
		cached_energy_labels_.push_back({
			.world_position = surface_label_pos,
			.text = surface_percent_text,
			.color = surface_label_color,
			.scale = 1.0f,
			.screen_position = glm::vec2(0.0f),
			.screen_position_valid = false
		});
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
		cache_energy_labels_from_emitters();
	}

	// Render cached text labels using pre-calculated screen positions
	if (text_render_callback_) {
		for (const auto& label : cached_energy_labels_) {
			if (label.screen_position_valid) {
				text_render_callback_(label.text, label.screen_position.x, label.screen_position.y, label.color);
			}
		}
	}
	else {
		// Fallback: render colored points where text would appear
		begin_points();
		for (const auto& label : cached_energy_labels_) {
			if (label.screen_position_valid) {
				// Offset the label position slightly above the actual point to avoid overlapping
				glm::vec3 label_pos = label.world_position + glm::vec3(0.0f, 0.0f, 0.05f);
				add_point(label_pos, label.color);
			}
		}
		end_points();
		draw_points();
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
		float b = 1.0f - t * 0.4f;  // From white (1.0) to very bright yellow (0.6)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.80f) {
		// 80-85% energy: very bright yellow to bright yellow
		float t = (normalized - 0.80f) / 0.05f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.6f - t * 0.2f;  // From very bright yellow (0.6) to bright yellow (0.4)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.70f) {
		// 70-80% energy: bright yellow to yellow
		float t = (normalized - 0.70f) / 0.10f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.4f - t * 0.2f;   // From bright yellow (0.4) to yellow (0.2)
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

glm::vec4 Renderer::get_layer_specific_energy_color(float energy, float min_energy, float max_energy,
													uint8_t tissue_id) {
	// Clamp and normalize energy like the original function
	energy = std::clamp(energy, min_energy, max_energy);
	float linear_normalized = (energy - min_energy) / (max_energy - min_energy);
	linear_normalized = std::clamp(linear_normalized, 0.0f, 1.0f);
	float normalized = std::pow(linear_normalized, 0.3f);
	normalized = std::clamp(normalized, 0.0f, 1.0f);

	// Define base colors for different material types
	glm::vec3 base_colors[6] = {
		glm::vec3(1.0f, 0.4f, 0.4f), // Red theme for material 0
		glm::vec3(0.4f, 0.4f, 1.0f), // Blue theme for material 1
		glm::vec3(0.4f, 1.0f, 0.4f), // Green theme for material 2
		glm::vec3(1.0f, 0.4f, 1.0f), // Magenta theme for material 3
		glm::vec3(1.0f, 1.0f, 0.4f), // Yellow theme for material 4
		glm::vec3(0.4f, 1.0f, 1.0f), // Cyan theme for material 5
	};

	glm::vec3 base_color = base_colors[tissue_id % 6];

	// Create energy gradient within the material's color theme
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

void Renderer::begin_voxel_instances() {
	// DON'T clear instances every frame - this breaks camera change detection!
	// Only clear when we actually need to rebuild (when voxel_instances_dirty_ is true)
	if (voxel_instances_dirty_) {
		voxel_instances_.clear();
	}
}

void Renderer::add_voxel_instance(const glm::vec3& position, const glm::vec4& color, float scale, float depth) {
	voxel_instances_.push_back({position, color, scale, depth});
}

void Renderer::end_voxel_instances(VoxelMode mode) {
	if (voxel_instances_.empty()) return;
	
	// PERFORMANCE FIX: Only skip processing if truly nothing changed
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
				glm::length(current_camera_pos - cached_camera_position_) > SIGNIFICANT_MOVEMENT ||
				glm::length(current_camera_target - cached_camera_target_) > SIGNIFICANT_MOVEMENT;
			
			if (significant_camera_movement && voxel_instances_.size() < 100000) { // Limit background sorting to reasonable sizes
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
					std::sort(background_sorted_voxels_.begin(), background_sorted_voxels_.end(), 
						[](const VoxelInstance& a, const VoxelInstance& b) {
							return a.depth > b.depth; // Back to front
						});
					
					// Mark as ready for GPU upload
					background_sort_ready_.store(true);
					background_sort_in_progress_.store(false);
				});
			}
		}
		
		// Current frame continues with existing data - no performance impact
		return;
	}
	
	// Store mode for consistency
	current_voxel_mode_ = mode;
	
	// PERFORMANCE vs QUALITY TRADE-OFF:
	// Use consistent sorting behavior for all modes
	bool should_sort_for_quality = false;
	
	if (voxel_instances_.size() < 50000) {
		should_sort_for_quality = true; // Sort small datasets
	}
	
	if (should_sort_for_quality) {
		// Sort by depth (back to front) for optimal transparency
		std::sort(voxel_instances_.begin(), voxel_instances_.end(), 
			[](const VoxelInstance& a, const VoxelInstance& b) {
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
}void Renderer::draw_voxel_instances() {
	if (voxel_instances_.empty() || !voxel_shader_program_) return;
	
	// Enable transparency - KEEP ORIGINAL APPEARANCE
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// Disable depth writing for transparent objects but keep depth testing
	glDepthMask(GL_FALSE);
	glEnable(GL_DEPTH_TEST);
	
	glUseProgram(voxel_shader_program_);
	
	// Set MVP matrix using cached location for performance
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(voxel_mvp_uniform_location_, 1, GL_FALSE, glm::value_ptr(mvp));
	
	// Bind VAO and draw instanced
	glBindVertexArray(voxel_cube_vao_);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(voxel_instances_.size()));
	
	// Restore OpenGL state
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glBindVertexArray(0);
}

bool Renderer::setup_voxel_instanced_rendering() {
	// Load shaders
	std::string vertex_source = load_shader_source("shaders/voxels.vert");
	std::string fragment_source = load_shader_source("shaders/voxels.frag");
	
	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load voxel instanced shaders" << std::endl;
		return false;
	}
	
	voxel_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!voxel_shader_program_) {
		return false;
	}

	// PERFORMANCE: Cache uniform location to avoid glGetUniformLocation every frame
	voxel_mvp_uniform_location_ = glGetUniformLocation(voxel_shader_program_, "uMVP");
	
	// Create unit cube geometry (centered at origin)
	float vertices[] = {
		// Front face (z = 0.5)
		-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // bottom-left
		 0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // bottom-right
		 0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // top-right
		 0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // top-right
		-0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // top-left
		-0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, // bottom-left
		
		// Back face (z = -0.5)
		-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // bottom-left
		 0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // top-right
		 0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // bottom-right
		 0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // top-right
		-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // bottom-left
		-0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, // top-left
		
		// Left face (x = -0.5)
		-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f, // top-right
		-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f, // bottom-left
		-0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f, // top-left
		-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f, // bottom-left
		-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f, // top-right
		-0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f, // bottom-right
		
		// Right face (x = 0.5)
		 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f, // top-left
		 0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f, // top-right
		 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f, // bottom-right
		 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f, // bottom-right
		 0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f, // bottom-left
		 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f, // top-left
		
		// Bottom face (y = -0.5)
		-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, // top-right
		 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, // bottom-left
		 0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, // top-left
		 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, // bottom-left
		-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, // top-right
		-0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, // bottom-right
		
		// Top face (y = 0.5)
		-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f, // top-left
		 0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f, // top-right
		 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f, // bottom-right
		 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f, // bottom-right
		-0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f, // bottom-left
		-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f  // top-left
	};
	
	// Create VAO and VBO for cube geometry
	glGenVertexArrays(1, &voxel_cube_vao_);
	glGenBuffers(1, &voxel_cube_vbo_);
	glGenBuffers(1, &voxel_instance_vbo_);
	
	glBindVertexArray(voxel_cube_vao_);
	
	// Setup cube geometry buffer
	glBindBuffer(GL_ARRAY_BUFFER, voxel_cube_vbo_);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	
	// Position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	
	// Normal attribute (location 1)  
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	
	// Setup instance data buffer
	glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);
	
	// Instance position (location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, position));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1); // One per instance
	
	// Instance color (location 3)
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1); // One per instance
	
	// Instance scale (location 4)
	glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(VoxelInstance), (void*)offsetof(VoxelInstance, scale));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1); // One per instance
	
	glBindVertexArray(0);
	return true;
}

bool Renderer::setup_point_instanced_rendering() {
	// Load shaders
	std::string vertex_source = load_shader_source("shaders/points_instanced.vert");
	std::string fragment_source = load_shader_source("shaders/points_instanced.frag");
	
	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load point instanced shaders" << std::endl;
		return false;
	}
	
	point_instanced_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!point_instanced_shader_program_) {
		return false;
	}

	// PERFORMANCE: Cache uniform location to avoid glGetUniformLocation every frame
	point_instanced_mvp_uniform_location_ = glGetUniformLocation(point_instanced_shader_program_, "uMVP");

	// Create and bind VAO for instanced point rendering
	glGenVertexArrays(1, &point_instanced_vao_);
	glBindVertexArray(point_instanced_vao_);

	// Generate VBOs for point vertex data and instance data
	glGenBuffers(1, &point_instanced_vbo_);
	glGenBuffers(1, &point_instance_vbo_);

	// Set up point vertex data (single point at origin)
	glBindBuffer(GL_ARRAY_BUFFER, point_instanced_vbo_);
	float point_vertex[3] = {0.0f, 0.0f, 0.0f}; // Single point at origin
	glBufferData(GL_ARRAY_BUFFER, sizeof(point_vertex), point_vertex, GL_STATIC_DRAW);

	// Vertex position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Set up instance data buffer (initially empty)
	glBindBuffer(GL_ARRAY_BUFFER, point_instance_vbo_);
	
	// Instance position attribute (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, position));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1); // One per instance

	// Instance color attribute (location 2)
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, color));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1); // One per instance

	// Instance size attribute (location 3)
	glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, size));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1); // One per instance

	glBindVertexArray(0);
	return true;
}

// Performance optimization: Cache energy range calculation to avoid expensive per-frame analysis
void Renderer::update_cached_energy_range(const Settings& settings) const {
	if (!simulator_) {
		energy_range_cached_ = false;
		return;
	}

	// PERFORMANCE FIX: Use proper cache invalidation instead of checking every frame
	if (energy_range_cached_) {
		return; // Use cached values - they're already invalidated by simulation version changes
	}

	// Recalculate energy range using optimized method
	std::vector<float> all_energies;
	
	// Collect energies from photon paths
	for (const Photon& photon : simulator_->get_paths()) {
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

	// Collect energies from voxels using the ORIGINAL method for accurate color mapping
	const auto& config = Config::get();
	if (config.nx() > 0 && config.ny() > 0 && config.nz() > 0) {
		int nx = config.nx();
		int ny = config.ny();
		int nz = config.nz();
		
		// Use full scan like the original to maintain color accuracy
		for (int iz = 0; iz < nz; iz++) {
			for (int iy = 0; iy < ny; iy++) {
				for (int ix = 0; ix < nx; ix++) {
					if (ix >= 0 && iy >= 0 && iz >= 0 && 
						static_cast<uint32_t>(ix) < static_cast<uint32_t>(nx) && 
						static_cast<uint32_t>(iy) < static_cast<uint32_t>(ny) && 
						static_cast<uint32_t>(iz) < static_cast<uint32_t>(nz)) {
						
						Voxel* voxel = simulator_->voxel_grid(static_cast<uint32_t>(ix), 
															 static_cast<uint32_t>(iy), 
															 static_cast<uint32_t>(iz));
						if (voxel) {
							float absorption = static_cast<float>(voxel->absorption);
							float emittance = static_cast<float>(voxel->total_emittance()); // Includes specular reflection
							
							// Use the current settings to determine energy calculation method
							float total_energy;
							if (settings.voxel_mode == VoxelMode::Absorption) {
								total_energy = absorption;
							}
							else if (settings.voxel_mode == VoxelMode::Emittance) {
								total_energy = emittance;
							}
							else if (settings.voxel_mode == VoxelMode::Layers) {
								total_energy = (voxel->material != nullptr) ? 1.0f : 0.0f;
							}
							else {
								total_energy = absorption + emittance;
							}

							if (total_energy > 0.0000001f) {
								all_energies.push_back(total_energy);
							}
						}
					}
				}
			}
		}
	}

	// Calculate range using the ORIGINAL percentile-based method for accurate colors
	cached_min_energy_ = 1e-8f; // Much lower threshold for emittance detection
	cached_max_energy_ = 0.01f; // Default fallback

	if (!all_energies.empty()) {
		// Use the ORIGINAL percentile-based scaling for better contrast
		std::ranges::sort(all_energies);

		// Use percentile-based scaling for better contrast (SAME as original)
		const size_t count = all_energies.size();
		const float p5 = all_energies[static_cast<size_t>(count * 0.05)];  // 5th percentile
		const float p95 = all_energies[static_cast<size_t>(count * 0.95)]; // 95th percentile

		cached_min_energy_ = std::max(p5, 1e-8f);                            // Lower minimum
		cached_max_energy_ = std::max(p95, cached_min_energy_ * 10.0f);      // Ensure reasonable range
	}

	// Update cache state
	energy_range_cached_ = true;
}

// Invalidate energy cache when simulation data changes
void Renderer::invalidate_energy_cache() const {
	energy_range_cached_ = false;
}

// PERFORMANCE OPTIMIZATION: Comprehensive cache invalidation
void Renderer::invalidate_all_caches() {
	// Invalidate energy-related caches
	invalidate_energy_cache();
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
}

// PERFORMANCE: Smart settings update - only invalidate path cache when path-related settings change
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
	// Load shaders
	std::string vertex_source = load_shader_source("shaders/lines_instanced.vert");
	std::string fragment_source = load_shader_source("shaders/lines_instanced.frag");
	
	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "Failed to load line instanced shaders" << std::endl;
		return false;
	}
	
	line_instanced_shader_program_ = create_shader_program(vertex_source, fragment_source);
	if (!line_instanced_shader_program_) {
		return false;
	}
	
	// PERFORMANCE: Cache uniform location to avoid glGetUniformLocation every frame
	line_instanced_mvp_uniform_location_ = glGetUniformLocation(line_instanced_shader_program_, "uMVP");
	
	// Create line geometry (just two endpoints: 0 and 1)
	float line_vertices[] = {
		0.0f, 0.0f, 0.0f, // Start point (will be interpolated with instance data)
		1.0f, 0.0f, 0.0f  // End point (will be interpolated with instance data)
	};
	
	// Create VAO and VBO for line geometry
	glGenVertexArrays(1, &line_instanced_vao_);
	glGenBuffers(1, &line_instanced_vbo_);
	glGenBuffers(1, &line_instance_vbo_);
	
	glBindVertexArray(line_instanced_vao_);
	
	// Set up line geometry (positions only)
	glBindBuffer(GL_ARRAY_BUFFER, line_instanced_vbo_);
	glBufferData(GL_ARRAY_BUFFER, sizeof(line_vertices), line_vertices, GL_STATIC_DRAW);
	
	// Position attribute (location 0) - just x component used for interpolation
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	
	// Set up instance data buffer (empty for now)
	glBindBuffer(GL_ARRAY_BUFFER, line_instance_vbo_);
	glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
	
	// Instance start position (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1); // One per instance
	
	// Instance end position (location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1); // One per instance
	
	// Instance start color (location 3)
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start_color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1); // One per instance
	
	// Instance end color (location 4)
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end_color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1); // One per instance
	
	glBindVertexArray(0);
	return true;
}

/**
 * Check if a point is inside the mesh geometry
 */
bool Renderer::is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const {
	glm::dvec3 dpoint(point.x, point.y, point.z);
	return simulator.is_point_inside_geometry(dpoint);
}
