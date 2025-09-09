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
#include "simulator/simulator.hpp"

Renderer::Renderer() = default;

Renderer::~Renderer() {
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
	std::cout << "Initializing modern OpenGL 4.5 renderer..." << std::endl;

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

	std::cout << "Consolidated renderer initialized successfully" << std::endl;
	return true;
}

void Renderer::setup_opengl() {
	// Query OpenGL version
	const GLubyte* version = glGetString(GL_VERSION);
	const GLubyte* renderer = glGetString(GL_RENDERER);
	std::cout << "OpenGL Version: " << version << std::endl;
	std::cout << "Renderer: " << renderer << std::endl;

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

	std::cout << "OpenGL 4.5 setup complete" << std::endl;
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
	
	// Invalidate energy cache if simulation data might have changed
	invalidate_energy_cache();

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
	for (const auto& layer : simulator.layers) {
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
	Range3 bounds = simulator.bounds;

	// Set camera target to center of bounds
	glm::vec3 center = to_float(bounds.center());
	camera_.set_target(center);

	// Set tissue elevation bounds for scroll wheel constraint (Y-axis)
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
	int nx = simulator_->config.nx();
	int ny = simulator_->config.ny();
	int nz = simulator_->config.nz();
	double voxsize = simulator_->config.vox_size();
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
						float emittance = static_cast<float>(voxel->emittance);

						float total_energy;
						if (settings.voxel_mode == VoxelMode::Absorption) {
							total_energy = absorption;
						}
						else if (settings.voxel_mode == VoxelMode::Emittance) {
							total_energy = emittance;
						}
						else if (settings.voxel_mode == VoxelMode::Layers) {
							// For Layers mode, use a constant value for all voxels with tissue
							total_energy = (voxel->tissue != nullptr) ? 1.0f : 0.0f;
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
						// Calculate center position using MCML's exact formula
						double x = simulator_->bounds.min_bounds.x + (voxsize * ix) + half_voxsize;
						double y = simulator_->bounds.min_bounds.y + (voxsize * iy) + half_voxsize;
						double z = simulator_->bounds.min_bounds.z + (voxsize * iz) + half_voxsize;

						// Check if this voxel is actually inside the mesh geometry AND has tissue
						bool is_in_geometry = simulator_->is_point_inside_geometry({x,y,z});
						bool has_tissue = (voxel->tissue != nullptr);
						
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
						
						// Only render voxels that have tissue AND intersect with the geometry
						if (!has_tissue || !voxel_intersects_geometry) {
							continue;
						}

						float fx = static_cast<float>(x);
						float fy = static_cast<float>(y);
						float fz = static_cast<float>(z);
						glm::vec3 voxel_pos(fx, fy, fz);

						// Energy-based coloring based on selected voxel mode (user request 4)
						float absorption = static_cast<float>(voxel->absorption);
						float emittance = static_cast<float>(voxel->emittance);

						// ADD SURFACE SCATTERING AS EMITTANCE at the surface entry voxel
						// Always detect the surface entry voxel, but only add emittance in Emittance mode
						if (!simulator_->sources.empty()) {
							const auto& source = simulator_->sources[0];
							glm::vec3 source_pos(static_cast<float>(source.origin.x),
												 static_cast<float>(source.origin.y),
												 static_cast<float>(source.origin.z));
							glm::vec3 source_dir(static_cast<float>(source.direction.x),
												 static_cast<float>(source.direction.y),
												 static_cast<float>(source.direction.z));

							// Calculate surface entry point
							float surface_y = 0.1f;
							if (!simulator_->layers.empty()) {
								surface_y = -1000.0f;
								for (const auto& layer : simulator_->layers) {
									for (const auto& triangle : layer.mesh) {
										surface_y = std::max({surface_y, static_cast<float>(triangle.v0().y),
															  static_cast<float>(triangle.v1().y),
															  static_cast<float>(triangle.v2().y)});
									}
								}
							}

							if (source_dir.y != 0.0f) {
								float t = (surface_y - source_pos.y) / source_dir.y;
								glm::vec3 surface_entry = source_pos + t * source_dir;

								// Check if this voxel contains the surface entry point
								float voxel_tolerance =
									1.5f * static_cast<float>(voxsize); // Increase tolerance for better detection
								float distance = glm::length(voxel_pos - surface_entry);

								if (distance < voxel_tolerance) {
									// Only add surface reflection as emittance in Emittance mode
									if (settings.voxel_mode == VoxelMode::Emittance) {
										float surface_scattering =
											static_cast<float>(simulator_->record.specular_reflection);
										emittance += surface_scattering;
									}
								}
							}
						}

						float total_energy;
						if (settings.voxel_mode == VoxelMode::Absorption) {
							total_energy = absorption;
						}
						else if (settings.voxel_mode == VoxelMode::Emittance) {
							total_energy = emittance;
						}
						else if (settings.voxel_mode == VoxelMode::Layers) {
							// For Layers mode, we ignore energy and just show tissue types
							total_energy = 1.0f;                   // Constant energy for all voxels with tissue
						}
						else {
							total_energy = absorption + emittance; // Fallback to combined
						}

						// Show even the tiniest energy interactions
						glm::vec4 color(0.0f, 0.0f, 0.0f, 0.0f); // Default transparent

						// Check if this voxel is actually inside the mesh geometry
						glm::dvec3 voxel_center = glm::dvec3(x, y, z);
						bool is_inside_geometry = simulator_->is_point_inside_geometry(voxel_center);

						// In Layers mode, show all voxels with tissue regardless of energy
						bool should_render_voxel = false;
						if (settings.voxel_mode == VoxelMode::Layers) {
							should_render_voxel = (voxel->tissue != nullptr) && is_inside_geometry;
						}
						else {
							should_render_voxel = (total_energy > 1e-8f || voxel->tissue != nullptr) && is_inside_geometry;
						}

						if (should_render_voxel) {
							// For Layers mode, use simplified rendering logic
							if (settings.voxel_mode == VoxelMode::Layers) {
								if (voxel->tissue != nullptr) {
									// Match the appearance of Absorption mode when there's no absorption
									// Use the same energy and alpha values as "no energy" voxels in absorption mode
									color = get_layer_specific_energy_color(1e-8f, min_energy, max_energy,
																			voxel->tissue->id);
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
								float max_dist = glm::length(to_float(simulator_->bounds.max_bounds));
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
																			voxel->tissue->id);
									color.a = alpha;
								}
								else if (voxel->tissue != nullptr) {
									// Voxel in medium but no recorded energy: very faint tissue-colored hint
									color = get_layer_specific_energy_color(1e-8f, min_energy, max_energy,
																			voxel->tissue->id);
									color.a = 0.05f; // Slightly more visible
								}
							}

							// Calculate distance to camera for depth sorting
							float distance = glm::length(voxel_pos - camera_pos);

							// Add to render list (only if it has visible color)
							voxels_to_render.push_back({voxel, voxel_pos, distance, color});
						}
					}
				}
			}
		}
	}

	// Sort voxels back-to-front for proper transparency blending (user request 4) - Modern C++20
	std::ranges::sort(voxels_to_render, [](const VoxelRenderData& a, const VoxelRenderData& b) noexcept {
		return a.distance_to_camera > b.distance_to_camera; // Furthest first
	});

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

	// Setup OpenGL state for transparent voxel rendering (match original)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE); // Disable depth writing for transparency
	glDisable(GL_CULL_FACE);

	// Use HIGH-PERFORMANCE instanced rendering instead of individual triangles
	begin_voxel_instances();

	// Get MCML grid parameters
	int nx = simulator_->config.nx();
	int ny = simulator_->config.ny();
	int nz = simulator_->config.nz();
	double voxsize = simulator_->config.vox_size();
	double half_voxsize = voxsize * 0.5;
	float voxel_scale = static_cast<float>(voxsize * 0.95); // Slightly smaller to show gaps

	// PERFORMANCE OPTIMIZATION: Use cached energy range instead of expensive per-frame analysis
	update_cached_energy_range(settings);

	// Set up instanced rendering for maximum performance
	begin_voxel_instances();
	
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
					if (!voxel || !voxel->tissue) continue; // Skip voxels with no tissue

					// Calculate energy based on current mode
					float absorption = static_cast<float>(voxel->absorption);
					float emittance = static_cast<float>(voxel->emittance);
					float total_energy;
					
					if (settings.voxel_mode == VoxelMode::Absorption) {
						total_energy = absorption;
					}
					else if (settings.voxel_mode == VoxelMode::Emittance) {
						total_energy = emittance;
					}
					else if (settings.voxel_mode == VoxelMode::Layers) {
						total_energy = 1.0f; // We already know voxel->tissue != nullptr
					}
					else {
						total_energy = absorption + emittance;
					}

					// IMPORTANT: Don't skip voxels with no energy! 
					// They should still render as very faint tissue-colored hints

					// Calculate world position using the same method as original code
					double x = simulator_->bounds.min_bounds.x + (voxsize * ix) + half_voxsize;
					double y = simulator_->bounds.min_bounds.y + (voxsize * iy) + half_voxsize;
					double z = simulator_->bounds.min_bounds.z + (voxsize * iz) + half_voxsize;

					glm::vec3 voxel_pos(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));

					// PERFORMANCE OPTIMIZATION: Skip expensive geometry tests
					// The voxel already has tissue assigned, so it's already been determined to be in geometry
					// during the voxelization process. No need to re-test every frame.

					// EXACT ORIGINAL COLOR CALCULATION - DO NOT CHANGE!
					glm::vec4 color(0.0f);

					if (settings.voxel_mode == VoxelMode::Layers) {
						// Layers mode: show all tissue voxels regardless of energy
						if (voxel->tissue != nullptr) {
							// Match the appearance of Absorption mode when there's no absorption
							// Use the same energy and alpha values as "no energy" voxels in absorption mode
							color = get_layer_specific_energy_color(1e-8f, cached_min_energy_, cached_max_energy_, voxel->tissue->id);
							color.a = 0.05f; // Same very faint alpha as no-absorption voxels
						}
					} else {
						// Original energy-based rendering logic for other modes

						// Use gamma correction for better mid-range contrast
						float normalized_energy = (total_energy - cached_min_energy_) / (cached_max_energy_ - cached_min_energy_);
						normalized_energy = std::clamp(normalized_energy, 0.0f, 1.0f);

						// Apply gamma correction for better contrast (gamma = 0.5 brightens mid-tones)
						float gamma_corrected = std::pow(normalized_energy, 0.5f);

						// Calculate distance from origin for depth-based alpha
						float max_dist = glm::length(to_float(simulator_->bounds.max_bounds));
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

						float alpha = min_alpha + (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.2f * normalized_dist);

						// Use layer-specific energy color mapping with dynamic range
						if (total_energy > 1e-8f) { // Lower threshold
							color = get_layer_specific_energy_color(total_energy, cached_min_energy_, cached_max_energy_, voxel->tissue->id);
							color.a = alpha;
						}
						else if (voxel->tissue != nullptr) {
							// Voxel in medium but no recorded energy: very faint tissue-colored hint
							color = get_layer_specific_energy_color(1e-8f, cached_min_energy_, cached_max_energy_, voxel->tissue->id);
							color.a = 0.05f; // Slightly more visible
						}
					}
					// Only add voxels that should be visible
					if (color.a > 0.0f) {
						// Calculate depth for sorting (distance from camera)
						glm::vec3 camera_pos = camera_.get_position();
						glm::vec3 camera_target = camera_.get_target();
						glm::vec3 camera_front = glm::normalize(camera_target - camera_pos);
						glm::vec3 view_dir = voxel_pos - camera_pos;
						float depth = glm::dot(view_dir, camera_front);
						
						add_voxel_instance(voxel_pos, color, voxel_scale, depth);
					}

				}
			}
		}
	}

	// Render ALL voxels in a single high-performance instanced draw call
	end_voxel_instances();
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
	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			auto current = path.head;
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

	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			// Draw connected line segments with energy gradient exactly like backup
			auto current = path.head;
			auto next = current ? current->next : nullptr;

			// Draw physically correct incident ray from source to tissue surface
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
				if (!simulator_->layers.empty()) {
					// Find the highest Y coordinate among all layers by examining triangle vertices
					surface_y = -1000.0f; // Start very low
					for (const auto& layer : simulator_->layers) {
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

					// If photon actually enters tissue, draw refracted ray from surface to first interaction
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
			current = path.head;
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

	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			auto current = path.head;
			int vertex_count = 0;

			// Count total vertices to identify key points properly
			auto temp = current;
			while (temp) {
				vertex_count++;
				temp = temp->next;
			}

			// Only add markers at specific key points: incident, scatter, exit (user request 2)
			auto path_current = path.head;
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
						// by checking if we're at the surface (z ≈ 0) or at tissue boundaries
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

	line_instances_.clear();
	std::vector<PointInstance> point_instances;

	// PERFORMANCE OPTIMIZATION: Use cached energy range instead of expensive per-frame analysis
	update_cached_energy_range(settings);

	// Create adaptive logarithmic mapping function using cached values
	auto adaptive_log_color = [this](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, cached_min_energy_, cached_max_energy_);
	};

	// PERFORMANCE: Cache surface calculation instead of computing every frame
	static float cached_surface_y = 0.1f;
	static bool surface_cached = false;
	
	if (!surface_cached && !simulator_->layers.empty()) {
		cached_surface_y = -1000.0f;
		for (const auto& layer : simulator_->layers) {
			for (const auto& triangle : layer.mesh) {
				cached_surface_y = std::max(cached_surface_y, static_cast<float>(triangle.v0().y));
				cached_surface_y = std::max(cached_surface_y, static_cast<float>(triangle.v1().y));
				cached_surface_y = std::max(cached_surface_y, static_cast<float>(triangle.v2().y));
			}
		}
		surface_cached = true;
	}

	// Generate all line instances
	for (const PhotonPath& path : simulator_->paths) {
			if (path.head) {
				// Generate connected line segments with energy gradient
				auto current = path.head;
				auto next = current ? current->next : nullptr;

				// Generate incident ray from source to tissue surface (cached surface)
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
						float t = (cached_surface_y - source_pos.y) / source_dir.y;
						glm::vec3 surface_entry = source_pos + t * source_dir;

						// Add incident ray
						glm::vec4 incident_color(1.0f, 1.0f, 1.0f, 1.0f);
						line_instances_.push_back({source_pos, surface_entry, incident_color});

						// Add refracted ray if needed
						if (first_interaction.y < cached_surface_y - 0.001f) {
							glm::vec4 refracted_color(0.9f, 0.9f, 1.0f, 0.8f);
							line_instances_.push_back({surface_entry, first_interaction, refracted_color});
						}
					}
				}

			while (current && next) {
				// PERFORMANCE OPTIMIZATION: Single line segment instead of 10 gradient segments
				float energy1 = static_cast<float>(current->value);
				float energy2 = static_cast<float>(next->value);
				float avg_energy = (energy1 + energy2) * 0.5f;

				glm::vec4 line_color = adaptive_log_color(avg_energy);
				line_color.a = 1.0f;

				glm::vec3 start(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
								static_cast<float>(current->position.z));
				glm::vec3 end(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
							  static_cast<float>(next->position.z));

				// Single line segment per path segment (10x fewer instances!)
				line_instances_.push_back({start, end, line_color});

				// Move to next segment
				current = next;
				next = current->next;
			}

			// Also generate emitted paths if they exist
			current = path.head;
			while (current) {
				if (current->emit) {
					float emit_energy = static_cast<float>(current->value);
					glm::vec4 emit_color = adaptive_log_color(emit_energy);
					emit_color.a = 1.0f;

					glm::vec3 start(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
									static_cast<float>(current->position.z));
					glm::vec3 emit_end(static_cast<float>(current->emit->position.x),
									   static_cast<float>(current->emit->position.y),
									   static_cast<float>(current->emit->position.z));
					
					line_instances_.push_back({start, emit_end, emit_color});
				}
				current = current->next;
			}
		}
	}

	// Render all line instances (MASSIVE performance improvement)
	if (!line_instances_.empty() && line_instanced_shader_program_) {
		glUseProgram(line_instanced_shader_program_);
		
		// Set MVP matrix
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		GLint mvp_location = glGetUniformLocation(line_instanced_shader_program_, "uMVP");
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));
		
		// Upload instance buffer (use GL_STATIC_DRAW for better performance)
		glBindBuffer(GL_ARRAY_BUFFER, line_instance_vbo_);
		glBufferData(GL_ARRAY_BUFFER, line_instances_.size() * sizeof(LineInstance), 
					 line_instances_.data(), GL_STATIC_DRAW);
		
		// Render
		glBindVertexArray(line_instanced_vao_);
		glDrawArraysInstanced(GL_LINES, 0, 2, static_cast<GLsizei>(line_instances_.size()));
		glBindVertexArray(0);
		
		glUseProgram(0);
	}

	// Collect scatter points for instanced rendering
	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			auto current = path.head;
			int vertex_count = 0;

			// Count total vertices to identify key points properly
			auto temp = current;
			while (temp) {
				vertex_count++;
				temp = temp->next;
			}

			// Only add markers at specific key points: incident, scatter, exit
			auto path_current = path.head;
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
					// Add to point instances instead of traditional rendering
					PointInstance point_instance;
					point_instance.position = pos;
					point_instance.color = marker_color;
					point_instance.size = 8.0f; // Point size in pixels
					point_instances.push_back(point_instance);
				}

				prev = path_current;
				path_current = path_current->next;
				current_index++;
			}
		}
	}

	// Render scatter points using instanced rendering
	if (!point_instances.empty()) {
		glUseProgram(point_instanced_shader_program_);
		
		// Set MVP matrix uniform
		glm::mat4 mvp = camera_.get_projection_matrix() * camera_.get_view_matrix();
		GLint mvp_loc = glGetUniformLocation(point_instanced_shader_program_, "uMVP");
		glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, glm::value_ptr(mvp));

		// Upload instance data
		glBindBuffer(GL_ARRAY_BUFFER, point_instance_vbo_);
		glBufferData(GL_ARRAY_BUFFER, point_instances.size() * sizeof(PointInstance), point_instances.data(), GL_DYNAMIC_DRAW);

		// Bind VAO and render
		glBindVertexArray(point_instanced_vao_);
		
		// Enable point size and blending for nice looking points
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
		glDrawArraysInstanced(GL_POINTS, 0, 1, static_cast<GLsizei>(point_instances.size()));
		
		glDisable(GL_BLEND);
		glDisable(GL_PROGRAM_POINT_SIZE);
		glBindVertexArray(0);
		glUseProgram(0);
	}
}

void Renderer::draw_emitters(const Settings& settings) {
	if (!simulator_ || simulator_->emitters.empty()) {
		return;
	}

	// Draw exit points as bright colored points at the true geometry boundary
	begin_points();
	
	for (const auto& emitter : simulator_->emitters) {
		// Convert position to float
		glm::vec3 exit_pos(
			static_cast<float>(emitter.position.x),
			static_cast<float>(emitter.position.y), 
			static_cast<float>(emitter.position.z)
		);
		
		// Color based on emitted weight - bright colors for high-energy exits
		float weight = static_cast<float>(emitter.weight);
		glm::vec4 exit_color;
		
		if (weight > 0.7f) {
			// High energy: bright white/yellow
			exit_color = glm::vec4(1.0f, 1.0f, 0.8f + 0.2f * weight, 1.0f);
		} else if (weight > 0.4f) {
			// Medium energy: orange to yellow
			float t = (weight - 0.4f) / 0.3f;
			exit_color = glm::vec4(1.0f, 0.6f + 0.4f * t, 0.2f * t, 1.0f);
		} else {
			// Low energy: red
			float t = weight / 0.4f;
			exit_color = glm::vec4(0.8f + 0.2f * t, 0.2f * t, 0.1f * t, 1.0f);
		}
		
		add_point(exit_pos, exit_color);
	}
	
	end_points();
	draw_points();
	
	// Optionally draw exit direction vectors
	if (settings.draw_paths) {  // Use existing draw_paths setting
		begin_lines();
		
		for (const auto& emitter : simulator_->emitters) {
			if (emitter.weight > 0.01) { // Only show significant exits
				glm::vec3 start_pos(
					static_cast<float>(emitter.position.x),
					static_cast<float>(emitter.position.y),
					static_cast<float>(emitter.position.z)
				);
				
				glm::vec3 end_pos(
					static_cast<float>(emitter.position.x + emitter.direction.x * 0.05),
					static_cast<float>(emitter.position.y + emitter.direction.y * 0.05),
					static_cast<float>(emitter.position.z + emitter.direction.z * 0.05)
				);
				
				// Direction line color - dimmer than the exit point
				float weight = static_cast<float>(emitter.weight);
				glm::vec4 direction_color;
				if (weight > 0.7f) {
					direction_color = glm::vec4(0.8f, 0.8f, 0.6f, 0.8f);
				} else {
					direction_color = glm::vec4(0.6f, 0.4f, 0.2f, 0.6f);
				}
				
				add_line(start_pos, end_pos, direction_color);
			}
		}
		
		end_lines();
		draw_lines();
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

	GLint mvp_location = glGetUniformLocation(line_shader_program_, "uMVP");
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

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

	GLint mvp_location = glGetUniformLocation(point_shader_program_, "uMVP");
	GLint size_location = glGetUniformLocation(point_shader_program_, "uPointSize");

	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));
	glUniform1f(size_location, 8.0f); // Smaller spheres as requested

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

	GLint mvp_location = glGetUniformLocation(triangle_shader_program_, "uMVP");
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

	// Set clipping plane uniforms
	GLint num_planes_location = glGetUniformLocation(triangle_shader_program_, "uNumClipPlanes");
	GLint clip_planes_location = glGetUniformLocation(triangle_shader_program_, "uClipPlanes");
	
	int num_planes = std::min(static_cast<int>(clipping_planes.size()), 6);
	glUniform1i(num_planes_location, num_planes);
	
	// Enable OpenGL clipping planes
	for (int i = 0; i < num_planes && i < 6; i++) {
		glEnable(GL_CLIP_DISTANCE0 + i);
	}
	
	if (num_planes > 0 && clip_planes_location != -1) {
		glUniform4fv(clip_planes_location, num_planes, glm::value_ptr(clipping_planes[0]));
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
	bool many_photons = (simulator_->paths.size() > 10);
	
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
	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			auto current = path.head;
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
					// Incidence point - check if at surface (reflected) or inside (transmitted into tissue)
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
						// Inside tissue - this is just the incident energy
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

					// Calculate tissue bottom boundary
					float bottom_y = -0.1f; // Default bottom
					if (!simulator_->layers.empty()) {
						bottom_y = 1000.0f; // Start very high
						for (const auto& layer : simulator_->layers) {
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

					// Calculate tissue top boundary as well
					float top_y = 0.1f;   // Default top
					if (!simulator_->layers.empty()) {
						top_y = -1000.0f; // Start very low
						for (const auto& layer : simulator_->layers) {
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
						// Inside tissue - just show percentage
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
							// Calculate tissue boundaries
							float top_y = 0.1f, bottom_y = -0.1f;
							if (!simulator_->layers.empty()) {
								top_y = -1000.0f;
								bottom_y = 1000.0f;
								for (const auto& layer : simulator_->layers) {
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

							if (at_top) {
								// At incident interface - reflected light
								if (energy_percent < 1.0f) {
									label_text = "Reflected: <1%";
								}
								else {
									label_text = std::format("Reflected: {}%", static_cast<int>(energy_percent));
								}
							}
							else if (at_bottom) {
								// At exit interface - transmitted light
								if (energy_percent < 1.0f) {
									label_text = "Transmitted: <1%";
								}
								else {
									label_text = std::format("Transmitted: {}%", static_cast<int>(energy_percent));
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
	if (simulator_->record.specular_reflection > 0.0) {
		// Calculate incident surface position
		glm::vec3 source_pos(0.05f, 0.0f, 0.05f); // From the incident ray calculations
		glm::vec3 source_dir(0.0f, 1.0f, 0.0f);   // Upward direction

		// Find topmost surface point
		float surface_y = 0.1f;
		if (!simulator_->layers.empty()) {
			surface_y = -1000.0f;
			for (const auto& layer : simulator_->layers) {
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
			float surface_scattering = static_cast<float>(simulator_->record.specular_reflection);
			float scatter_percent = surface_scattering * 100.0f;

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
			surface_label.color = get_layer_specific_energy_color(surface_scattering, 0.0f, 1.0f, 0);
			surface_label.scale = 1.0f;

			cached_energy_labels_.push_back(surface_label);
		}
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
	float current_distance = camera_.get_distance();
	float current_azimuth = camera_.get_azimuth();
	float current_elevation = camera_.get_elevation();

	bool camera_changed = false;
	if (glm::length(current_position - last_camera_position_) > 1e-4f
		|| std::abs(current_distance - last_camera_distance_) > 1e-4f
		|| std::abs(current_azimuth - last_camera_azimuth_) > 1e-4f
		|| std::abs(current_elevation - last_camera_elevation_) > 1e-4f) {
		camera_changed = true;
	}

	// Update screen positions only when camera changes
	if (camera_changed || camera_state_changed_) {
		last_camera_position_ = current_position;
		last_camera_distance_ = current_distance;
		last_camera_azimuth_ = current_azimuth;
		last_camera_elevation_ = current_elevation;
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
	// Non-linear normalization to better visualize energy distribution
	// Most voxels have very low energy, so we need to expand that range visually

	// Clamp energy to valid range
	energy = std::clamp(energy, min_energy, max_energy);

	// First do linear normalization
	float linear_normalized = (energy - min_energy) / (max_energy - min_energy);
	linear_normalized = std::clamp(linear_normalized, 0.0f, 1.0f);

	// Apply power function to expand low-energy visualization
	// Using power of 0.3 significantly expands the low-energy range
	float normalized = std::pow(linear_normalized, 0.3f);

	// Alternative: logarithmic mapping for even more low-energy expansion
	// float normalized = std::log10(linear_normalized * 9.0f + 1.0f); // log10(1) to log10(10) = 0 to 1

	normalized = std::clamp(normalized, 0.0f, 1.0f);

	// Create more distinct color zones with emphasis on lower energies
	// Now the low energy ranges get much more visual space

	if (normalized > 0.85f) {
		// Very high energy: pure white (still top ~5% but harder to reach)
		return glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	}
	else if (normalized > 0.70f) {
		// High energy: white to bright yellow with smooth transition
		float t = (normalized - 0.70f) / 0.15f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 1.0f - t * 0.2f; // From white (1.0) to bright yellow (0.8)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.55f) {
		// Medium-high energy: bright yellow to yellow
		float t = (normalized - 0.55f) / 0.15f;
		float r = 1.0f;
		float g = 1.0f;
		float b = 0.8f - t * 0.3f; // From bright yellow (0.8) to yellow (0.5)
		return glm::vec4(r, g, b, 1.0f);
	}
	else if (normalized > 0.40f) {
		// Medium energy: yellow to orange
		float t = (normalized - 0.40f) / 0.15f;
		float r = 1.0f;
		float g = 1.0f - t * 0.3f; // From 1.0 (yellow) to 0.7 (orange)
		float b = 0.5f - t * 0.5f; // From 0.5 (yellow) to 0.0 (orange)
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

	// Define base colors for different tissue types
	glm::vec3 base_colors[6] = {
		glm::vec3(1.0f, 0.4f, 0.4f), // Red theme for tissue 0
		glm::vec3(0.4f, 0.4f, 1.0f), // Blue theme for tissue 1
		glm::vec3(0.4f, 1.0f, 0.4f), // Green theme for tissue 2
		glm::vec3(1.0f, 0.4f, 1.0f), // Magenta theme for tissue 3
		glm::vec3(1.0f, 1.0f, 0.4f), // Yellow theme for tissue 4
		glm::vec3(0.4f, 1.0f, 1.0f), // Cyan theme for tissue 5
	};

	glm::vec3 base_color = base_colors[tissue_id % 6];

	// Create energy gradient within the tissue's color theme
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
	voxel_instances_.clear();
}

void Renderer::add_voxel_instance(const glm::vec3& position, const glm::vec4& color, float scale, float depth) {
	voxel_instances_.push_back({position, color, scale, depth});
}

void Renderer::end_voxel_instances() {
	if (voxel_instances_.empty()) return;
	
	// Sort by depth (back to front) for proper transparency
	std::sort(voxel_instances_.begin(), voxel_instances_.end(), 
		[](const VoxelInstance& a, const VoxelInstance& b) {
			return a.depth > b.depth; // Back to front
		});

	// Upload instance data to GPU
	glBindBuffer(GL_ARRAY_BUFFER, voxel_instance_vbo_);
	glBufferData(GL_ARRAY_BUFFER,
				 voxel_instances_.size() * sizeof(VoxelInstance), 
				 voxel_instances_.data(), 
				 GL_DYNAMIC_DRAW);
}void Renderer::draw_voxel_instances() {
	if (voxel_instances_.empty() || !voxel_shader_program_) return;
	
	// Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// Disable depth writing for transparent objects but keep depth testing
	glDepthMask(GL_FALSE);
	glEnable(GL_DEPTH_TEST);
	
	glUseProgram(voxel_shader_program_);
	
	// Set MVP matrix
	GLint mvp_location = glGetUniformLocation(voxel_shader_program_, "uMVP");
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));
	
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

	// Check if we need to recalculate (data has changed)
	size_t current_path_count = simulator_->paths.size();
	bool paths_changed = (current_path_count != last_path_count_);
	
	// For now, assume voxel data version changes when simulation runs
	// In a more advanced implementation, you'd track voxel data changes
	bool voxels_changed = !energy_range_cached_;

	if (energy_range_cached_ && !paths_changed && !voxels_changed) {
		return; // Use cached values
	}

	// Recalculate energy range using the ORIGINAL method for accurate color mapping
	std::vector<float> all_energies;
	
	// Collect energies from photon paths (SAME as original)
	for (const PhotonPath& path : simulator_->paths) {
		if (path.head) {
			auto current = path.head;
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
	if (simulator_->config.nx() > 0 && simulator_->config.ny() > 0 && simulator_->config.nz() > 0) {
		int nx = simulator_->config.nx();
		int ny = simulator_->config.ny();
		int nz = simulator_->config.nz();
		
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
							float emittance = static_cast<float>(voxel->emittance);
							
							// Use the current settings to determine energy calculation method
							float total_energy;
							if (settings.voxel_mode == VoxelMode::Absorption) {
								total_energy = absorption;
							}
							else if (settings.voxel_mode == VoxelMode::Emittance) {
								total_energy = emittance;
							}
							else if (settings.voxel_mode == VoxelMode::Layers) {
								total_energy = (voxel->tissue != nullptr) ? 1.0f : 0.0f;
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
	last_path_count_ = current_path_count;
	last_voxel_data_version_++; // Increment version
}

// Invalidate energy cache when simulation data changes
void Renderer::invalidate_energy_cache() const {
	energy_range_cached_ = false;
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
	
	// Instance color (location 3)
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1); // One per instance
	
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
