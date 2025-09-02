#include "renderer.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "../simulator/simulator.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

Renderer::Renderer() :
	line_vao_(0), line_vbo_(0), point_vao_(0), point_vbo_(0), triangle_vao_(0), triangle_vbo_(0),
	line_shader_program_(0), point_shader_program_(0), triangle_shader_program_(0),
	is_arc_camera_mode_(true), // Start in Orbit mode
	viewport_width_(800), viewport_height_(600), simulator_(nullptr) {
}

Renderer::~Renderer() {
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

	// Set clear color to dark blue like original renderer (matches the image)
	glClearColor(0.1f, 0.1f, 0.3f, 1.0f);

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

void Renderer::render(Simulator* simulator) {
	if (!simulator)
		return;

	// Store simulator reference for use in drawing functions
	simulator_ = simulator;

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

	// Draw everything based on current settings (user request 4 - re-enable ImGui menu)
	if (settings_.draw_bounds) {
		draw_bounds_modern(simulator);
	}

	// Always draw voxels with current energy-based coloring (keep as default)
	draw_voxels_modern(settings_);

	if (settings_.draw_paths) {
		draw_paths_modern(settings_); // Includes incident photon and scatter markers
	}

	// Always draw current photon positions and directions
	draw_photons_modern(settings_);

	// Draw energy labels as billboards (user request 3)
	draw_energy_labels(settings_);

	// Draw layers (tissue boundaries)
	draw_layers_modern(simulator);
}

void Renderer::draw_coordinate_axes() {
	// Draw coordinate axes with muted colors to not interfere with data
	begin_lines();

	// X axis - Muted Red
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec4(0.8f, 0.3f, 0.3f, 0.6f));

	// Y axis - Muted Green
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec4(0.3f, 0.8f, 0.3f, 0.6f));

	// Z axis - Muted Blue
	add_line(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f), glm::vec4(0.3f, 0.3f, 0.8f, 0.6f));

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

void Renderer::draw_bounds_modern(Simulator* simulator) {
	if (!simulator)
		return;

	begin_lines();

	// Get bounds
	Range3 bounds = simulator->bounds;
	float minX = static_cast<float>(bounds.x_min);
	float maxX = static_cast<float>(bounds.x_max);
	float minY = static_cast<float>(bounds.y_min);
	float maxY = static_cast<float>(bounds.y_max);
	float minZ = static_cast<float>(bounds.z_min);
	float maxZ = static_cast<float>(bounds.z_max);

	glm::vec4 bounds_color(0.7f, 0.7f, 0.7f, 0.6f); // Subtle gray bounds like backup

	// Bottom face edges
	add_line(glm::vec3(minX, minY, minZ), glm::vec3(maxX, minY, minZ), bounds_color);
	add_line(glm::vec3(maxX, minY, minZ), glm::vec3(maxX, maxY, minZ), bounds_color);
	add_line(glm::vec3(maxX, maxY, minZ), glm::vec3(minX, maxY, minZ), bounds_color);
	add_line(glm::vec3(minX, maxY, minZ), glm::vec3(minX, minY, minZ), bounds_color);

	// Top face edges
	add_line(glm::vec3(minX, minY, maxZ), glm::vec3(maxX, minY, maxZ), bounds_color);
	add_line(glm::vec3(maxX, minY, maxZ), glm::vec3(maxX, maxY, maxZ), bounds_color);
	add_line(glm::vec3(maxX, maxY, maxZ), glm::vec3(minX, maxY, maxZ), bounds_color);
	add_line(glm::vec3(minX, maxY, maxZ), glm::vec3(minX, minY, maxZ), bounds_color);

	// Vertical edges
	add_line(glm::vec3(minX, minY, minZ), glm::vec3(minX, minY, maxZ), bounds_color);
	add_line(glm::vec3(maxX, minY, minZ), glm::vec3(maxX, minY, maxZ), bounds_color);
	add_line(glm::vec3(maxX, maxY, minZ), glm::vec3(maxX, maxY, maxZ), bounds_color);
	add_line(glm::vec3(minX, maxY, minZ), glm::vec3(minX, maxY, maxZ), bounds_color);

	end_lines();
	draw_lines();
}

void Renderer::update_camera() {
	// Update camera matrices for viewport
	float aspect = static_cast<float>(viewport_width_) / static_cast<float>(viewport_height_);
	camera_.set_aspect_ratio(aspect);
}

void Renderer::update_camera_target(Simulator* simulator) {
	if (!simulator)
		return;

	Range3 bounds = simulator->bounds;

	// Set camera target to center of bounds
	glm::vec3 center(static_cast<float>((bounds.x_min + bounds.x_max) * 0.5),
					 static_cast<float>((bounds.y_min + bounds.y_max) * 0.5),
					 static_cast<float>((bounds.z_min + bounds.z_max) * 0.5));
	camera_.set_target(center);

	// Set tissue elevation bounds for scroll wheel constraint (Y-axis)
	camera_.set_elevation_bounds(static_cast<float>(bounds.y_min), static_cast<float>(bounds.y_max));
}

void Renderer::set_viewport(int width, int height) {
	viewport_width_ = width;
	viewport_height_ = height;
	glViewport(0, 0, width, height);
}

// Input handlers (simplified for now)
void Renderer::handle_key_input(int key, int scancode, int action, int mods) {
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
}

void Renderer::handle_mouse_button(int button, int action, int mods) {
	camera_.handle_mouse_button(button, action);
}

void Renderer::handle_mouse_scroll(float xoffset, float yoffset) {
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

void Renderer::draw_voxels_modern(const Settings& settings) {
	if (!simulator_) {
		return;
	}

	// Only draw voxels if voxel mode is not None (user request 4)
	if (settings.voxel_mode == VoxelMode::None) {
		return;
	}

	// Draw actual voxels using MCML's EXACT grid system like backup
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_FALSE);   // Don't write to depth buffer for transparent voxels
	glDisable(GL_CULL_FACE); // Ensure all faces render regardless of angle (user request 2)

	begin_triangles();

	// Get MCML grid parameters - EXACTLY like backup
	int nx = simulator_->config.nx;
	int ny = simulator_->config.ny;
	int nz = simulator_->config.nz;
	double voxsize = simulator_->config.vox_size;
	double half_voxsize = voxsize * 0.5;

	// Collect voxels with distance for depth sorting (user request 4)
	struct VoxelRenderData {
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
				int voxel_index = iz * nx * ny + iy * nx + ix;

				if (voxel_index >= 0 && static_cast<size_t>(voxel_index) < simulator_->voxels.size()) {
					Voxel* voxel = simulator_->voxels[static_cast<size_t>(voxel_index)];

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
	float min_energy = 0.0000001f; // Minimum threshold for visibility
	float max_energy = 0.01f;      // Default fallback

	if (!all_energies.empty()) {
		std::sort(all_energies.begin(), all_energies.end());

		// Use percentile-based scaling for better contrast
		size_t count = all_energies.size();
		float p5 = all_energies[static_cast<size_t>(count * 0.05)];  // 5th percentile
		float p95 = all_energies[static_cast<size_t>(count * 0.95)]; // 95th percentile

		min_energy = std::max(p5, 0.0000001f);
		max_energy = std::max(p95, min_energy * 10.0f);              // Ensure reasonable range
	}

	// Second pass: collect all voxels with their render data using adaptive scaling
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				// Calculate the voxel index using MCML's exact formula
				int voxel_index = iz * nx * ny + iy * nx + ix;

				if (voxel_index >= 0 && static_cast<size_t>(voxel_index) < simulator_->voxels.size()) {
					Voxel* voxel = simulator_->voxels[static_cast<size_t>(voxel_index)];

					if (voxel) {
						// Calculate center position using MCML's exact formula
						double x = simulator_->bounds.x_min + (voxsize * ix) + half_voxsize;
						double y = simulator_->bounds.y_min + (voxsize * iy) + half_voxsize;
						double z = simulator_->bounds.z_min + (voxsize * iz) + half_voxsize;

						float fx = static_cast<float>(x);
						float fy = static_cast<float>(y);
						float fz = static_cast<float>(z);
						glm::vec3 voxel_pos(fx, fy, fz);

						// Energy-based coloring based on selected voxel mode (user request 4)
						float absorption = static_cast<float>(voxel->absorption);
						float emittance = static_cast<float>(voxel->emittance);

						float total_energy;
						if (settings.voxel_mode == VoxelMode::Absorption) {
							total_energy = absorption;
						}
						else if (settings.voxel_mode == VoxelMode::Emittance) {
							total_energy = emittance;
						}
						else {
							total_energy = absorption + emittance; // Fallback to combined
						}

						// Show even the tiniest energy interactions
						glm::vec4 color(0.0f, 0.0f, 0.0f, 0.0f); // Default transparent

						if (total_energy > 0.0000001 || voxel->tissue != nullptr) {
							// Calculate percentage of total energy for this voxel
							float energy_percentage = (total_accumulated_energy > 0.0f)
														  ? (total_energy / total_accumulated_energy) * 100.0f
														  : 0.0f;

							// Use gamma correction for better mid-range contrast
							float normalized_energy = (total_energy - min_energy) / (max_energy - min_energy);
							normalized_energy = std::clamp(normalized_energy, 0.0f, 1.0f);

							// Apply gamma correction for better contrast (gamma = 0.5 brightens mid-tones)
							float gamma_corrected = std::pow(normalized_energy, 0.5f);

							// Calculate distance from origin for depth-based alpha
							float max_dist = glm::length(glm::vec3(static_cast<float>(simulator_->bounds.x_max),
																   static_cast<float>(simulator_->bounds.y_max),
																   static_cast<float>(simulator_->bounds.z_max)));
							float dist_from_origin = glm::length(voxel_pos);
							float normalized_dist = glm::clamp(dist_from_origin / max_dist, 0.0f, 1.0f);

							// More visible alpha scaling based on gamma-corrected energy
							float min_alpha = 0.15f; // More visible base
							float max_alpha = 0.7f;  // Higher maximum visibility
							float alpha =
								min_alpha + (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.3f * normalized_dist);

							// Use adaptive energy color mapping with dynamic range
							if (total_energy > 0.0000001f) {
								color = get_adaptive_energy_color(total_energy, min_energy, max_energy);
								color.a = alpha;
							}
							else if (voxel->tissue != nullptr) {
								// Voxel in medium but no recorded energy: very faint blue
								color = glm::vec4(0.4f, 0.4f, 0.8f, 0.02f);
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

	// Sort voxels back-to-front for proper transparency blending (user request 4)
	std::sort(voxels_to_render.begin(), voxels_to_render.end(), [](const VoxelRenderData& a, const VoxelRenderData& b) {
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

			// Define the 8 corners of the voxel cube exactly like backup
			glm::vec3 vertices[8] = {
				glm::vec3(pos.x - half, pos.y - half, pos.z - half), // 0
				glm::vec3(pos.x + half, pos.y - half, pos.z - half), // 1
				glm::vec3(pos.x + half, pos.y + half, pos.z - half), // 2
				glm::vec3(pos.x - half, pos.y + half, pos.z - half), // 3
				glm::vec3(pos.x - half, pos.y - half, pos.z + half), // 4
				glm::vec3(pos.x + half, pos.y - half, pos.z + half), // 5
				glm::vec3(pos.x + half, pos.y + half, pos.z + half), // 6
				glm::vec3(pos.x - half, pos.y + half, pos.z + half)  // 7
			};

			// Draw all 6 quad faces using triangles to emulate GL_QUADS exactly like backup
			// Front face (z+)
			add_triangle(vertices[4], vertices[5], vertices[6], voxel_data.color);
			add_triangle(vertices[4], vertices[6], vertices[7], voxel_data.color);

			// Back face (z-)
			add_triangle(vertices[1], vertices[0], vertices[3], voxel_data.color);
			add_triangle(vertices[1], vertices[3], vertices[2], voxel_data.color);

			// Left face (x-)
			add_triangle(vertices[0], vertices[4], vertices[7], voxel_data.color);
			add_triangle(vertices[0], vertices[7], vertices[3], voxel_data.color);

			// Right face (x+)
			add_triangle(vertices[5], vertices[1], vertices[2], voxel_data.color);
			add_triangle(vertices[5], vertices[2], vertices[6], voxel_data.color);

			// Top face (y+)
			add_triangle(vertices[3], vertices[7], vertices[6], voxel_data.color);
			add_triangle(vertices[3], vertices[6], vertices[2], voxel_data.color);

			// Bottom face (y-)
			add_triangle(vertices[0], vertices[1], vertices[5], voxel_data.color);
			add_triangle(vertices[0], vertices[5], vertices[4], voxel_data.color);
		}

		// Upload vertex data to GPU
		end_triangles();

		// Render current batch
		draw_triangles();
	}

	// Restore OpenGL state
	glDepthMask(GL_TRUE);   // Re-enable depth writing
	glEnable(GL_CULL_FACE); // Re-enable face culling
}

void Renderer::draw_paths_modern(const Settings& settings) {
	if (!simulator_)
		return;

	// First pass: analyze energy distribution for adaptive logarithmic mapping
	std::vector<float> all_energies;
	for (const Graph& path : simulator_->paths) {
		if (path.head) {
			Vertex* current = path.head;
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
		auto minmax = std::minmax_element(all_energies.begin(), all_energies.end());
		min_energy = *minmax.first;
		max_energy = *minmax.second;
	}

	// Create adaptive logarithmic mapping function using helper method
	auto adaptive_log_color = [this, min_energy, max_energy](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, min_energy, max_energy);
	};

	// Draw photon path histories with adaptive energy-based coloring
	begin_lines();

	for (const Graph& path : simulator_->paths) {
		if (path.head) {
			// Draw connected line segments with energy gradient exactly like backup
			Vertex* current = path.head;
			Vertex* next = current ? current->next : nullptr;

			// First, draw incident photon above the medium (user request 1)
			if (current) {
				glm::vec3 first_point(static_cast<float>(current->position.x), static_cast<float>(current->position.y),
									  static_cast<float>(current->position.z));

				// Get the actual source direction from the simulator
				glm::vec3 source_direction(0.0f, -1.0f, -0.5f); // From config: direction = 0, -1, -0.5
				source_direction = normalize(source_direction);

				// Draw incident photon coming from the source direction
				// Start point should be BEFORE the incidence point, opposite to the source direction
				glm::vec3 incident_start = first_point - source_direction * 0.5f;
				glm::vec4 incident_color(1.0f, 1.0f, 1.0f, 1.0f); // Bright white for incident photon
				add_line(incident_start, first_point, incident_color);

				// TEMPORARILY DISABLED - If there's reflection at the surface, draw the reflected ray going back up
				if (false && next) {
					glm::vec3 next_point(static_cast<float>(next->position.x), static_cast<float>(next->position.y),
										 static_cast<float>(next->position.z));

					// If we're at the surface (z ≈ 0), add a reflected ray going upward
					if (abs(first_point.z) < 0.1f) { // Near surface
						// FORCE reflected ray to point clearly outward - no complex calculation
						glm::vec3 reflected_dir(0.0f, -1.0f, 0.8f); // Same Y as incident, but positive Z (outward)
						reflected_dir = normalize(reflected_dir);

						glm::vec3 reflected_end = first_point + reflected_dir * 0.5f;
						glm::vec4 reflected_color(1.0f, 0.0f, 0.0f, 1.0f); // BRIGHT RED to make it obvious
						add_line(first_point, reflected_end, reflected_color);
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

				// Special handling for first segment (from incident point)
				if (current == path.head && abs(start.z) < 0.1f) {
					// This is the first segment from the surface incident point
					// Instead of drawing the transmitted ray into the medium,
					// draw a scattered ray going outward (this is what you want to see!)
					glm::vec3 scattered_direction = normalize(end - start); // Original direction into medium
					scattered_direction.z = -scattered_direction.z;         // Flip Z to point outward
					glm::vec3 scattered_end = start + scattered_direction * 0.5f;

					// Use BRIGHT GREEN color to identify this ray clearly
					glm::vec4 scattered_color(0.0f, 1.0f, 0.0f, 1.0f); // BRIGHT GREEN
					add_line(start, scattered_end, scattered_color);

					// Skip the original transmitted segment - don't draw it
				}
				else {
					// Normal path segment - create gradient by adding multiple line segments with interpolated colors
					const int gradient_segments = 10;
					for (int i = 0; i < gradient_segments; i++) {
						float t1 = static_cast<float>(i) / static_cast<float>(gradient_segments);
						float t2 = static_cast<float>(i + 1) / static_cast<float>(gradient_segments);

						glm::vec3 seg_start = start + t1 * (end - start);
						glm::vec3 seg_end = start + t2 * (end - start);

						glm::vec4 seg_color = start_color * (1.0f - (t1 + t2) * 0.5f) + end_color * ((t1 + t2) * 0.5f);
						seg_color.a = 1.0f;

						add_line(seg_start, seg_end, seg_color);
					}
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

					// Use the original direction - trust the simulator's calculation
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

	for (const Graph& path : simulator_->paths) {
		if (path.head) {
			Vertex* current = path.head;
			int vertex_count = 0;

			// Count total vertices to identify key points properly
			Vertex* temp = current;
			while (temp) {
				vertex_count++;
				temp = temp->next;
			}

			// Only add markers at specific key points: incident, scatter, exit (user request 2)
			Vertex* path_current = path.head;
			Vertex* prev = nullptr;
			Vertex* next = nullptr;
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

void Renderer::draw_photons_modern(const Settings& settings) {
	if (!simulator_)
		return;

	// Draw current photon positions as energy-colored points exactly like backup
	begin_points();
	for (const auto& photon : simulator_->photons) {
		if (photon.alive && photon.weight > 0.001) {
			// Energy-based coloring: White = full energy, Red = medium energy, Dark red = low energy
			float weight = static_cast<float>(photon.weight);
			float energy = std::max(0.1f, weight); // Clamp to avoid invisible photons

			glm::vec4 photon_color;
			if (energy > 0.7f) {
				// High energy: bright white-yellow
				photon_color = glm::vec4(1.0f, 1.0f, 0.8f + 0.2f * energy, 1.0f);
			}
			else if (energy > 0.4f) {
				// Medium energy: orange-red gradient
				float t = (energy - 0.4f) / 0.3f;
				photon_color = glm::vec4(1.0f, 0.5f + 0.5f * t, 0.2f * t, 1.0f);
			}
			else {
				// Low energy: red to dark red
				float t = energy / 0.4f;
				photon_color = glm::vec4(0.5f + 0.5f * t, 0.1f * t, 0.1f * t, 1.0f);
			}

			glm::vec3 photon_pos(static_cast<float>(photon.position.x), static_cast<float>(photon.position.y),
								 static_cast<float>(photon.position.z));
			add_point(photon_pos, photon_color);
		}
	}
	end_points();
	draw_points();

	// Draw photon direction indicators exactly like backup
	begin_lines();
	for (const auto& photon : simulator_->photons) {
		if (photon.alive && photon.weight > 0.001) {
			// Same energy-based color as above
			float weight = static_cast<float>(photon.weight);
			float energy = std::max(0.1f, weight);

			glm::vec4 direction_color;
			if (energy > 0.7f) {
				direction_color = glm::vec4(1.0f, 1.0f, 0.8f + 0.2f * energy, 1.0f);
			}
			else if (energy > 0.4f) {
				float t = (energy - 0.4f) / 0.3f;
				direction_color = glm::vec4(1.0f, 0.5f + 0.5f * t, 0.2f * t, 1.0f);
			}
			else {
				float t = energy / 0.4f;
				direction_color = glm::vec4(0.5f + 0.5f * t, 0.1f * t, 0.1f * t, 1.0f);
			}

			// Draw direction vector exactly like backup
			glm::vec3 start_pos(static_cast<float>(photon.position.x), static_cast<float>(photon.position.y),
								static_cast<float>(photon.position.z));
			glm::vec3 end_pos(
				static_cast<float>(photon.position.x + photon.direction.x * 0.08), // Direction indicator like backup
				static_cast<float>(photon.position.y + photon.direction.y * 0.08),
				static_cast<float>(photon.position.z + photon.direction.z * 0.08));

			add_line(start_pos, end_pos, direction_color);
		}
	}
	end_lines();
	draw_lines();
}

void Renderer::draw_layers_modern(Simulator* simulator) {
	if (!simulator)
		return;

	begin_lines();

	// Draw tissue layer boundaries as very subtle wireframe outlines exactly like backup
	for (size_t i = 0; i < simulator->layers.size(); ++i) {
		// Very subtle gray wireframe for layer boundaries - EXACT backup colors
		glm::vec4 layer_color(0.4f, 0.4f, 0.4f, 1.0f); // Neutral gray like backup

		// Draw layer as wireframe outline (assuming layers are horizontal) - LIKE BACKUP
		Range3 bounds = simulator->bounds;
		float layerZ = static_cast<float>(
			bounds.z_min
			+ (bounds.z_max - bounds.z_min) * (static_cast<float>(i) / static_cast<float>(simulator->layers.size())));

		// Draw rectangular wireframe using line loop equivalent (4 lines forming a rectangle)
		glm::vec3 corner1(static_cast<float>(bounds.x_min), static_cast<float>(bounds.y_min), layerZ);
		glm::vec3 corner2(static_cast<float>(bounds.x_max), static_cast<float>(bounds.y_min), layerZ);
		glm::vec3 corner3(static_cast<float>(bounds.x_max), static_cast<float>(bounds.y_max), layerZ);
		glm::vec3 corner4(static_cast<float>(bounds.x_min), static_cast<float>(bounds.y_max), layerZ);

		// Draw the 4 edges of the rectangle (emulating GL_LINE_LOOP exactly like backup)
		add_line(corner1, corner2, layer_color);
		add_line(corner2, corner3, layer_color);
		add_line(corner3, corner4, layer_color);
		add_line(corner4, corner1, layer_color);
	}

	end_lines();
	draw_lines();
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
	if (triangle_vertices_.empty() || !triangle_shader_program_)
		return;

	glUseProgram(triangle_shader_program_);

	GLint mvp_location = glGetUniformLocation(triangle_shader_program_, "uMVP");
	glm::mat4 mvp = camera_.get_mvp_matrix();
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

	glBindVertexArray(triangle_vao_);
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(triangle_vertices_.size()));
	glBindVertexArray(0);
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
		char info_log[512];
		glGetProgramInfoLog(program, 512, nullptr, info_log);
		std::cerr << "Program linking failed: " << info_log << std::endl;
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
		char info_log[512];
		glGetShaderInfoLog(shader, 512, nullptr, info_log);
		std::cerr << "Shader compilation failed: " << info_log << std::endl;
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

void Renderer::draw_energy_labels(const Settings& settings) {
	if (!simulator_)
		return;

	// Collect energy labels from photon paths
	std::vector<EnergyLabel> labels;

	for (const Graph& path : simulator_->paths) {
		if (path.head) {
			Vertex* current = path.head;
			int vertex_count = 0;

			// Count vertices
			Vertex* temp = current;
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
					// Incidence point
					should_label = true;
					float energy_percent = static_cast<float>(current->value * 100.0);
					label_text = std::to_string(static_cast<int>(energy_percent)) + "%";
					label_color = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f); // White for incident
				}
				else if (current_index == vertex_count - 1) {
					// End point
					should_label = true;
					float energy_percent = static_cast<float>(current->value * 100.0);
					label_text = std::to_string(static_cast<int>(energy_percent)) + "%";

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
							label_text = std::to_string(static_cast<int>(energy_percent)) + "%";

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

					// Scale based on camera distance
					label.scale = camera_.get_distance() * 0.02f;

					labels.push_back(label);
				}

				current = current->next;
				current_index++;
			}
		}
	}

	// Render labels as billboards using point sprites with offset text rendering
	// For now, we'll render them as simple colored points at a larger size
	begin_points();

	for (const auto& label : labels) {
		// Offset the label position slightly above the actual point to avoid overlapping
		glm::vec3 label_pos = label.world_position + glm::vec3(0.0f, 0.0f, 0.05f);
		add_point(label_pos, label.color);
	}

	end_points();
	draw_points();

	// TODO: Add actual text rendering using texture atlases or bitmap fonts
	// For now, the colored points indicate where energy labels would appear
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
