#include "renderer.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "camera.hpp"
#include "config.hpp"
#include "settings.hpp"
#include "simulator/simulator.hpp"
#include "structs/range3.hpp"

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

	// Check if simulator has changed and invalidate energy label cache if needed
	if (simulator_ != simulator) {
		invalidate_energy_label_cache();
	}

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

	// Update energy label screen positions once per frame (fixes label positioning stability)
	update_energy_label_screen_positions();

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

	// Invalidate energy label cache when viewport changes to ensure proper coordinate conversion
	invalidate_energy_label_cache();
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
	// Mark camera state as changed to update label positions
	camera_state_changed_ = true;
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
	float min_energy = 1e-8f; // Much lower threshold for emittance detection
	float max_energy = 0.01f; // Default fallback

	if (!all_energies.empty()) {
		std::sort(all_energies.begin(), all_energies.end());

		// Use percentile-based scaling for better contrast
		size_t count = all_energies.size();
		float p5 = all_energies[static_cast<size_t>(count * 0.05)];  // 5th percentile
		float p95 = all_energies[static_cast<size_t>(count * 0.95)]; // 95th percentile

		min_energy = std::max(p5, 1e-8f);                            // Lower minimum
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
										surface_y = std::max({surface_y, static_cast<float>(triangle.v0.y),
															  static_cast<float>(triangle.v1.y),
															  static_cast<float>(triangle.v2.y)});
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
										float surface_scattering = static_cast<float>(simulator_->record.rs);
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
						else {
							total_energy = absorption + emittance; // Fallback to combined
						}

						// Show even the tiniest energy interactions
						glm::vec4 color(0.0f, 0.0f, 0.0f, 0.0f);                // Default transparent

						if (total_energy > 1e-8f || voxel->tissue != nullptr) { // Much lower threshold
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
							float min_alpha = 0.2f; // Higher minimum for emittance visibility
							float max_alpha = 0.9f; // Higher maximum visibility

							// Special case for emittance mode: boost visibility significantly
							if (settings.voxel_mode == VoxelMode::Emittance) {
								min_alpha = 0.6f; // Much more visible
								max_alpha = 1.0f; // Full opacity for high emittance
							}

							float alpha =
								min_alpha + (max_alpha - min_alpha) * gamma_corrected * (1.0f - 0.2f * normalized_dist);

							// Use layer-specific energy color mapping with dynamic range
							if (total_energy > 1e-8f) { // Lower threshold
								color = get_layer_specific_energy_color(total_energy, min_energy, max_energy,
																		voxel->tissue->id);
								color.a = alpha;
							}
							else if (voxel->tissue != nullptr) {
								// Voxel in medium but no recorded energy: very faint tissue-colored hint
								color =
									get_layer_specific_energy_color(1e-8f, min_energy, max_energy, voxel->tissue->id);
								color.a = 0.05f; // Slightly more visible
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
							float y0 = static_cast<float>(triangle.v0.y);
							float y1 = static_cast<float>(triangle.v1.y);
							float y2 = static_cast<float>(triangle.v2.y);
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

	// For now, disable tissue interface visualization completely
	// You can enable it later with a UI toggle if desired
	// draw_tissue_interfaces(simulator);
}

void Renderer::draw_tissue_interfaces(Simulator* simulator) {
	if (!simulator || simulator->voxels.empty()) {
		return;
	}

	// Much more subtle approach - only draw a few key boundary indicators
	begin_lines();

	// Very subtle colors - barely visible
	std::vector<glm::vec4> tissue_colors = {
		glm::vec4(0.3f, 0.3f, 0.3f, 0.1f), // Very faint gray for tissue 0
		glm::vec4(0.4f, 0.2f, 0.2f, 0.1f), // Very faint red for tissue 1
		glm::vec4(0.2f, 0.4f, 0.2f, 0.1f), // Very faint green for tissue 2
		glm::vec4(0.2f, 0.2f, 0.4f, 0.1f), // Very faint blue for tissue 3
		glm::vec4(0.4f, 0.4f, 0.2f, 0.1f), // Very faint yellow for tissue 4
		glm::vec4(0.4f, 0.2f, 0.4f, 0.1f), // Very faint magenta for tissue 5
	};

	// Get MCML grid parameters
	int nx = simulator->config.nx;
	int ny = simulator->config.ny;
	int nz = simulator->config.nz;
	double voxsize = simulator->config.vox_size;
	double half_voxsize = voxsize * 0.5;

	// Only draw horizontal lines at Y boundaries where tissues actually change (layer interfaces)
	// Skip most voxels to reduce visual clutter dramatically
	for (int iz = 0; iz < nz; iz += 4) {         // Only every 4th slice
		for (int iy = 0; iy < ny; ++iy) {
			for (int ix = 0; ix < nx; ix += 4) { // Only every 4th voxel
				size_t voxel_index = static_cast<size_t>(iz) * static_cast<size_t>(nx) * static_cast<size_t>(ny)
									 + static_cast<size_t>(iy) * static_cast<size_t>(nx) + static_cast<size_t>(ix);

				if (voxel_index >= simulator->voxels.size())
					continue;

				Voxel* current_voxel = simulator->voxels[voxel_index];
				if (!current_voxel || !current_voxel->tissue)
					continue;

				// Only check Y+1 neighbor for horizontal layer boundaries
				if (iy + 1 < ny) {
					size_t neighbor_index = static_cast<size_t>(iz) * static_cast<size_t>(nx) * static_cast<size_t>(ny)
											+ static_cast<size_t>(iy + 1) * static_cast<size_t>(nx)
											+ static_cast<size_t>(ix);
					if (neighbor_index < simulator->voxels.size()) {
						Voxel* neighbor = simulator->voxels[neighbor_index];
						if (!neighbor || !neighbor->tissue || neighbor->tissue->id != current_voxel->tissue->id) {
							// Found a tissue boundary - draw just a small horizontal line indicator
							double x = simulator->bounds.x_min + (voxsize * ix) + half_voxsize;
							double y = simulator->bounds.y_min + (voxsize * iy) + half_voxsize;
							double z = simulator->bounds.z_min + (voxsize * iz) + half_voxsize;

							float interface_y = static_cast<float>(y + half_voxsize);
							glm::vec4 tissue_color = tissue_colors[current_voxel->tissue->id % tissue_colors.size()];
							tissue_color.a = 0.2f; // Slightly more visible at boundaries

							// Just draw a small cross mark, not full wireframe
							float mark_size = static_cast<float>(voxsize * 0.3);
							glm::vec3 center(static_cast<float>(x), interface_y, static_cast<float>(z));

							add_line(center - glm::vec3(mark_size, 0, 0), center + glm::vec3(mark_size, 0, 0),
									 tissue_color);
							add_line(center - glm::vec3(0, 0, mark_size), center + glm::vec3(0, 0, mark_size),
									 tissue_color);
						}
					}
				}
			}
		}
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

void Renderer::cache_energy_labels() {
	cached_energy_labels_.clear();
	energy_labels_cached_ = false;

	if (!simulator_)
		return;

	// Collect energy labels from photon paths
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
							label_text = "Reflected: " + std::to_string(static_cast<int>(energy_percent)) + "%";
						}
					}
					else {
						// Inside tissue - this is just the incident energy
						if (energy_percent < 1.0f) {
							label_text = "<1%";
						}
						else {
							label_text = std::to_string(static_cast<int>(energy_percent)) + "%";
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
								float y0 = static_cast<float>(triangle.v0.y);
								float y1 = static_cast<float>(triangle.v1.y);
								float y2 = static_cast<float>(triangle.v2.y);
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
								float y0 = static_cast<float>(triangle.v0.y);
								float y1 = static_cast<float>(triangle.v1.y);
								float y2 = static_cast<float>(triangle.v2.y);
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
							label_text = "Reflected: " + std::to_string(static_cast<int>(energy_percent)) + "%";
						}
					}
					else if (at_bottom) {
						// At exit interface - this is transmitted light
						if (energy_percent < 1.0f) {
							label_text = "Transmitted: <1%";
						}
						else {
							label_text = "Transmitted: " + std::to_string(static_cast<int>(energy_percent)) + "%";
						}
					}
					else {
						// Inside tissue - just show percentage
						if (energy_percent < 1.0f) {
							label_text = "<1%";
						}
						else {
							label_text = std::to_string(static_cast<int>(energy_percent)) + "%";
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
										float y0 = static_cast<float>(triangle.v0.y);
										float y1 = static_cast<float>(triangle.v1.y);
										float y2 = static_cast<float>(triangle.v2.y);
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
									label_text = "Reflected: " + std::to_string(static_cast<int>(energy_percent)) + "%";
								}
							}
							else if (at_bottom) {
								// At exit interface - transmitted light
								if (energy_percent < 1.0f) {
									label_text = "Transmitted: <1%";
								}
								else {
									label_text =
										"Transmitted: " + std::to_string(static_cast<int>(energy_percent)) + "%";
								}
							}
							else {
								// Internal scattering - just percentage
								if (energy_percent < 1.0f) {
									label_text = "<1%";
								}
								else {
									label_text = std::to_string(static_cast<int>(energy_percent)) + "%";
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
	if (simulator_->record.rs > 0.0) {
		// Calculate incident surface position
		glm::vec3 source_pos(0.05f, 0.0f, 0.05f); // From the incident ray calculations
		glm::vec3 source_dir(0.0f, 1.0f, 0.0f);   // Upward direction

		// Find topmost surface point
		float surface_y = 0.1f;
		if (!simulator_->layers.empty()) {
			surface_y = -1000.0f;
			for (const auto& layer : simulator_->layers) {
				for (const auto& triangle : layer.mesh) {
					float y0 = static_cast<float>(triangle.v0.y);
					float y1 = static_cast<float>(triangle.v1.y);
					float y2 = static_cast<float>(triangle.v2.y);
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
			float surface_scattering = static_cast<float>(simulator_->record.rs);
			float scatter_percent = surface_scattering * 100.0f;

			std::string label_text;
			if (scatter_percent < 1.0f) {
				label_text = "Reflected: <1%";
			}
			else {
				label_text = "Reflected: " + std::to_string(static_cast<int>(scatter_percent)) + "%";
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

void Renderer::draw_energy_labels(const Settings& settings) {
	if (!simulator_ || !settings.draw_path_labels)
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
