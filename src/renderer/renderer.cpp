#include "renderer.hpp"

#include <windows.h>

#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "structs/graph.hpp"
#include "structs/layer.hpp"
#include "structs/range3.hpp"
#include "structs/vertex.hpp"
#include "utilities/utilities.hpp"
#include "window.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

// Static member definitions - Simple rendering
std::unique_ptr<GLWindow> Renderer::gl_window_ = nullptr;
glm::mat4 Renderer::projection_matrix_ = glm::mat4(1.0f);
glm::mat4 Renderer::view_matrix_ = glm::mat4(1.0f);
glm::mat4 Renderer::model_matrix_ = glm::mat4(1.0f);
glm::mat4 Renderer::mvp_matrix_ = glm::mat4(1.0f);

// Interactive camera controls
float Renderer::camera_distance_ = 5.0f;
float Renderer::camera_rotation_x_ = 30.0f;
float Renderer::camera_rotation_y_ = 45.0f;

Mouse Renderer::mouse_ = Mouse();
Camera Renderer::camera_ = Camera();
WindowSize Renderer::window_size_ = WindowSize();
Settings Renderer::settings_ = Settings();
Simulator* Renderer::simulator_ = nullptr;

// Simple OpenGL initialization
void Renderer::initialize() {
	// Initialize OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	// Enable blending for transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Set clear color
	glClearColor(0.1f, 0.1f, 0.2f, 1.0f);

	// Set up default perspective
	setup_perspective(1200, 800, 45.0f, 0.1f, 1000.0f);

	std::cout << "Simple OpenGL renderer initialized successfully!" << std::endl;
}

void Renderer::setup_perspective(int width, int height, float fov, float near_plane, float far_plane) {
	float aspect = static_cast<float>(width) / static_cast<float>(height);
	projection_matrix_ = glm::perspective(glm::radians(fov), aspect, near_plane, far_plane);
	update_mvp();

	// Also set legacy OpenGL matrices
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projection_matrix_));
}

void Renderer::set_view_matrix(const glm::vec3& eye, const glm::vec3& center, const glm::vec3& up) {
	view_matrix_ = glm::lookAt(eye, center, up);
	update_mvp();

	// Also set legacy OpenGL matrices
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(view_matrix_));
}

void Renderer::update_mvp() {
	mvp_matrix_ = projection_matrix_ * view_matrix_ * model_matrix_;
}

void Renderer::draw_bounds() {
	if (!simulator_) {
		return;
	}

	// Draw simple bounding box using legacy OpenGL
	glColor3f(1.0f, 1.0f, 1.0f); // White wireframe
	glBegin(GL_LINES);

	Range3 bounds = simulator_->bounds;
	float minX = static_cast<float>(bounds.x_min);
	float maxX = static_cast<float>(bounds.x_max);
	float minY = static_cast<float>(bounds.y_min);
	float maxY = static_cast<float>(bounds.y_max);
	float minZ = static_cast<float>(bounds.z_min);
	float maxZ = static_cast<float>(bounds.z_max);

	// Bottom face
	glVertex3f(minX, minY, minZ);
	glVertex3f(maxX, minY, minZ);
	glVertex3f(maxX, minY, minZ);
	glVertex3f(maxX, maxY, minZ);
	glVertex3f(maxX, maxY, minZ);
	glVertex3f(minX, maxY, minZ);
	glVertex3f(minX, maxY, minZ);
	glVertex3f(minX, minY, minZ);

	// Top face
	glVertex3f(minX, minY, maxZ);
	glVertex3f(maxX, minY, maxZ);
	glVertex3f(maxX, minY, maxZ);
	glVertex3f(maxX, maxY, maxZ);
	glVertex3f(maxX, maxY, maxZ);
	glVertex3f(minX, maxY, maxZ);
	glVertex3f(minX, maxY, maxZ);
	glVertex3f(minX, minY, maxZ);

	// Vertical edges
	glVertex3f(minX, minY, minZ);
	glVertex3f(minX, minY, maxZ);
	glVertex3f(maxX, minY, minZ);
	glVertex3f(maxX, minY, maxZ);
	glVertex3f(maxX, maxY, minZ);
	glVertex3f(maxX, maxY, maxZ);
	glVertex3f(minX, maxY, minZ);
	glVertex3f(minX, maxY, maxZ);

	glEnd();
}

void Renderer::draw_voxels() {
	if (!simulator_) {
		return;
	}

	// Draw actual voxels using MCML's EXACT grid system
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_FALSE); // Don't write to depth buffer for transparent voxels

	// Get MCML grid parameters
	int32_t nx = simulator_->config.nx;
	int32_t ny = simulator_->config.ny;
	int32_t nz = simulator_->config.nz;
	double vox_size = simulator_->config.vox_size;
	double half_voxsize = vox_size * 0.5;

	// Draw ALL voxels in the MCML grid, not just ones with recorded energy
	for (int32_t iz = 0; iz < nz; iz++) {
		for (int32_t iy = 0; iy < ny; iy++) {
			for (int32_t ix = 0; ix < nx; ix++) {
				// Calculate the voxel index using MCML's exact formula
				int32_t voxel_index = iz * nx * ny + iy * nx + ix;

				if (voxel_index >= 0 && static_cast<size_t>(voxel_index) < simulator_->voxels.size()) {
					Voxel* voxel = simulator_->voxels[static_cast<size_t>(voxel_index)];

					if (voxel) {
						// Calculate center position using MCML's exact formula
						double x = simulator_->bounds.x_min + (vox_size * ix) + half_voxsize;
						double y = simulator_->bounds.y_min + (vox_size * iy) + half_voxsize;
						double z = simulator_->bounds.z_min + (vox_size * iz) + half_voxsize;

						float fx = static_cast<float>(x);
						float fy = static_cast<float>(y);
						float fz = static_cast<float>(z);

						// Energy-based coloring - show EVERY voxel, even with tiny energy
						float absorption = static_cast<float>(voxel->absorption);
						float emittance = static_cast<float>(voxel->emittance);
						float total_energy = absorption + emittance;

						// Show even the tiniest energy interactions
						if (total_energy > 0.0000001 || voxel->tissue != nullptr) {
							float normalized = std::min(total_energy / 0.01f, 1.0f); // Very sensitive scale

							// Guaranteed minimum visibility for ALL voxels in the medium
							float min_alpha = 0.05f; // Base visibility
							float max_alpha = 0.5f;  // Maximum visibility
							float alpha = min_alpha + (max_alpha - min_alpha) * normalized;

							// Heat map with enhanced visibility for low energies
							if (total_energy > 0.01f) {
								// High energy: red
								glColor4f(1.0f, 0.2f, 0.2f, alpha);
							}
							else if (total_energy > 0.001f) {
								// Medium energy: yellow-orange
								glColor4f(1.0f, 0.8f, 0.2f, alpha);
							}
							else if (total_energy > 0.0001f) {
								// Low energy: green
								glColor4f(0.2f, 1.0f, 0.4f, alpha);
							}
							else if (total_energy > 0.0000001f) {
								// Very low energy: cyan (this should catch photon paths!)
								glColor4f(0.2f, 0.8f, 1.0f, alpha);
							}
							else if (voxel->tissue != nullptr) {
								// Voxel in medium but no recorded energy: very faint blue
								glColor4f(0.4f, 0.4f, 0.8f, 0.02f);
							}

							// Draw voxel cube - nearly full size for better coverage
							float half = static_cast<float>(half_voxsize * 0.95);

							glBegin(GL_QUADS);

							// Front face (z+)
							glVertex3f(fx - half, fy - half, fz + half);
							glVertex3f(fx + half, fy - half, fz + half);
							glVertex3f(fx + half, fy + half, fz + half);
							glVertex3f(fx - half, fy + half, fz + half);

							// Back face (z-)
							glVertex3f(fx - half, fy - half, fz - half);
							glVertex3f(fx - half, fy + half, fz - half);
							glVertex3f(fx + half, fy + half, fz - half);
							glVertex3f(fx + half, fy - half, fz - half);

							// Left face (x-)
							glVertex3f(fx - half, fy - half, fz - half);
							glVertex3f(fx - half, fy - half, fz + half);
							glVertex3f(fx - half, fy + half, fz + half);
							glVertex3f(fx - half, fy + half, fz - half);

							// Right face (x+)
							glVertex3f(fx + half, fy - half, fz - half);
							glVertex3f(fx + half, fy + half, fz - half);
							glVertex3f(fx + half, fy + half, fz + half);
							glVertex3f(fx + half, fy - half, fz + half);

							// Bottom face (y-)
							glVertex3f(fx - half, fy - half, fz - half);
							glVertex3f(fx + half, fy - half, fz - half);
							glVertex3f(fx + half, fy - half, fz + half);
							glVertex3f(fx - half, fy - half, fz + half);

							// Top face (y+)
							glVertex3f(fx - half, fy + half, fz - half);
							glVertex3f(fx - half, fy + half, fz + half);
							glVertex3f(fx + half, fy + half, fz + half);
							glVertex3f(fx + half, fy + half, fz - half);

							glEnd();
						}
					}
				}
			}
		}
	}

	glDepthMask(GL_TRUE); // Re-enable depth writing
	glDisable(GL_BLEND);
}

void Renderer::draw_layers() {
	if (!simulator_) {
		return;
	}

	// Draw tissue layer boundaries as very subtle wireframe outlines (no more green tint!)
	glDisable(GL_BLEND); // No transparency to avoid visual clutter
	glLineWidth(1.0f);

	for (size_t i = 0; i < simulator_->layers.size(); ++i) {
		// Very subtle gray wireframe for layer boundaries
		glColor3f(0.4f, 0.4f, 0.4f); // Neutral gray

		// Draw layer as wireframe outline (assuming layers are horizontal)
		Range3 bounds = simulator_->bounds;
		float layerZ = static_cast<float>(
			bounds.z_min
			+ (bounds.z_max - bounds.z_min) * (static_cast<float>(i) / static_cast<float>(simulator_->layers.size())));

		glBegin(GL_LINE_LOOP); // Wireframe outline instead of filled quad
		glVertex3f(static_cast<GLfloat>(bounds.x_min), static_cast<GLfloat>(bounds.y_min), layerZ);
		glVertex3f(static_cast<GLfloat>(bounds.x_max), static_cast<GLfloat>(bounds.y_min), layerZ);
		glVertex3f(static_cast<GLfloat>(bounds.x_max), static_cast<GLfloat>(bounds.y_max), layerZ);
		glVertex3f(static_cast<GLfloat>(bounds.x_min), static_cast<GLfloat>(bounds.y_max), layerZ);
		glEnd();
	}
}

void Renderer::draw_photons() {
	if (!simulator_) {
		return;
	}

	// Draw photon paths showing energy loss through color changes
	glLineWidth(3.0f);

	// First, draw current photon positions as energy-colored points
	glPointSize(8.0f);
	glBegin(GL_POINTS);
	for (const auto& photon : simulator_->photons) {
		if (photon.alive && photon.weight > 0.001) {
			// Energy-based coloring: White = full energy, Red = medium energy, Dark red = low energy
			float weight = static_cast<float>(photon.weight);
			float energy = std::max(0.1f, weight); // Clamp to avoid invisible photons

			if (energy > 0.7f) {
				// High energy: bright white-yellow
				glColor3f(1.0f, 1.0f, 0.8f + 0.2f * energy);
			}
			else if (energy > 0.4f) {
				// Medium energy: orange-red gradient
				float t = (energy - 0.4f) / 0.3f;
				glColor3f(1.0f, 0.5f + 0.5f * t, 0.2f * t);
			}
			else {
				// Low energy: red to dark red
				float t = energy / 0.4f;
				glColor3f(0.5f + 0.5f * t, 0.1f * t, 0.1f * t);
			}

			glVertex3f(static_cast<GLfloat>(photon.position.x), static_cast<GLfloat>(photon.position.y),
					   static_cast<GLfloat>(photon.position.z));
		}
	}
	glEnd();

	// Draw photon direction indicators
	glBegin(GL_LINES);
	for (const auto& photon : simulator_->photons) {
		if (photon.alive && photon.weight > 0.001) {
			// Same energy-based color as above
			float weight = static_cast<float>(photon.weight);
			float energy = std::max(0.1f, weight);

			if (energy > 0.7f) {
				glColor3f(1.0f, 1.0f, 0.8f + 0.2f * energy);
			}
			else if (energy > 0.4f) {
				float t = (energy - 0.4f) / 0.3f;
				glColor3f(1.0f, 0.5f + 0.5f * t, 0.2f * t);
			}
			else {
				float t = energy / 0.4f;
				glColor3f(0.5f + 0.5f * t, 0.1f * t, 0.1f * t);
			}

			// Draw direction vector
			Point3 end_pos;
			end_pos.x = photon.position.x + photon.direction.x * 0.08; // Direction indicator
			end_pos.y = photon.position.y + photon.direction.y * 0.08;
			end_pos.z = photon.position.z + photon.direction.z * 0.08;

			glVertex3f(static_cast<GLfloat>(photon.position.x), static_cast<GLfloat>(photon.position.y),
					   static_cast<GLfloat>(photon.position.z));
			glVertex3f(static_cast<GLfloat>(end_pos.x), static_cast<GLfloat>(end_pos.y),
					   static_cast<GLfloat>(end_pos.z));
		}
	}
	glEnd();

	// Draw photon path histories with energy-based coloring
	glLineWidth(2.0f);
	glBegin(GL_LINES);

	for (const auto& path : simulator_->paths) {
		if (path.head) {
			// Draw connected line segments with energy gradient
			Vertex* current = path.head;
			Vertex* next = current ? current->next : nullptr;

			while (current && next) {
				// Color based on vertex value (energy/weight at that point)
				float energy1 = std::max(0.1f, static_cast<float>(current->value));
				float energy2 = std::max(0.1f, static_cast<float>(next->value));

				// Start color (current vertex)
				if (energy1 > 0.7f) {
					glColor3f(1.0f, 1.0f, 0.9f);
				}
				else if (energy1 > 0.4f) {
					float t = (energy1 - 0.4f) / 0.3f;
					glColor3f(1.0f, 0.6f + 0.4f * t, 0.3f * t);
				}
				else {
					float t = energy1 / 0.4f;
					glColor3f(0.6f + 0.4f * t, 0.2f * t, 0.2f * t);
				}
				glVertex3f(static_cast<GLfloat>(current->x), static_cast<GLfloat>(current->y),
						   static_cast<GLfloat>(current->z));

				// End color (next vertex)
				if (energy2 > 0.7f) {
					glColor3f(1.0f, 1.0f, 0.9f);
				}
				else if (energy2 > 0.4f) {
					float t = (energy2 - 0.4f) / 0.3f;
					glColor3f(1.0f, 0.6f + 0.4f * t, 0.3f * t);
				}
				else {
					float t = energy2 / 0.4f;
					glColor3f(0.6f + 0.4f * t, 0.2f * t, 0.2f * t);
				}
				glVertex3f(static_cast<GLfloat>(next->x), static_cast<GLfloat>(next->y), static_cast<GLfloat>(next->z));

				// Move to next segment
				current = next;
				next = current->next;
			}

			// Also draw emitted paths if they exist
			current = path.head;
			while (current) {
				if (current->emit) {
					glColor3f(1.0f, 0.8f, 0.0f); // Bright yellow for emitted photons
					glVertex3f(static_cast<GLfloat>(current->x), static_cast<GLfloat>(current->y),
							   static_cast<GLfloat>(current->z));
					glVertex3f(static_cast<GLfloat>(current->emit->x), static_cast<GLfloat>(current->emit->y),
							   static_cast<GLfloat>(current->emit->z));
				}
				current = current->next;
			}
		}
	}

	glEnd();
	glLineWidth(1.0f);
}

// Simple renderer initialization
void Renderer::initialize(int /*argc*/, char* /*argv*/[]) {
	gl_window_ = std::make_unique<GLWindow>();
	if (!gl_window_->create("Photron - Monte Carlo Voxel Renderer")) {
		std::cout << "Failed to create OpenGL window - continuing in headless mode" << std::endl;
		gl_window_.reset();
		return;
	}

	std::cout << "OpenGL window created successfully!" << std::endl;

	// Initialize OpenGL renderer
	try {
		initialize();
		std::cout << "OpenGL renderer initialized successfully!" << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Renderer initialization failed: " << e.what() << std::endl;
	}

	// Initialize camera parameters
	camera_distance_ = 5.0f;
	camera_rotation_x_ = 30.0f;
	camera_rotation_y_ = 45.0f;
}

void Renderer::render(Simulator* simulator) {
	simulator_ = simulator;

	if (!gl_window_) {
		std::cout << "Headless mode: Simulation complete, visual rendering disabled" << std::endl;
		std::cout << "(OpenGL window creation failed - this is normal in some environments)" << std::endl;
		return;
	}

	// Set up camera
	update_camera();

	int frame_count = 0;

	// Main rendering loop
	while (!gl_window_->should_close()) {
		frame_count++;

		// Poll for window events and handle mouse input
		gl_window_->poll_events();
		handle_mouse_input();

		// Update camera based on mouse input
		update_camera();

		// Clear screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Draw scene in proper order: opaque objects first, then transparent
		draw_bounds();  // Opaque wireframe
		draw_layers();  // Opaque wireframes
		draw_voxels();  // Transparent voxels (with depth mask disabled)
		draw_photons(); // Opaque photon paths (drawn last to be visible)

		// Swap buffers to display
		gl_window_->swap_buffers();
	}
}

void Renderer::update_camera() {
	if (!simulator_) {
		return;
	}

	Range3 bounds = simulator_->bounds;
	double center_x = (bounds.x_min + bounds.x_max) / 2.0;
	double center_y = (bounds.y_min + bounds.y_max) / 2.0;
	double center_z = (bounds.z_min + bounds.z_max) / 2.0;
	double size = std::max({bounds.x_max - bounds.x_min, bounds.y_max - bounds.y_min, bounds.z_max - bounds.z_min});

	// Use spherical coordinates for camera position
	float radius_x = glm::radians(camera_rotation_x_);
	float radius_y = glm::radians(camera_rotation_y_);

	float distance = camera_distance_ * static_cast<float>(size);
	glm::vec3 eye(center_x + distance * cos(radius_x) * cos(radius_y), center_y + distance * sin(radius_x),
				  center_z + distance * cos(radius_x) * sin(radius_y));

	glm::vec3 center(center_x, center_y, center_z);
	glm::vec3 up(0.0f, 1.0f, 0.0f);

	// Set up view matrix
	view_matrix_ = glm::lookAt(eye, center, up);

	// Set up projection matrix
	projection_matrix_ = glm::perspective(glm::radians(45.0f), 1200.0f / 800.0f, 0.1f, 1000.0f);

	// Load matrices into OpenGL
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projection_matrix_));

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(view_matrix_));
}

void Renderer::handle_mouse_input() {
	if (!gl_window_) {
		return;
	}

	const auto& mouse_state = gl_window_->get_mouse_state();

	// Left mouse button: rotate camera
	if (mouse_state.left_button_down) {
		camera_rotation_y_ += static_cast<float>(mouse_state.delta_x) * 0.5f;
		camera_rotation_x_ += static_cast<float>(mouse_state.delta_y) * 0.5f;

		// Clamp vertical rotation
		camera_rotation_x_ = std::max(-89.0f, std::min(89.0f, camera_rotation_x_));

		// Wrap horizontal rotation
		while (camera_rotation_y_ > 360.0f) {
			camera_rotation_y_ -= 360.0f;
		}
		while (camera_rotation_y_ < 0.0f) {
			camera_rotation_y_ += 360.0f;
		}
	}

	// Right mouse button: zoom (like in VolRec)
	if (mouse_state.right_button_down) {
		camera_distance_ *= (1.0f + static_cast<float>(mouse_state.delta_y) * 0.01f);
		camera_distance_ = std::max(0.1f, std::min(20.0f, camera_distance_));
	}

	// Mouse wheel: also zoom (alternative method)
	if (mouse_state.scroll_delta != 0.0f) {
		camera_distance_ *= (1.0f - mouse_state.scroll_delta * 0.1f);
		camera_distance_ = std::max(0.1f, std::min(20.0f, camera_distance_));
	}
}
