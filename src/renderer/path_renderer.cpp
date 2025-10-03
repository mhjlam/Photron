/**
 * @file path_renderer.cpp
 * @brief Complete implementation of self-contained photon path rendering system
 *
 * Transfers all photon path visualization functionality from Renderer::draw_paths()
 * into a specialized, self-contained rendering system with full OpenGL resource
 * management and sophisticated caching for optimal performance.
 */

#include "path_renderer.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>

#include "renderer/camera.hpp"
#include "renderer/shader.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"

PathRenderer::PathRenderer() {
	if (!initialize()) {
		std::cerr << "PathRenderer: Failed to initialize OpenGL resources" << std::endl;
	}
}

PathRenderer::~PathRenderer() {
	// Clean up OpenGL resources
	if (lines_vao_) {
		glDeleteVertexArrays(1, &lines_vao_);
	}
	if (lines_vbo_) {
		glDeleteBuffers(1, &lines_vbo_);
	}
	if (lines_instance_vbo_) {
		glDeleteBuffers(1, &lines_instance_vbo_);
	}
	if (lines_shader_) {
		glDeleteProgram(lines_shader_);
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
	if (points_shader_) {
		glDeleteProgram(points_shader_);
	}
}

bool PathRenderer::initialize() {
	return setup_line_instanced_rendering() && setup_point_instanced_rendering();
}

void PathRenderer::render_paths(const Settings& settings, const Simulator& simulator, const Camera& camera) {
	if (!settings.draw_paths) {
		std::cout << "PathRenderer: draw_paths disabled in settings" << std::endl;
		return;
	}

	// Use cached energy range for better performance
	update_cached_energy_range(simulator);

	// Create adaptive logarithmic mapping function using cached values
	auto adaptive_log_color = [this](float energy) -> glm::vec4 {
		return get_adaptive_energy_color(energy, cached_min_energy_, cached_max_energy_);
	};

	// Use incremental caching instead of rebuilding every frame
	size_t current_photon_count = simulator.photons.size();

	// Check if we need to rebuild cache completely or just add new photons
	if (!path_instances_cached_ || current_photon_count < cached_photon_count_) {
		// Complete rebuild needed (first time or photons were removed)
		cached_line_instances_.clear();
		cached_point_instances_.clear();
		cached_photon_count_ = 0;

		// Reset buffer upload flags when cache is invalidated
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;
	}
	else if (current_photon_count > cached_photon_count_) {
		// Incremental update - only process new photons
		// Keep existing cache, just mark buffers for re-upload
		line_buffer_uploaded_ = false;
		point_buffer_uploaded_ = false;
	}

	// Only process photons if we have new ones to add
	if (current_photon_count > cached_photon_count_) {
		// Update cached surface calculation
		update_cached_surface(simulator);

		// Collect line instances for photon paths
		collect_path_line_instances(simulator, adaptive_log_color);

		// Collect emitter direction vectors
		collect_emitter_vectors(simulator, adaptive_log_color);

		// Collect scatter and emitter points
		collect_path_point_instances(simulator, adaptive_log_color);

		// Update cached photon count and mark cache as up to date
		cached_photon_count_ = current_photon_count;
		path_instances_cached_ = true;
	}

	// Render cached instances
	render_line_instances(camera);
	render_point_instances(camera);
}

void PathRenderer::invalidate_cache() {
	path_instances_cached_ = false;
	cached_line_instances_.clear();
	cached_point_instances_.clear();
	cached_photon_count_ = 0;
	line_buffer_uploaded_ = false;
	point_buffer_uploaded_ = false;
	energy_range_cached_ = false;
	surface_cached_ = false;
}

void PathRenderer::update_cached_energy_range(const Simulator& simulator) {
	if (!energy_range_cached_ || simulator.photons.size() != cached_photon_count_) {
		cached_min_energy_ = 1.0f;
		cached_max_energy_ = 0.0f;

		for (const Photon& photon : simulator.photons) {
			if (photon.path_head) {
				auto current = photon.path_head;
				while (current) {
					float energy = static_cast<float>(current->value);
					cached_min_energy_ = std::min(cached_min_energy_, energy);
					cached_max_energy_ = std::max(cached_max_energy_, energy);
					current = current->next;
				}
			}
		}

		// Ensure valid range
		if (cached_min_energy_ >= cached_max_energy_) {
			cached_min_energy_ = 0.001f;
			cached_max_energy_ = 1.0f;
		}

		energy_range_cached_ = true;
	}
}

void PathRenderer::update_cached_surface(const Simulator& simulator) {
	if (!surface_cached_) {
		const auto& layers = simulator.get_all_layers();
		if (!layers.empty()) {
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
	}
}

void PathRenderer::collect_path_line_instances(const Simulator& simulator,
											   const std::function<glm::vec4(float)>& adaptive_log_color) {
	// Process only new photons (incremental caching)
	const auto& paths = simulator.photons;
	for (size_t i = cached_photon_count_; i < paths.size(); ++i) {
		const Photon& photon = paths[i];
		if (photon.path_head) {
			// Generate connected line segments with energy gradient
			auto current = photon.path_head;
			auto next = current ? current->next : nullptr;

			// Generate incident ray from source to material surface (cached surface)
			if (current && !simulator.sources.empty()) {
				glm::vec3 first_interaction(static_cast<float>(current->position.x),
											static_cast<float>(current->position.y),
											static_cast<float>(current->position.z));

				const auto& source = simulator.sources[0];
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
}
void PathRenderer::collect_emitter_vectors(const Simulator& simulator,
										   const std::function<glm::vec4(float)>& adaptive_log_color) {
	// Add all emitter direction vectors to instanced rendering cache with deduplication
	if (!simulator.emitters.empty()) {
		// Group emitters by origin position and direction (with margin for similar directions)
		std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::shared_ptr<Emitter>>> direction_groups;
		const double POSITION_PRECISION = 100.0; // Group positions within 0.01 units
		const double DIRECTION_PRECISION = 50.0; // Group directions within ~0.02 radians (~1.1 degrees)

		for (const auto& emitter : simulator.emitters) {
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

			// Use actual averaged energy for consistent color mapping with photon paths
			float averaged_energy = static_cast<float>(total_weight / grouped_emitters.size());
			glm::vec4 direction_color = adaptive_log_color(averaged_energy);
			direction_color.a = 0.8f; // Slightly transparent for distinction

			cached_line_instances_.push_back({start_pos, end_pos, direction_color, direction_color});
		}
	}

	// Add specular reflection emitter for incident photon's surface reflection
	// This is integrated into the direction grouping above to avoid duplication
	double specular_reflection = simulator.get_metrics().aggregate_medium_energy_data(simulator).specular_reflection;
	if (specular_reflection > 0.0 && !simulator.sources.empty()) {
		const auto& source = simulator.sources[0];

		// Check if there are any regular emitters at the same position with similar direction
		bool found_similar_emitter = false;
		glm::dvec3 specular_pos = source.intersect;
		glm::dvec3 specular_dir = glm::normalize(source.specular_direction);

		const double POSITION_TOLERANCE = 0.01; // Same as POSITION_PRECISION above
		const double ANGLE_TOLERANCE = 0.02;    // Same as DIRECTION_PRECISION above

		for (const auto& emitter : simulator.emitters) {
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
			glm::vec4 specular_color = adaptive_log_color(specular_energy);
			specular_color.a = 1.0f; // Full opacity for better visibility

			cached_line_instances_.push_back({reflection_start, reflection_end, specular_color, specular_color});
		}
	}
}

void PathRenderer::collect_path_point_instances(const Simulator& simulator,
												const std::function<glm::vec4(float)>& adaptive_log_color) {
	// First, collect scatter points from photon paths
	for (const Photon& photon : simulator.photons) {
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

	// Second, add emitter exit points (these are the accurate surface boundary points)
	if (!simulator.emitters.empty()) {
		// Add emitter points to the point instances vector
		for (const auto& emitter : simulator.emitters) {
			// Use emitter->position which contains the corrected surface intersection coordinates
			glm::vec3 exit_pos(static_cast<float>(emitter->position.x),
							   static_cast<float>(emitter->position.y),
							   static_cast<float>(emitter->position.z));

			// Use actual emitter weight for proper energy-based coloring
			// This ensures emitters match the energy mapping of photon paths
			float emitter_energy = static_cast<float>(emitter->weight);
			glm::vec4 exit_color = adaptive_log_color(emitter_energy);
			exit_color.a = 1.0f; // Full opacity for emitter points

			// Add emitter points to point instances
			PointInstance point_instance;
			point_instance.position = exit_pos;
			point_instance.color = exit_color;
			point_instance.size = 6.0f; // Same size as scatter points
			cached_point_instances_.push_back(point_instance);
		}

		// Add surface specular reflection as an emitter point
		double surface_specular_reflection =
			simulator.get_metrics().aggregate_medium_energy_data(simulator).specular_reflection;
		if (surface_specular_reflection > 0.0 && !simulator.sources.empty()) {
			const auto& source = simulator.sources[0];

			// Use actual source intersection point
			glm::vec3 surface_entry(static_cast<float>(source.intersect.x),
									static_cast<float>(source.intersect.y),
									static_cast<float>(source.intersect.z));

			// Use actual specular reflection energy for consistent color mapping
			float specular_energy = static_cast<float>(surface_specular_reflection);
			glm::vec4 surface_point_color = adaptive_log_color(specular_energy);
			surface_point_color.a = 1.0f; // Full opacity for surface reflection point

			PointInstance surface_point;
			surface_point.position = surface_entry;
			surface_point.color = surface_point_color;
			surface_point.size = 8.0f;    // Slightly larger for surface reflection point
			cached_point_instances_.push_back(surface_point);
		}
	}
}

void PathRenderer::render_line_instances(const Camera& camera) {
	// Render cached line instances
	if (!cached_line_instances_.empty() && lines_shader_) {
		// Set line width for visibility (match original Renderer)
		glLineWidth(3.0f);

		// Ensure proper OpenGL state for line rendering (match original Renderer)
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS); // Match original Renderer depth function
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glUseProgram(lines_shader_);
		current_shader_program_ = lines_shader_;

		// Use cached uniform location, compute MVP directly
		glm::mat4 mvp = camera.get_projection_matrix() * camera.get_view_matrix();
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
		current_shader_program_ = 0;

		// Reset line width to default
		glLineWidth(1.0f);
	}
}

void PathRenderer::render_point_instances(const Camera& camera) {
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
		glm::mat4 mvp = camera.get_projection_matrix() * camera.get_view_matrix();
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
		current_shader_program_ = 0;

		// Restore OpenGL state once after rendering
		disable_blending(); // Use state management helper
		glDisable(GL_PROGRAM_POINT_SIZE);
	}
}
// ========================================
// OPENGL SETUP AND UTILITY METHODS
// ========================================

bool PathRenderer::setup_line_instanced_rendering() {
	// Load line shaders
	std::string vertex_source = Shader::load_shader_source("shaders/lines.vert");
	std::string fragment_source = Shader::load_shader_source("shaders/lines.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "PathRenderer: Failed to load line shaders" << std::endl;
		return false;
	}

	// Create line shader program using Shader class methods
	GLuint vertex_shader = Shader::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = Shader::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader)
			glDeleteShader(vertex_shader);
		if (fragment_shader)
			glDeleteShader(fragment_shader);
		std::cerr << "PathRenderer: Failed to compile line shaders" << std::endl;
		return false;
	}

	lines_shader_ = glCreateProgram();
	glAttachShader(lines_shader_, vertex_shader);
	glAttachShader(lines_shader_, fragment_shader);
	glLinkProgram(lines_shader_);

	GLint success;
	glGetProgramiv(lines_shader_, GL_LINK_STATUS, &success);
	if (!success) {
		std::array<char, 512> info_log {};
		glGetProgramInfoLog(lines_shader_, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "PathRenderer: Line shader linking failed: " << info_log.data() << std::endl;
		glDeleteProgram(lines_shader_);
		lines_shader_ = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	if (lines_shader_ == 0) {
		return false;
	}

	// Cache uniform locations
	line_instanced_mvp_uniform_location_ = glGetUniformLocation(lines_shader_, "uMVP");

	if (line_instanced_mvp_uniform_location_ == -1) {
		std::cerr << "PathRenderer: Failed to find uMVP uniform in lines shader" << std::endl;
	}

	// Create VAO and VBOs
	glGenVertexArrays(1, &lines_vao_);
	glGenBuffers(1, &lines_vbo_);
	glGenBuffers(1, &lines_instance_vbo_);

	glBindVertexArray(lines_vao_);

	// Setup base line geometry (parameter values for shader interpolation)
	glBindBuffer(GL_ARRAY_BUFFER, lines_vbo_);
	float line_vertices[] = {
		0.0f,
		0.0f,
		0.0f, // t=0 (start point parameter)
		1.0f,
		0.0f,
		0.0f  // t=1 (end point parameter)
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(line_vertices), line_vertices, GL_STATIC_DRAW);

	// Vertex position attribute (location 0) - shader uses aPosition.x for interpolation
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Setup instance data buffer
	glBindBuffer(GL_ARRAY_BUFFER, lines_instance_vbo_);

	// Instance start position (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// Instance end position (location 2)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// Instance start color (location 3)
	glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, start_color));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	// Instance end color (location 4)
	glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(LineInstance), (void*)offsetof(LineInstance, end_color));
	glEnableVertexAttribArray(4);
	glVertexAttribDivisor(4, 1);

	glBindVertexArray(0);
	return true;
}

bool PathRenderer::setup_point_instanced_rendering() {
	// Load point shaders
	std::string vertex_source = Shader::load_shader_source("shaders/points.vert");
	std::string fragment_source = Shader::load_shader_source("shaders/points.frag");

	if (vertex_source.empty() || fragment_source.empty()) {
		std::cerr << "PathRenderer: Failed to load point shaders" << std::endl;
		return false;
	}

	// Create point shader program using Shader class methods
	GLuint vertex_shader = Shader::compile_shader(vertex_source, GL_VERTEX_SHADER);
	GLuint fragment_shader = Shader::compile_shader(fragment_source, GL_FRAGMENT_SHADER);

	if (vertex_shader == 0 || fragment_shader == 0) {
		if (vertex_shader)
			glDeleteShader(vertex_shader);
		if (fragment_shader)
			glDeleteShader(fragment_shader);
		std::cerr << "PathRenderer: Failed to compile point shaders" << std::endl;
		return false;
	}

	points_shader_ = glCreateProgram();
	glAttachShader(points_shader_, vertex_shader);
	glAttachShader(points_shader_, fragment_shader);
	glLinkProgram(points_shader_);

	GLint success;
	glGetProgramiv(points_shader_, GL_LINK_STATUS, &success);
	if (!success) {
		std::array<char, 512> info_log {};
		glGetProgramInfoLog(points_shader_, static_cast<GLsizei>(info_log.size()), nullptr, info_log.data());
		std::cerr << "PathRenderer: Point shader linking failed: " << info_log.data() << std::endl;
		glDeleteProgram(points_shader_);
		points_shader_ = 0;
	}

	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	if (points_shader_ == 0) {
		return false;
	}

	// Cache uniform locations
	point_instanced_mvp_uniform_location_ = glGetUniformLocation(points_shader_, "uMVP");

	// Create VAO and VBOs
	glGenVertexArrays(1, &points_vao_);
	glGenBuffers(1, &points_vbo_);
	glGenBuffers(1, &points_instance_vbo_);

	glBindVertexArray(points_vao_);

	// Setup base point geometry (single point at origin)
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo_);
	float point_vertex[3] = {0.0f, 0.0f, 0.0f};
	glBufferData(GL_ARRAY_BUFFER, sizeof(point_vertex), point_vertex, GL_STATIC_DRAW);

	// Vertex position attribute (location 0)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Setup instance data buffer
	glBindBuffer(GL_ARRAY_BUFFER, points_instance_vbo_);

	// Instance position (location 1)
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, position));
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// Instance color (location 2)
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, color));
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	// Instance size (location 3)
	glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(PointInstance), (void*)offsetof(PointInstance, size));
	glEnableVertexAttribArray(3);
	glVertexAttribDivisor(3, 1);

	glBindVertexArray(0);
	return true;
}

void PathRenderer::enable_blending() const {
	glEnable(GL_BLEND);
}

void PathRenderer::disable_blending() const {
	glDisable(GL_BLEND);
}

glm::vec4 PathRenderer::get_adaptive_energy_color(float energy, float min_energy, float max_energy) const {
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
