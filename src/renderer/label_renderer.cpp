/**
 * @file label_renderer.cpp
 * @brief Implementation of specialized energy label and text overlay system
 */

#include "label_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <tuple>

#include "math/math.hpp"
#include "renderer/camera.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"

LabelRenderer::LabelRenderer() = default;

bool LabelRenderer::initialize() {
	// Label renderer doesn't need OpenGL resources
	return true;
}

void LabelRenderer::render(const Simulator& simulator, const Settings& settings, const Camera& camera) {
	if (!settings.draw_labels) {
		return;
	}

	// Cache labels if not already cached
	if (!labels_cached_) {
		cache_energy_labels(simulator);
	}

	// Always update screen positions to ensure labels follow camera movement
	update_screen_positions(camera);

	// Render cached labels using callback
	if (text_render_callback_) {
		for (const auto& label : cached_labels_) {
			if (label.screen_position_valid) {
				text_render_callback_(label.text, label.screen_position.x, label.screen_position.y, label.color);
			}
		}
	}
}

void LabelRenderer::set_viewport(int width, int height) {
	viewport_width_ = width;
	viewport_height_ = height;
	camera_changed_ = true; // Force screen position recalculation
}

void LabelRenderer::invalidate_cache() {
	labels_cached_ = false;
	cached_labels_.clear();
	camera_changed_ = true;
}

void LabelRenderer::auto_manage_labels(Settings& settings, size_t photon_count) {
	bool many_photons = (photon_count > 10);

	// Auto-disable when crossing the 10 photon threshold
	if (many_photons && !auto_disabled_labels_) {
		settings.draw_labels = false;
		auto_disabled_labels_ = true;
	}
	// Reset auto-disable flag when photon count drops back down
	else if (!many_photons && auto_disabled_labels_) {
		auto_disabled_labels_ = false;
	}
}

void LabelRenderer::cache_energy_labels(const Simulator& simulator) {
	cached_labels_.clear();
	labels_cached_ = false;

	// ORIGINAL EMITTER-BASED ENERGY LABEL SYSTEM RESTORED
	// Group emitters by position to avoid overlapping labels like the original

	const auto& emitters = simulator.emitters;
	if (emitters.empty()) {
		labels_cached_ = true;
		return;
	}

	// Group emitters by approximate position to avoid overlapping labels (ORIGINAL LOGIC)
	std::map<std::tuple<int, int, int>, std::vector<std::shared_ptr<Emitter>>> position_groups;

	for (const auto& emitter : emitters) {
		if (emitter->weight < 0.001)
			continue; // Skip very low energy exits

		// Group by position rounded to nearest 0.01 units (ORIGINAL POSITIONING)
		int pos_x = static_cast<int>(std::round(emitter->position.x * 100));
		int pos_y = static_cast<int>(std::round(emitter->position.y * 100));
		int pos_z = static_cast<int>(std::round(emitter->position.z * 100));
		auto pos_key = std::make_tuple(pos_x, pos_y, pos_z);

		position_groups[pos_key].push_back(emitter);
	}

	// Create labels for each position group (ORIGINAL COMBINING LOGIC)
	for (const auto& [pos_key, emitters_at_pos] : position_groups) {
		if (emitters_at_pos.empty())
			continue;

		// Use the first emitter's position as the label position (ORIGINAL POSITIONING)
		const auto& representative = emitters_at_pos[0];
		glm::vec3 label_pos(static_cast<float>(representative->position.x),
							static_cast<float>(representative->position.y),
							static_cast<float>(representative->position.z));

		// Sum up energy and classify by exit type (ORIGINAL CLASSIFICATION)
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

		// Create label with proper classification and NORMALIZED percentages (ORIGINAL LOGIC)
		std::string label_text;
		glm::vec4 label_color(1.0f, 1.0f, 1.0f, 1.0f);

		double total_energy = total_reflected_energy + total_transmitted_energy + total_unclassified_energy;

		// ORIGINAL: Normalize by total number of photons for percentage calculation
		double total_initial_energy = static_cast<double>(simulator.photons.size());
		double energy_percent = (total_energy / total_initial_energy) * 100.0;

		// Format percentage, showing "<1%" instead of "0%" for very small values (ORIGINAL FORMAT)
		int rounded_percent = static_cast<int>(energy_percent);
		std::string percent_text =
			(rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);

		if (total_reflected_energy > total_transmitted_energy && total_reflected_energy > total_unclassified_energy) {
			// Predominantly reflected (ORIGINAL BRIGHT GREEN)
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.8f, 0.2f, 1.0f); // Bright green for reflection
		}
		else if (total_transmitted_energy > total_reflected_energy
				 && total_transmitted_energy > total_unclassified_energy) {
			// Predominantly transmitted (ORIGINAL BRIGHT BLUE)
			label_text = percent_text;
			label_color = glm::vec4(0.2f, 0.6f, 1.0f, 1.0f); // Bright blue for transmission
		}
		else {
			// Mixed or unclassified (ORIGINAL GRAY)
			label_text = percent_text;
			label_color = glm::vec4(0.8f, 0.8f, 0.8f, 1.0f); // Gray for mixed/unclassified
		}

		cached_labels_.push_back({.world_position = label_pos,
								  .text = label_text,
								  .color = label_color,
								  .scale = 1.0f,
								  .screen_position = glm::vec2(0.0f),
								  .screen_position_valid = false});
	}

	// Handle specular surface reflection if present (ORIGINAL SPECULAR LOGIC)
	auto energy_data = simulator.get_metrics().aggregate_medium_energy_data(simulator);
	double specular_reflection = energy_data.specular_reflection;
	if (specular_reflection > 0.0 && !simulator.sources.empty()) {
		const Source& source = simulator.sources[0];
		double surface_refraction = energy_data.surface_refraction;
		double energy_percent = (specular_reflection / surface_refraction) * 100.0;

		// Format percentage like original
		int rounded_percent = static_cast<int>(energy_percent);
		std::string surface_percent_text =
			(rounded_percent == 0 && energy_percent > 0.0) ? "<1%" : std::format("{}%", rounded_percent);

		// Position label at the tip of the specular reflection vector (ORIGINAL POSITIONING)
		glm::vec3 surface_label_pos(static_cast<float>(source.intersect.x + source.specular_direction.x * 0.1),
									static_cast<float>(source.intersect.y + source.specular_direction.y * 0.1),
									static_cast<float>(source.intersect.z + source.specular_direction.z * 0.1));

		// Use bright purple color for surface specular reflection (ORIGINAL BRIGHT PURPLE)
		glm::vec4 surface_label_color = glm::vec4(0.8f, 0.3f, 1.0f, 1.0f);

		cached_labels_.push_back({.world_position = surface_label_pos,
								  .text = surface_percent_text,
								  .color = surface_label_color,
								  .scale = 1.0f,
								  .screen_position = glm::vec2(0.0f),
								  .screen_position_valid = false});
	}

	labels_cached_ = true;
}

void LabelRenderer::update_screen_positions(const Camera& camera) {
	for (auto& label : cached_labels_) {
		label.screen_position = world_to_screen(label.world_position, camera);
		label.screen_position_valid = (label.screen_position.x >= 0 && label.screen_position.x < viewport_width_
									   && label.screen_position.y >= 0 && label.screen_position.y < viewport_height_);
	}
}

glm::vec2 LabelRenderer::world_to_screen(const glm::vec3& world_pos, const Camera& camera) const {
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
	float screen_x = (ndc.x + 1.0f) * 0.5f * viewport_width_;
	float screen_y = (1.0f - ndc.y) * 0.5f * viewport_height_; // Flip Y for screen coordinates

	return glm::vec2(screen_x, screen_y);
}

glm::vec3 LabelRenderer::get_voxel_world_position(const Voxel& voxel, const Simulator& simulator) const {
	// Get voxel center position in world coordinates
	// This assumes uniform voxel size and grid-to-world coordinate mapping

	// Get combined bounds from simulator
	auto bounds = simulator.get_combined_bounds();

	// Calculate voxel size (assuming uniform cubic voxels)
	// Note: This should match the voxel size used in the simulator
	double voxel_size = 0.1; // Default voxel size - should be retrieved from config

	// Calculate world position from grid coordinates
	glm::vec3 world_pos;
	world_pos.x = static_cast<float>(bounds.min_bounds.x + (voxel.ix() + 0.5) * voxel_size);
	world_pos.y = static_cast<float>(bounds.min_bounds.y + (voxel.iy() + 0.5) * voxel_size);
	world_pos.z = static_cast<float>(bounds.min_bounds.z + (voxel.iz() + 0.5) * voxel_size);

	return world_pos;
}
