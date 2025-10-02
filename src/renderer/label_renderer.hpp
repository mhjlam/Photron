/**
 * @file label_renderer.hpp
 * @brief Specialized energy label and text overlay system
 *
 * Handles rendering of energy percentage labels, screen-space text positioning,
 * and billboard text rendering through callback integration with ImGui or
 * other text rendering systems.
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "renderer/settings.hpp"

// Forward declarations
class Voxel;

// Forward declarations
class Simulator;
class Camera;

/**
 * @class LabelRenderer
 * @brief Specialized renderer for energy labels and text overlays
 *
 * Handles all aspects of label rendering including:
 * - Energy percentage calculation and formatting
 * - 3D to 2D screen space projection
 * - Adaptive label grouping based on density
 * - Text rendering through external callbacks
 * - Label visibility management
 */
class LabelRenderer
{
public:
	/**
	 * @brief Construct a new LabelRenderer
	 */
	LabelRenderer();

	/**
	 * @brief Destroy the LabelRenderer
	 */
	~LabelRenderer() = default;

	/**
	 * @brief Initialize label rendering system
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize();

	/**
	 * @brief Render energy labels with current settings
	 * @param simulator Simulator containing energy data
	 * @param settings Current rendering settings
	 * @param camera Camera for screen space projection
	 */
	void render(const Simulator& simulator, const Settings& settings, const Camera& camera);

	/**
	 * @brief Update viewport dimensions for screen space calculations
	 * @param width New viewport width
	 * @param height New viewport height
	 */
	void set_viewport(int width, int height);

	/**
	 * @brief Set callback function for text rendering integration
	 * @param callback Function accepting (text, x, y, color) parameters for text display
	 */
	template<typename Callable>
		requires std::invocable<Callable, const std::string&, float, float, const glm::vec4&>
	void set_text_render_callback(Callable&& callback) {
		text_render_callback_ = std::forward<Callable>(callback);
	}

	/**
	 * @brief Invalidate cached label data
	 */
	void invalidate_cache();

	/**
	 * @brief Automatically manage label visibility based on photon count
	 * @param settings Settings object to modify for label visibility
	 * @param photon_count Current number of photons in simulation
	 */
	void auto_manage_labels(Settings& settings, size_t photon_count);

private:
	/**
	 * @brief Energy label structure for billboard text rendering
	 */
	struct EnergyLabel
	{
		glm::vec3 world_position;
		std::string text;
		glm::vec4 color;
		float scale {1.0f};
		glm::vec2 screen_position;          // Cached screen position
		bool screen_position_valid {false}; // Whether screen position is current
	};

	// Label generation
	void cache_energy_labels(const Simulator& simulator);

	// Screen space projection
	void update_screen_positions(const Camera& camera);
	glm::vec2 world_to_screen(const glm::vec3& world_pos, const Camera& camera) const;

	// Voxel position utilities for original energy label system
	glm::vec3 get_voxel_world_position(const Voxel& voxel, const Simulator& simulator) const;

	// Text rendering callback
	std::function<void(const std::string&, float, float, const glm::vec4&)> text_render_callback_;

	// Label cache
	std::vector<EnergyLabel> cached_labels_;
	mutable bool labels_cached_ {false};
	mutable bool camera_changed_ {true};

	// Auto-management state
	bool auto_disabled_labels_ {false};

	// Viewport
	int viewport_width_ {800};
	int viewport_height_ {600};
};
