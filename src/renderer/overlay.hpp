/**
 * @file overlay.hpp
 * @brief ImGui-based user interface overlay for simulation control
 *
 * Provides comprehensive GUI controls for Monte Carlo simulation parameters,
 * real-time visualization settings, and interactive parameter adjustment.
 */

#pragma once

#include <functional>
#include <optional>
#include <string>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "renderer/settings.hpp"
#include "simulator/metrics.hpp"

// Forward declarations
class Simulator;

// World text overlay system
struct WorldText
{
	std::string text;
	float screen_x, screen_y;
	glm::vec4 color;
};

class Overlay
{
public:
	/**
	 * @brief Construct new Overlay with default settings
	 */
	Overlay();

	/**
	 * @brief Destroy Overlay and clean up ImGui resources
	 */
	~Overlay();

	/**
	 * @brief Initialize ImGui context and rendering backend
	 *
	 * Sets up ImGui for OpenGL rendering with GLFW input integration.
	 * Must be called after valid OpenGL context is established.
	 *
	 * @param window GLFW window handle for input integration
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize(GLFWwindow* window);

	/**
	 * @brief Shutdown ImGui and clean up resources
	 *
	 * Should be called before OpenGL context destruction.
	 */
	void shutdown();

	/**
	 * @brief Render basic overlay without simulation data
	 *
	 * Renders minimal UI when no simulator is available,
	 * typically just file operations and basic controls.
	 */
	void render();

	/**
	 * @brief Render complete overlay with simulation integration
	 *
	 * Renders full interface including simulation controls, progress
	 * monitoring, energy statistics, and real-time data display.
	 *
	 * @param simulator Reference to active simulator instance
	 */
	void render_with_simulator(Simulator& simulator);

	// Settings access

	/**
	 * @brief Get current rendering settings (read-only)
	 * @return const Settings& Current settings configuration
	 */
	const Settings& get_settings() const { return settings_; }

	/**
	 * @brief Get current rendering settings (mutable)
	 * @return Settings& Settings configuration for modification
	 */
	Settings& get_settings() { return settings_; }

	/**
	 * @brief Set callback for configuration file loading operations
	 * @param callback Function to call with selected config file path
	 */
	void set_open_config_callback(std::function<void(const std::string&)> callback) {
		open_config_callback_ = callback;
	}

	/**
	 * @brief Set callback for starting simulation
	 * @param callback Function to call when user requests simulation start
	 */
	void set_run_simulation_callback(std::function<void()> callback) { run_simulation_callback_ = callback; }

	/**
	 * @brief Set callback for restarting simulation with same parameters
	 * @param callback Function to call when user requests simulation restart
	 */
	void set_rerun_simulation_callback(std::function<void()> callback) { rerun_simulation_callback_ = callback; }

	/**
	 * @brief Set callback for saving simulation results to file
	 * @param callback Function to call with target file path for results
	 */
	void set_save_results_callback(std::function<void(const std::string&)> callback) {
		save_results_callback_ = callback;
	}

	/**
	 * @brief Set callback for direct save operation (auto-generated filename)
	 * @param callback Function to call for automatic result saving
	 */
	void set_direct_save_results_callback(std::function<void()> callback) { direct_save_results_callback_ = callback; }

	/**
	 * @brief Set callback for camera view reset operations
	 * @param callback Function to call when user requests view reset
	 */
	void set_reset_view_callback(std::function<void()> callback) { reset_view_callback_ = callback; }

	/**
	 * @brief Set callback for camera mode change notifications
	 * @param callback Function to call with new FPS mode state (true = FPS, false = orbital)
	 */
	void set_camera_mode_changed_callback(std::function<void(bool)> callback) {
		camera_mode_changed_callback_ = callback;
	}

	/**
	 * @brief Set callback for closing current project
	 * @param callback Function to call when user requests to close project and return to config selection
	 */
	void set_close_project_callback(std::function<void()> callback) { close_project_callback_ = callback; }

	/**
	 * @brief Enable or disable UI interaction (used during simulation)
	 * @param enabled True to enable UI, false to disable during simulation
	 */
	void set_ui_enabled(bool enabled) { ui_enabled_ = enabled; }

	/**
	 * @brief Check if UI interaction is currently enabled
	 * @return bool True if UI is responsive, false if disabled
	 */
	bool is_ui_enabled() const { return ui_enabled_; }

	/**
	 * @brief Check if file selection dialog is currently displayed
	 * @return bool True if file dialog is open, false otherwise
	 */
	bool is_file_dialog_open() const { return show_file_dialog_; }

	/**
	 * @brief Open configuration file dialog immediately
	 */
	void open_config_dialog();

	/**
	 * @brief Schedule configuration dialog to open on next frame
	 */
	void open_config_dialog_deferred();

	/**
	 * @brief Add text overlay at screen coordinates for 3D world rendering
	 * @param text Text content to display
	 * @param screen_x Screen X coordinate for text placement
	 * @param screen_y Screen Y coordinate for text placement
	 * @param color RGBA color vector for text rendering
	 */
	void add_world_text(const std::string& text, float screen_x, float screen_y, const glm::vec4& color);

	/**
	 * @brief Clear all queued world text overlays
	 */
	void clear_world_text();

	/**
	 * @brief Display temporary feedback message for save operations
	 * @param message Feedback text to display to user
	 */
	void show_save_feedback(const std::string& message);

private:
	void render_main_menu_bar(bool has_simulator = true);
	void render_control_panel(Simulator* simulator = nullptr);
	void render_file_dialog();
	void render_world_text_overlays();
	void render_save_feedback();
	void handle_keyboard_shortcuts(bool has_simulator = false);

	Settings settings_; ///< Current rendering and interaction settings
	bool ui_enabled_;   ///< Flag controlling overall UI responsiveness
	bool fullscreen_mode_; ///< Flag to hide all UI elements (fullscreen mode)

	// Save feedback state
	bool show_save_feedback_;           ///< Flag to display save operation status
	std::string save_feedback_message_; ///< Message text for save feedback display
	float save_feedback_timer_;         ///< Timer for auto-hiding save feedback

	// File dialog state
	bool show_file_dialog_;         ///< Flag to display file selection dialog
	bool deferred_open_config_;     ///< Flag for deferred dialog opening
	char file_path_buffer_[512];    ///< C-style buffer for ImGui file dialog
	std::string file_path_;         ///< Modern string for file path storage
	std::string current_directory_; ///< Current working directory for file operations

	// UI interaction callbacks
	std::function<void(const std::string&)> open_config_callback_;  ///< Callback for configuration file loading
	std::function<void()> run_simulation_callback_;                 ///< Callback for starting simulation
	std::function<void()> rerun_simulation_callback_;               ///< Callback for restarting simulation
	std::function<void(const std::string&)> save_results_callback_; ///< Callback for saving results to file
	std::function<void()> direct_save_results_callback_;            ///< Callback for direct save operation
	std::function<void()> reset_view_callback_;                     ///< Callback for camera view reset
	std::function<void(bool)> camera_mode_changed_callback_;        ///< Callback for camera mode changes
	std::function<void()> close_project_callback_;                  ///< Callback for closing current project
	std::vector<WorldText> world_text_queue_;                       ///< Queue of 3D world text overlays to render

	// Performance optimization: Energy data caching
	mutable std::optional<Metrics::EnergyDisplayData> cached_energy_data_; ///< Cached energy statistics to avoid per-frame recalculation
	mutable uint64_t cached_simulation_version_ {0}; ///< Last simulation version for cache invalidation

	// Helper method to get cached energy data with event-driven updates
	const Metrics::EnergyDisplayData& get_cached_energy_data(Simulator* simulator) const;
};
