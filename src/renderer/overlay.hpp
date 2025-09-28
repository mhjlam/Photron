#pragma once

#include <functional>
#include <string>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "renderer/settings.hpp"

class Simulator;

enum class FileDialogMode
{
	LoadConfig,
	SaveResults
};

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
	Overlay();
	~Overlay();

	bool initialize(GLFWwindow* window);
	void shutdown();
	void render();
	void render_with_simulator(Simulator& simulator);

	// Getters for settings
	const Settings& get_settings() const { return settings_; }
	Settings& get_settings() { return settings_; }

	// Callback setters for file operations
	void set_open_config_callback(std::function<void(const std::string&)> callback) {
		open_config_callback_ = callback;
	}
	void set_run_simulation_callback(std::function<void()> callback) { run_simulation_callback_ = callback; }
	void set_rerun_simulation_callback(std::function<void()> callback) { rerun_simulation_callback_ = callback; }
	void set_save_results_callback(std::function<void(const std::string&)> callback) {
		save_results_callback_ = callback;
	}
	void set_direct_save_results_callback(std::function<void()> callback) {
		direct_save_results_callback_ = callback;
	}
	void set_reset_view_callback(std::function<void()> callback) { reset_view_callback_ = callback; }
	void set_camera_mode_changed_callback(std::function<void(bool)> callback) {
		camera_mode_changed_callback_ = callback;
	}

	// Enable/disable UI (for simulation running)
	void set_ui_enabled(bool enabled) { ui_enabled_ = enabled; }
	bool is_ui_enabled() const { return ui_enabled_; }
	
	// Check if file dialog is currently open
	bool is_file_dialog_open() const { return show_file_dialog_; }
	void open_config_dialog(); // Method to auto-open config dialog
	void open_config_dialog_deferred(); // Method to auto-open config dialog in next frame

	// Text overlay system for 3D world text rendering
	void add_world_text(const std::string& text, float screen_x, float screen_y, const glm::vec4& color);
	void clear_world_text();
	
	// Save feedback system
	void show_save_feedback(const std::string& message);

private:
	void render_main_menu_bar(bool has_simulator = true);
	void render_control_panel(Simulator* simulator = nullptr);
	void render_file_dialog();
	void render_world_text_overlays();
	void render_save_feedback();
	void handle_keyboard_shortcuts(bool has_simulator = false);

	Settings settings_;
	bool ui_enabled_;

	// Save feedback
	bool show_save_feedback_;
	std::string save_feedback_message_;
	float save_feedback_timer_;
	
	// File dialog state
	bool show_file_dialog_;
	bool deferred_open_config_; // Flag for deferred dialog opening
	FileDialogMode file_dialog_mode_;
	char file_path_buffer_[512];
	std::string file_path_; // Modern C++ string for file path
	std::string current_directory_;

	// Callbacks
	std::function<void(const std::string&)> open_config_callback_;
	std::function<void()> run_simulation_callback_;
	std::function<void()> rerun_simulation_callback_;
	std::function<void(const std::string&)> save_results_callback_;
	std::function<void()> direct_save_results_callback_;
	std::function<void()> reset_view_callback_;
	std::function<void(bool)> camera_mode_changed_callback_;
	std::vector<WorldText> world_text_queue_;
};
