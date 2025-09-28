#pragma once

#include <memory>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

class Overlay;
class Renderer;
class Simulator;
class Metrics;

class App
{
public:
	App();
	~App();

	bool initialize(int argc, char* argv[]);
	void run();
	
	// Access to shared metrics for external components
	std::shared_ptr<Metrics> get_shared_metrics() const { return shared_metrics_; }
	void shutdown();

	// Getters
	GLFWwindow* get_window() const { return window_; }
	int get_window_width() const { return window_width_; }
	int get_window_height() const { return window_height_; }
	Simulator* get_simulator() const { return simulator_.get(); }

	// Static method to get executable directory for output files
	static std::string get_executable_directory();
	static void set_executable_path(const char* argv0);
	
	// Unified output path generation (replaces multiple scattered methods)
	static std::string get_output_path(const std::string& filename, const std::string& subdir = "");
	static std::string get_results_path(const std::string& filename);
	static std::string get_debug_path(const std::string& filename);

private:
	void setup_callbacks();
	void setup_overlay_callbacks();
	void update();
	void render();
	
	// Separate initialization methods
	bool initialize_simulator();
	bool initialize_gui();
	void run_simulation_with_progress();

	// Results saving methods
	void save_results_as_json(const std::string& filepath);
	void save_results_as_text(const std::string& filepath);

	// GLFW callbacks
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos);
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

	GLFWwindow* window_;
	int window_width_;
	int window_height_;

	std::unique_ptr<Simulator> simulator_;
	std::unique_ptr<Renderer> renderer_;
	std::unique_ptr<Overlay> overlay_;
	std::shared_ptr<Metrics> shared_metrics_;  // Shared metrics for all components

	bool should_close_;
	std::string config_file_;
	bool gui_mode_; // Track if running in GUI mode (vs headless)
	bool force_csv_output_; // Force CSV output even in GUI mode (from command line)
	bool first_simulation_run_; // Track if this is the first simulation run
	
	// Mouse tracking for click vs drag detection
	bool left_mouse_pressed_;
	double mouse_press_x_, mouse_press_y_;
	static constexpr double DRAG_THRESHOLD = 5.0; // pixels
	
private:
	static std::string executable_directory_;
};
