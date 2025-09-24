#pragma once

#include <memory>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

class Overlay;
class Renderer;
class Simulator;

class App
{
public:
	App();
	~App();

	bool initialize(int argc, char* argv[]);
	void run();
	void shutdown();

	// Getters
	GLFWwindow* get_window() const { return window_; }
	int get_window_width() const { return window_width_; }
	int get_window_height() const { return window_height_; }
	Simulator* get_simulator() const { return simulator_.get(); }

	// Static method to get executable directory for output files
	static std::string get_executable_directory();
	static void set_executable_path(const char* argv0);
	static std::string get_output_path(const std::string& filename);

private:
	void setup_callbacks();
	void setup_overlay_callbacks();
	void update();
	void render();
	
	// Separate initialization methods
	bool initialize_simulator(bool verbose_mode);
	bool initialize_gui();
	void run_simulation_with_progress();

	// Results saving methods
	void save_results_as_json(const std::string& filepath);
	void save_results_as_text(const std::string& filepath);

	// Helper methods to aggregate metrics from all mediums
	double aggregate_path_length() const;
	double aggregate_scatter_events() const;
	double aggregate_average_step_size() const;
	double aggregate_diffusion_distance() const;
	
	// Helper methods to aggregate energy conservation values from all mediums
	double aggregate_total_absorption() const;
	double aggregate_total_reflection() const;
	double aggregate_total_transmission() const;
	double aggregate_total_diffusion() const;
	double aggregate_surface_reflection() const;
	double aggregate_surface_refraction() const;

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

	bool should_close_;
	std::string config_file_;
	
private:
	static std::string executable_directory_;
};
