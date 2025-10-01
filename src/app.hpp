/**
 * @file app.hpp
 * @brief Main application class that coordinates between simulator and renderer
 *
 * The App class serves as the central controller for the Photron application,
 * managing initialization, the main loop, and coordination between the Monte Carlo
 * simulation engine and the OpenGL-based 3D visualization renderer.
 */

#pragma once

#include <memory>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// Forward declarations
class Overlay;
class Renderer;
class Simulator;
class Metrics;

/**
 * @class App
 * @brief Central application controller for Photron Monte Carlo photon transport simulation
 *
 * The App class manages the complete application lifecycle including:
 * - Command line parsing and configuration loading
 * - OpenGL/GLFW window management
 * - Coordination between simulation and rendering subsystems
 * - GUI overlay management with ImGui
 * - File I/O operations for results and configuration
 * - Event handling and user input processing
 *
 * The application supports both interactive GUI mode and headless batch processing.
 * In GUI mode, users can load configurations, run simulations, and visualize results
 * in real-time. In headless mode, simulations run without graphics for batch processing.
 */
class App
{
public:
	/**
	 * @brief Construct a new App object with default settings
	 *
	 * Initializes window dimensions, creates shared metrics object for inter-component
	 * communication, and sets up default application state.
	 */
	App();

	/**
	 * @brief Destroy the App object and clean up resources
	 *
	 * Ensures proper shutdown sequence and resource deallocation.
	 */
	~App();

	/**
	 * @brief Initialize the application with command line arguments
	 *
	 * Parses command line options, loads configuration files, initializes OpenGL
	 * context, and sets up all subsystems (simulator, renderer, overlay).
	 *
	 * @param argc Number of command line arguments
	 * @param argv Array of command line argument strings
	 * @return true if initialization succeeded, false otherwise
	 */
	bool initialize(int argc, char* argv[]);

	/**
	 * @brief Run the main application loop
	 *
	 * In GUI mode, runs the GLFW event loop with rendering and user interaction.
	 * In headless mode, runs simulation directly and exits.
	 */
	void run();

	/**
	 * @brief Get shared metrics object for inter-component communication
	 *
	 * The shared metrics object allows the simulator and renderer to communicate
	 * simulation progress and results without tight coupling.
	 *
	 * @return std::shared_ptr<Metrics> Shared metrics instance
	 */
	std::shared_ptr<Metrics> get_shared_metrics() const { return shared_metrics_; }

	/**
	 * @brief Clean shutdown of all application subsystems
	 *
	 * Properly destroys overlay, renderer, simulator, resets configuration,
	 * and terminates GLFW context.
	 */
	void shutdown();

	// Getters

	/**
	 * @brief Get the GLFW window handle
	 * @return GLFWwindow* Window handle or nullptr if not initialized
	 */
	GLFWwindow* get_window() const { return window_; }

	/**
	 * @brief Get current window width in pixels
	 * @return int Window width
	 */
	int get_window_width() const { return window_width_; }

	/**
	 * @brief Get current window height in pixels
	 * @return int Window height
	 */
	int get_window_height() const { return window_height_; }

	/**
	 * @brief Get pointer to the simulator instance
	 * @return Simulator* Simulator instance or nullptr if not initialized
	 */
	Simulator* get_simulator() const { return simulator_.get(); }

	// Static utility methods for file path management

	/**
	 * @brief Get the directory containing the executable
	 *
	 * Used to construct relative paths for shader files, configuration files,
	 * and output directories. Essential for portable file access.
	 *
	 * @return std::string Absolute path to executable directory
	 */
	static std::string get_executable_directory();

	/**
	 * @brief Set executable path from argv[0] for path resolution
	 *
	 * Must be called during initialization to enable proper file path resolution.
	 *
	 * @param argv0 The argv[0] parameter from main()
	 */
	static void set_executable_path(const char* argv0);

	/**
	 * @brief Generate unified output path for files with optional subdirectory
	 *
	 * Creates paths relative to executable directory with consistent structure.
	 * Automatically creates subdirectories if they don't exist.
	 *
	 * @param filename Name of the output file
	 * @param subdir Optional subdirectory (default: root output directory)
	 * @return std::string Complete file path for output
	 */
	static std::string get_output_path(const std::string& filename, const std::string& subdir = "");

	/**
	 * @brief Generate path for simulation results files
	 *
	 * Convenience wrapper for results-specific output directory.
	 *
	 * @param filename Name of the results file
	 * @return std::string Complete file path in results directory
	 */
	static std::string get_results_path(const std::string& filename);

	/**
	 * @brief Generate path for debug output files
	 *
	 * Convenience wrapper for debug-specific output directory.
	 *
	 * @param filename Name of the debug file
	 * @return std::string Complete file path in debug directory
	 */
	static std::string get_debug_path(const std::string& filename);

private:
	// Initialization methods

	/**
	 * @brief Set up GLFW event callbacks for user input handling
	 */
	void setup_callbacks();

	/**
	 * @brief Configure overlay callbacks for GUI interactions
	 *
	 * Sets up callbacks for file operations, simulation control,
	 * and camera management through the ImGui interface.
	 */
	void setup_overlay_callbacks();

	/**
	 * @brief Initialize the Monte Carlo simulator subsystem
	 *
	 * Creates simulator instance, configures shared metrics, and loads
	 * configuration file if provided.
	 *
	 * @return true if simulator initialization succeeded
	 */
	bool initialize_simulator();

	/**
	 * @brief Initialize GUI subsystem (renderer and overlay)
	 *
	 * Sets up OpenGL context, creates renderer and ImGui overlay,
	 * configures callbacks and initial state.
	 *
	 * @return true if GUI initialization succeeded
	 */
	bool initialize_gui();

	// Main loop methods

	/**
	 * @brief Update application state each frame
	 *
	 * Handles per-frame logic updates before rendering.
	 */
	void update();

	/**
	 * @brief Render the current frame
	 *
	 * Coordinates 3D scene rendering and ImGui overlay rendering
	 * with proper OpenGL state management.
	 */
	void render();

	// Simulation control

	/**
	 * @brief Run simulation with progress tracking for headless mode
	 *
	 * Executes Monte Carlo simulation without GUI progress overlay,
	 * with console output for progress monitoring.
	 */
	void run_simulation_with_progress();

	// GLFW event callbacks

	/**
	 * @brief Handle window framebuffer resize events
	 * @param window GLFW window handle
	 * @param width New framebuffer width
	 * @param height New framebuffer height
	 */
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height);

	/**
	 * @brief Handle keyboard input events
	 * @param window GLFW window handle
	 * @param key GLFW key code
	 * @param scancode Platform-specific scan code
	 * @param action GLFW action (press, release, repeat)
	 * @param mods Modifier key flags
	 */
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

	/**
	 * @brief Handle cursor position events for camera control
	 * @param window GLFW window handle
	 * @param xpos Cursor X coordinate
	 * @param ypos Cursor Y coordinate
	 */
	static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos);

	/**
	 * @brief Handle mouse button events for camera and UI interaction
	 * @param window GLFW window handle
	 * @param button Mouse button code
	 * @param action GLFW action (press or release)
	 * @param mods Modifier key flags
	 */
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

	/**
	 * @brief Handle mouse scroll events for camera zoom
	 * @param window GLFW window handle
	 * @param xoffset Horizontal scroll offset
	 * @param yoffset Vertical scroll offset
	 */
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

	// Core subsystem instances
	GLFWwindow* window_;                      		///< GLFW window handle
	int window_width_;                        		///< Current window width in pixels
	int window_height_;                       		///< Current window height in pixels

	std::unique_ptr<Simulator> simulator_;    		///< Monte Carlo photon transport engine
	std::unique_ptr<Renderer> renderer_;      		///< OpenGL 3D visualization renderer
	std::unique_ptr<Overlay> overlay_;        		///< ImGui interface overlay
	std::shared_ptr<Metrics> shared_metrics_; 		///< Shared simulation metrics for inter-component communication

	// Application state
	bool should_close_;         					///< Flag to terminate main loop
	std::string config_file_;   					///< Path to loaded configuration file
	bool gui_mode_;             					///< True for GUI mode, false for headless
	bool force_csv_output_;     					///< Force CSV output even in GUI mode (from command line)
	bool first_simulation_run_; 					///< Track if this is the first simulation run for output control

	// Mouse interaction state
	bool left_mouse_pressed_;                     	///< Left mouse button currently pressed
	double mouse_press_x_, mouse_press_y_;        	///< Mouse position when button was pressed
	static constexpr double DRAG_THRESHOLD = 5.0; 	///< Minimum pixel distance to distinguish click from drag

	// Static file path management
	static std::string executable_directory_; 		///< Cached executable directory path
};
