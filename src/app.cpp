#include "app.hpp"

// Add includes for complete type definitions
#include "renderer/overlay.hpp"
#include "renderer/renderer.hpp"
#include "renderer/camera.hpp"
#include "renderer/settings.hpp"
#include "simulator/simulator.hpp"
#include "simulator/metrics.hpp"

#include "common/file_utils.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include "app.hpp"

#include "common/results_exporter.hpp"
#include "common/error_handler.hpp"
#include "common/result.hpp"
#include "common/error_types.hpp"

#include <cxxopts.hpp>

#include "simulator/logger.hpp"

#include "renderer/overlay.hpp"
#include "renderer/renderer.hpp"
#include "simulator/simulator.hpp"
#include "simulator/metrics.hpp"
#include "simulator/config.hpp"

// Static pointer for callbacks
static App* app_instance = nullptr;

// Static member definition
std::string App::executable_directory_;

// Static method implementations
void App::set_executable_path(const char* argv0) {
	std::filesystem::path exe_path(argv0);
	executable_directory_ = exe_path.parent_path().string();
	if (executable_directory_.empty()) {
		executable_directory_ = ".";
	}
}

std::string App::get_executable_directory() {
	return executable_directory_;
}

std::string App::get_output_path(const std::string& filename, const std::string& subdir) {
	std::filesystem::path base_dir = std::filesystem::path(executable_directory_);
	
	// Build path with optional subdirectory
	if (subdir.empty()) {
		base_dir = base_dir / "out";
	} else {
		base_dir = base_dir / subdir;
	}
	
	// Create the directory if it doesn't exist
	std::filesystem::create_directories(base_dir);
	
	return (base_dir / filename).string();
}

std::string App::get_results_path(const std::string& filename) {
	return get_output_path(filename, "results");
}

std::string App::get_debug_path(const std::string& filename) {
	return get_output_path(filename, "out");
}

App::App() : window_(nullptr), window_width_(1200), window_height_(800), should_close_(false),
              gui_mode_(true), force_csv_output_(false), first_simulation_run_(true), 
              left_mouse_pressed_(false), mouse_press_x_(0.0), mouse_press_y_(0.0) {
	app_instance = this;
	// Create shared metrics that will be used by all components
	shared_metrics_ = std::make_shared<Metrics>();
}

App::~App() {
	shutdown();
	app_instance = nullptr;
}

bool App::initialize(int argc, char* argv[]) {
	// Extract executable directory for output file paths
	if (argc > 0) {
		set_executable_path(argv[0]);
	}
	
	// Setup command line options
	cxxopts::Options options("Photron", "Photron - Monte Carlo Photon Transport Renderer");
	
	options.add_options()
		("c,config", "Configuration file path (optional for GUI mode)", cxxopts::value<std::string>())
		("headless", "Run simulation without GUI and generate output")
		("o,out", "Generate output files in GUI mode (first run only)")
		("l,log", "Enable debug logging messages")
		("h,help", "Show help message");

	// Allow positional arguments
	options.parse_positional({"config"});
	options.positional_help("[config_file]");

	try {
		auto result = options.parse(argc, argv);

		// Handle help
		if (result.count("help")) {
			std::cout << options.help() << std::endl;
			should_close_ = true;
			return true;
		}

		// Check for config file (now optional)
		bool has_config_file = result.count("config") > 0;
		if (has_config_file) {
			config_file_ = result["config"].as<std::string>();
		}
		bool log_mode = result.count("log") > 0;
		bool headless_mode = result.count("headless") > 0;
		bool force_csv = result.count("out") > 0;

		gui_mode_ = !headless_mode; // Set GUI mode (opposite of headless)
		force_csv_output_ = force_csv;
		first_simulation_run_ = true;

		// If headless mode is requested without a config file, that's an error
		if (headless_mode && !has_config_file) {
			std::cerr << "Error: Headless mode requires a configuration file" << std::endl;
			std::cerr << std::endl << options.help() << std::endl;
			return false;
		}

		// Always initialize Config (with defaults or from file)
		if (has_config_file) {
			if (!Config::initialize(config_file_)) {
				std::cerr << "Error: Failed to initialize configuration. Exiting." << std::endl;
				return false;
			}
		} else {
			Config::initialize(); // Initialize with defaults
		}

		// Set log mode if requested
		if (log_mode) {
			Config::get().set_log(true);
		}

		// Configure ErrorHandler with logging state
		ErrorHandler::instance().set_logging_enabled(Config::get().log());

		// Initialize simulator and run simulation only if we have a config file
		if (has_config_file) {
			// Initialize simulator first (always, regardless of headless mode)
			if (!initialize_simulator()) {
				std::cerr << "Failed to initialize simulator" << std::endl;
				return false;
			}

			// Run simulation with progress feedback
			run_simulation_with_progress();

			// If headless mode, we're done - no GUI needed
			if (headless_mode) {
				should_close_ = true;
				return true;
			}
		}

		// Initialize GUI (either after simulation or directly if no config file)
			if (!initialize_gui()) {
				ErrorHandler::instance().report_error("Failed to initialize GUI");
				return false;
			}		// If no config file was provided, auto-open the config dialog
		if (!has_config_file) {
			overlay_->open_config_dialog();
		}

		return true;
	}
	catch (const cxxopts::exceptions::exception& e) {
		std::cerr << "Error parsing command line: " << e.what() << std::endl;
		std::cerr << std::endl << options.help() << std::endl;
		return false;
	}
}

void App::run() {
	// In headless mode, just exit immediately
	if (!window_) {
		return;
	}
	
	while (!glfwWindowShouldClose(window_) && !should_close_) {
		glfwPollEvents();
		update();
		render();
		glfwSwapBuffers(window_);
	}
}

void App::shutdown() {
	if (overlay_) {
		overlay_->shutdown();
		overlay_.reset();
	}

	if (renderer_) {
		renderer_.reset();
	}

	if (simulator_) {
		simulator_.reset();
	}
	
	// Reset config service
	Config::shutdown();

	if (window_) {
		glfwDestroyWindow(window_);
		window_ = nullptr;
		glfwTerminate();
	}
}

void App::setup_callbacks() {
	glfwSetFramebufferSizeCallback(window_, framebuffer_size_callback);
	glfwSetKeyCallback(window_, key_callback);
	glfwSetCursorPosCallback(window_, cursor_pos_callback);
	glfwSetMouseButtonCallback(window_, mouse_button_callback);
	glfwSetScrollCallback(window_, scroll_callback);
}

void App::setup_overlay_callbacks() {
	// Set up file operation callbacks
	overlay_->set_open_config_callback([this](const std::string& filepath) {
		// Wrap entire callback in try-catch to prevent crashes
		try {
			// Disable UI while loading
			overlay_->set_ui_enabled(false);

			// Update config file path
			config_file_ = filepath;
			
			// Reinitialize Config with the new file
			Config::shutdown();
			if (!Config::initialize(config_file_)) {
				// Config parsing failed - reset simulator state and re-enable UI
				simulator_.reset(); // Reset simulator so GUI shows "No configuration loaded"
				ErrorHandler::instance().report_error("Failed to load configuration file: " + filepath);
				// Auto-open the config dialog so user can try loading another file (deferred to next frame)
				overlay_->open_config_dialog_deferred();
				overlay_->set_ui_enabled(true);
				return;
			}

			// Configure ErrorHandler with new logging state
			ErrorHandler::instance().set_logging_enabled(Config::get().log());

			// Create simulator if it doesn't exist
			if (!simulator_) {
				// Create simulator without initializing it yet
				simulator_ = std::make_unique<Simulator>();
				// Provide shared metrics to simulator
				simulator_->set_shared_metrics(shared_metrics_);
			}

			// Load new configuration into simulator (only call this once)
			auto init_result = simulator_->initialize(filepath.c_str());
			if (!init_result.is_ok()) {
				ErrorHandler::instance().report_error("Failed to initialize simulator: " + ErrorMessage::format(init_result.error()));
				overlay_->set_ui_enabled(true);
				return;
			}

		// Run simulation immediately
		auto sim_result = simulator_->simulate();
		if (!sim_result.is_ok()) {
			ErrorHandler::instance().report_error("Simulation failed: " + ErrorMessage::format(sim_result.error()));
			overlay_->set_ui_enabled(true);
			return;
		}
		bool generate_csv = !gui_mode_ || (force_csv_output_ && first_simulation_run_);
		simulator_->report(generate_csv);
		if (first_simulation_run_) first_simulation_run_ = false; // Mark first run as complete

		// Reset camera to default position when new config is loaded
		if (renderer_) {
			// Always switch to Orbit mode before resetting to ensure proper camera direction reset
			renderer_->set_camera_mode(true); // true = Orbit mode
			renderer_->reset_camera();

			// Update the UI to reflect the mode change
			if (overlay_) {
				Settings& settings = overlay_->get_settings();
				settings.camera_mode = CameraMode::Orbit;
			}
		}

		// Re-enable UI
		overlay_->set_ui_enabled(true);
		
		} catch (const std::exception& e) {
			// Catch any unhandled exceptions to prevent crashes
			// Use direct cerr to avoid potential recursive exceptions in ErrorHandler
			std::cerr << "Exception in config loading callback: " << e.what() << std::endl;
			overlay_->set_ui_enabled(true);
		} catch (...) {
			// Catch any other unknown exceptions
			std::cerr << "Unknown exception in config loading callback" << std::endl;
			overlay_->set_ui_enabled(true);
		}
	});

	overlay_->set_run_simulation_callback([this]() {
		// Disable UI while running
		overlay_->set_ui_enabled(false);

		// Run a single additional photon instead of the full simulation
		simulator_->simulate_single_photon();

		// Invalidate all caches since simulation data has changed
		if (renderer_) {
			renderer_->invalidate_all_caches();
		}

		// Re-enable UI
		overlay_->set_ui_enabled(true);
	});

	overlay_->set_rerun_simulation_callback([this]() {
		// Disable UI while running
		overlay_->set_ui_enabled(false);

		// Clear previous results and rerun simulation cleanly
		std::cout << "Rerunning simulation (clearing previous results)" << std::endl;
		if (Config::get().log()) {
			Logger::instance().log_info("Rerunning simulation (clearing previous results)");
		}

		// Re-initialize with the same config file to clear data
		if (!config_file_.empty()) {
			auto init_result = simulator_->initialize(config_file_.c_str());
			if (!init_result.is_ok()) {
				ErrorHandler::instance().report_error("Failed to reinitialize simulator: " + ErrorMessage::format(init_result.error()));
				overlay_->set_ui_enabled(true); // Re-enable UI on error
				return;
			}
			
			auto sim_result = simulator_->simulate();
			if (!sim_result.is_ok()) {
				ErrorHandler::instance().report_error("Rerun simulation failed: " + ErrorMessage::format(sim_result.error()));
				overlay_->set_ui_enabled(true); // Re-enable UI on error
				return;
			}
			
			simulator_->report(!gui_mode_); // Only generate CSV in headless mode (no CSV for GUI restarts)
		}
		else {
			std::cout << "No config file loaded. Please load a configuration first." << std::endl;
		}

		// Invalidate all caches since simulation data has been reset
		if (renderer_) {
			renderer_->invalidate_all_caches();
		}

		// Re-enable UI
		overlay_->set_ui_enabled(true);
	});

	overlay_->set_save_results_callback([this](const std::string& filepath) {
		// Implement JSON results export
		if (filepath.find(".json") != std::string::npos) {
			ResultsExporter::instance().export_json(*simulator_, filepath);
		}
		else {
			ResultsExporter::instance().export_text(*simulator_, filepath);
		}
	});
	
	overlay_->set_direct_save_results_callback([this]() {
		// Direct save to out directory
		if (simulator_) {
			// Show save progress feedback immediately with initial message
			overlay_->show_save_feedback("Saving simulation results...");
			
			// Small delay to ensure the progress bar is visible
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			
			simulator_->report(true); // Always generate all files (CSV + log) when user explicitly saves
			
			// Update the message after save is complete
			overlay_->show_save_feedback("Simulation results saved to out/ directory");
		}
	});

	overlay_->set_reset_view_callback([this]() {
		if (renderer_) {
			// Always switch to Orbit mode before resetting to ensure proper camera direction reset
			renderer_->set_camera_mode(true); // true = Orbit mode
			renderer_->reset_camera();

			// Update the UI to reflect the mode change
			if (overlay_) {
				Settings& settings = overlay_->get_settings();
				settings.camera_mode = CameraMode::Orbit;
			}
		}
	});

	overlay_->set_camera_mode_changed_callback([this](bool is_arc_mode) {
		if (renderer_) {
			renderer_->set_camera_mode(is_arc_mode);
		}
	});

	// Set up callback for when renderer changes camera mode (e.g., WASD auto-switch)
	if (renderer_) {
		renderer_->set_camera_mode_change_callback([this](bool is_arc_mode) {
			if (overlay_) {
				Settings& settings = overlay_->get_settings();
				settings.camera_mode = is_arc_mode ? CameraMode::Orbit : CameraMode::Free;
			}
		});

		// Set up text rendering callback for energy labels
		renderer_->set_text_render_callback([this](const std::string& text, float x, float y, const glm::vec4& color) {
			if (overlay_) {
				overlay_->add_world_text(text, x, y, color);
			}
		});
	}
}

void App::update() {
	// Update camera and other systems
	if (renderer_) {
		renderer_->update();

		// Handle cursor visibility after all updates
		bool should_hide_cursor = renderer_->should_capture_mouse();
		glfwSetInputMode(window_, GLFW_CURSOR, should_hide_cursor ? GLFW_CURSOR_HIDDEN : GLFW_CURSOR_NORMAL);
	}
}

void App::render() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Clear previous frame's world text before rendering new frame
	if (overlay_) {
		overlay_->clear_world_text();
	}

	// Update renderer settings from overlay before rendering
	if (renderer_ && overlay_) {
		Settings& settings = overlay_->get_settings();
		renderer_->auto_manage_energy_labels(settings);
		renderer_->set_settings(settings);
	}

	// Render 3D scene FIRST
	if (renderer_ && simulator_) {
		renderer_->render(*simulator_);
	}

	// Render ImGui overlay AFTER 3D scene (so it appears on top)
	if (overlay_) {
		// Reset OpenGL state for ImGui
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		if (simulator_) {
			overlay_->render_with_simulator(*simulator_);
		} else {
			overlay_->render();
		}

		// Restore depth testing for next frame
		glEnable(GL_DEPTH_TEST);
	}
}

// GLFW callbacks
void App::framebuffer_size_callback(GLFWwindow*, int width, int height) {
	if (app_instance) {
		app_instance->window_width_ = width;
		app_instance->window_height_ = height;
		glViewport(0, 0, width, height);

		if (app_instance->renderer_) {
			app_instance->renderer_->set_viewport(width, height);
		}
	}
}

void App::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS) {
		switch (key) {
			case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
		}
	}

	// Handle camera input for WASD only when file dialog is not open
	// This allows smooth movement while interacting with UI elements but disables when dialog is open
	if (app_instance && app_instance->renderer_ && app_instance->overlay_) {
		// Don't handle camera input if file dialog is open
		if (!app_instance->overlay_->is_file_dialog_open()) {
			app_instance->renderer_->handle_key_input(key, scancode, action, mods);
		}
	}
}

void App::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos) {
	// Don't handle camera input if ImGui wants the mouse
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (app_instance && app_instance->renderer_ && app_instance->overlay_) {
		// Don't handle camera input if file dialog is open
		if (!app_instance->overlay_->is_file_dialog_open()) {
			
			// Check if we're dragging in free camera mode - reset and switch to orbit immediately
			if (app_instance->left_mouse_pressed_ && !app_instance->renderer_->is_arc_camera_mode()) {
				double dx = xpos - app_instance->mouse_press_x_;
				double dy = ypos - app_instance->mouse_press_y_;
				double distance = std::sqrt(dx * dx + dy * dy);
				
				// If we've dragged beyond threshold, reset to orbit and continue dragging
				if (distance >= app_instance->DRAG_THRESHOLD) {
					app_instance->renderer_->set_camera_mode(true); // true = Orbit mode
					app_instance->renderer_->reset_camera();
					
					// Update overlay settings to reflect the change
					auto& settings = app_instance->overlay_->get_settings();
					settings.camera_mode = CameraMode::Orbit;
					
					// Clear the pressed flag so we don't trigger this again
					app_instance->left_mouse_pressed_ = false;
					
					// Continue processing the mouse movement for orbit camera
				}
			}
			
			app_instance->renderer_->handle_mouse_move(static_cast<float>(xpos), static_cast<float>(ypos));

			// Handle mouse capture for FPS mode (but don't set cursor here - done in update loop)
			bool should_capture = app_instance->renderer_->should_capture_mouse();
			if (should_capture) {
				int width, height;
				glfwGetWindowSize(window, &width, &height);
				float center_x = width / 2.0f;
				float center_y = height / 2.0f;

				// Set the center position for proper mouse handling
				if (app_instance->renderer_) {
					app_instance->renderer_->get_camera().set_mouse_center(center_x, center_y);
				}

				// Center the cursor for FPS mode
				glfwSetCursorPos(window, center_x, center_y);
			}
		}
	}
}

void App::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	// Don't handle camera input if ImGui wants the mouse
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (app_instance && app_instance->renderer_ && app_instance->overlay_) {
		// Don't handle camera input if file dialog is open
		if (!app_instance->overlay_->is_file_dialog_open()) {
			
			if (button == GLFW_MOUSE_BUTTON_LEFT) {
				if (action == GLFW_PRESS) {
					// Record mouse press position for drag detection
					double xpos, ypos;
					glfwGetCursorPos(window, &xpos, &ypos);
					app_instance->left_mouse_pressed_ = true;
					app_instance->mouse_press_x_ = xpos;
					app_instance->mouse_press_y_ = ypos;
				}
				else if (action == GLFW_RELEASE) {
					// Clear the pressed flag on release
					app_instance->left_mouse_pressed_ = false;
				}
			}
			
			// Always pass through to normal mouse handling
			app_instance->renderer_->handle_mouse_button(button, action, mods);
		}
	}
}

void App::scroll_callback(GLFWwindow*, double xoffset, double yoffset) {
	// Don't handle camera input if ImGui wants the mouse
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (app_instance && app_instance->renderer_) {
		app_instance->renderer_->handle_mouse_scroll(static_cast<float>(xoffset), static_cast<float>(yoffset));
	}
}


bool App::initialize_simulator() {
	// Initialize simulator (Config should already be initialized by this point)
	simulator_ = std::make_unique<Simulator>();
	// Provide shared metrics to simulator
	simulator_->set_shared_metrics(shared_metrics_);
	
	// No progress callbacks needed - simulation runs optimally without GUI interference
	
	auto result = simulator_->initialize(config_file_.c_str());
	if (!result.is_ok()) {
		// Use ErrorHandler to report the structured error
		ErrorHandler::instance().report_error("Simulator initialization failed: " + ErrorMessage::format(result.error()));
		return false;
	}
	
	return true;
}

bool App::initialize_gui() {
	// Initialize GLFW
	if (!glfwInit()) {
		std::cerr << "Failed to initialize GLFW" << std::endl;
		return false;
	}

	// Set OpenGL version and profile
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Create window
	window_ = glfwCreateWindow(window_width_, window_height_, "Photron - Monte Carlo Photon Transport", nullptr, nullptr);
	if (!window_) {
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return false;
	}

	// Center window on screen
	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	if (mode) {
		int center_x = (mode->width - window_width_) / 2;
		int center_y = (mode->height - window_height_) / 2;
		glfwSetWindowPos(window_, center_x, center_y);
	}

	glfwMakeContextCurrent(window_);
	glfwSwapInterval(1); // Enable vsync

	// Initialize GLEW
	if (glewInit() != GLEW_OK) {
		std::cerr << "Failed to initialize GLEW" << std::endl;
		glfwTerminate();
		return false;
	}

	// Setup callbacks
	setup_callbacks();

	// Initialize renderer
	renderer_ = std::make_unique<Renderer>();
	if (!renderer_->initialize()) {
		std::cerr << "Failed to initialize renderer" << std::endl;
		return false;
	}

	// Set correct initial viewport size (framebuffer callback only triggers on resize)
	renderer_->set_viewport(window_width_, window_height_);

	// Initialize overlay
	overlay_ = std::make_unique<Overlay>();
	if (!overlay_->initialize(window_)) {
		std::cerr << "Failed to initialize overlay" << std::endl;
		return false;
	}

	// Set up overlay callbacks
	setup_overlay_callbacks();
	
	return true;
}

void App::run_simulation_with_progress() {
	if (config_file_.empty()) {
		ErrorHandler::instance().report_error(ErrorMessage::format(ConfigError::FileNotFound, "run_simulation_with_progress called without config file"));
		return;
	}
	
	if (!simulator_) {
		ErrorHandler::instance().report_error(ErrorMessage::format(SimulationError::InitializationFailed, "run_simulation_with_progress called without simulator"));
		return;
	}

	std::cout << "=== Starting Monte Carlo Simulation ===" << std::endl;
	std::cout << "Configuration: " << config_file_ << std::endl;
	std::cout << "Number of photons: " << Config::get().num_photons() << std::endl;
	std::cout << std::endl;
	
	// NOTE: Don't show simulation progress overlay here - this is for headless mode
	// Progress bar is only shown when simulation is triggered from GUI
	
	if (Config::get().log()) {
		Logger::instance().log_info("=== Starting Monte Carlo Simulation ===");
		Logger::instance().log_info("Configuration: " + config_file_);
		Logger::instance().log_info("Number of photons: " + std::to_string(Config::get().num_photons()));
	}
	
	// Run simulation
	auto sim_result = simulator_->simulate();
	if (!sim_result.is_ok()) {
		ErrorHandler::instance().report_error("Simulation failed: " + ErrorMessage::format(sim_result.error()));
		return;
	}
	
	bool generate_csv = !gui_mode_ || (force_csv_output_ && first_simulation_run_);
	simulator_->report(generate_csv);
	if (first_simulation_run_) first_simulation_run_ = false; // Mark first run as complete
	
	// NOTE: Don't hide simulation progress overlay here - it wasn't shown
	
	if (Config::get().log()) {
		Logger::instance().log_info("=== Simulation Complete ===");
		Logger::instance().log_info("Initializing GUI...");
	}
}
