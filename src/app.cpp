#include "app.hpp"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui.h>

#include "renderer/overlay.hpp"
#include "renderer/renderer.hpp"
#include "simulator/simulator.hpp"

// Static pointer for callbacks
static App* app_instance = nullptr;

App::App() : window_(nullptr), window_width_(1200), window_height_(800), should_close_(false) {
	app_instance = this;
}

App::~App() {
	shutdown();
	app_instance = nullptr;
}

bool App::initialize(int argc, char* argv[]) {
	// Check command-line arguments
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <config_file> [--verbose]" << std::endl;
		std::cerr << "  --verbose  : Enable verbose initialization messages" << std::endl;
		return false;
	}

	config_file_ = argv[1];
	
	// Check for verbose flag
	bool verbose_mode = false;
	for (int i = 2; i < argc; ++i) {
		if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v") {
			verbose_mode = true;
			break;
		}
	}

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
	window_ =
		glfwCreateWindow(window_width_, window_height_, "Photron - Monte Carlo Photon Transport", nullptr, nullptr);
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

	std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;

	// Setup callbacks
	setup_callbacks();

	// Initialize simulator
	simulator_ = std::make_unique<Simulator>();
	
	// Set verbose mode if requested
	if (verbose_mode) {
		simulator_->config.set_verbose(true);
	}
	
	if (!simulator_->initialize(config_file_.c_str())) {
		std::cerr << "Failed to initialize simulator" << std::endl;
		return false;
	}

	// Run simulation
	std::cout << "Running photon transport simulation..." << std::endl;
	simulator_->simulate();
	simulator_->report();

	// Stop if rendering is disabled
	if (argc > 2 && std::string(argv[2]) == "headless") {
		should_close_ = true;
		return true;
	}

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

void App::run() {
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

	if (window_) {
		glfwDestroyWindow(window_);
		window_ = nullptr;
	}

	glfwTerminate();
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
		// Disable UI while loading
		overlay_->set_ui_enabled(false);

		// Load new configuration
		simulator_->initialize(filepath.c_str());
		config_file_ = filepath;

		// Run simulation immediately
		std::cout << "Running simulation with config: " << filepath << std::endl;
		simulator_->simulate();
		simulator_->report();

		// Re-enable UI
		overlay_->set_ui_enabled(true);
	});

	overlay_->set_run_simulation_callback([this]() {
		// Disable UI while running
		overlay_->set_ui_enabled(false);

		// Run a single additional photon instead of the full simulation
		simulator_->simulate_single_photon();

		// Invalidate energy label cache since simulation data has changed
		if (renderer_) {
			renderer_->invalidate_energy_label_cache();
		}

		// Re-enable UI
		overlay_->set_ui_enabled(true);
	});

	overlay_->set_rerun_simulation_callback([this]() {
		// Disable UI while running
		overlay_->set_ui_enabled(false);

		// Clear previous results and rerun simulation cleanly
		std::cout << "Rerunning simulation (clearing previous results)..." << std::endl;

		// Re-initialize with the same config file to clear data
		if (!config_file_.empty()) {
			simulator_->initialize(config_file_.c_str());
			simulator_->simulate();
			simulator_->report();
		}
		else {
			std::cout << "No config file loaded. Please load a configuration first." << std::endl;
		}

		// Invalidate energy label cache since simulation data has been reset
		if (renderer_) {
			renderer_->invalidate_energy_label_cache();
		}

		// Re-enable UI
		overlay_->set_ui_enabled(true);
	});

	overlay_->set_save_results_callback([this](const std::string& filepath) {
		std::cout << "Saving simulation results to: " << filepath << std::endl;

		// Implement JSON results export
		if (filepath.find(".json") != std::string::npos) {
			save_results_as_json(filepath);
		}
		else {
			save_results_as_text(filepath);
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

		overlay_->render_with_simulator(*simulator_);

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

	// Always handle camera input for WASD even when ImGui has keyboard focus
	// This allows smooth movement while interacting with UI elements
	if (app_instance && app_instance->renderer_) {
		app_instance->renderer_->handle_key_input(key, scancode, action, mods);
	}
}

void App::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos) {
	// Don't handle camera input if ImGui wants the mouse
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (app_instance && app_instance->renderer_) {
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

void App::mouse_button_callback(GLFWwindow*, int button, int action, int mods) {
	// Don't handle camera input if ImGui wants the mouse
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (app_instance && app_instance->renderer_) {
		app_instance->renderer_->handle_mouse_button(button, action, mods);
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

void App::save_results_as_json(const std::string& filepath) {
	if (!simulator_) {
		std::cerr << "No simulator available for saving results." << std::endl;
		return;
	}

	std::ofstream file(filepath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing: " << filepath << std::endl;
		return;
	}

	// Save simulation results as JSON
	file << "{\n";
	file << "  \"config_file\": \"" << config_file_ << "\",\n";
	file << "  \"simulation_parameters\": {\n";
	file << "    \"num_photons\": " << simulator_->config.num_photons() << ",\n";
	file << "    \"voxel_size\": " << simulator_->config.vox_size() << ",\n";
	file << "    \"grid_size\": [" << simulator_->config.nx() << ", " << simulator_->config.ny() << ", "
		 << simulator_->config.nz() << "]\n";
	file << "  },\n";
	file << "  \"results\": {\n";
	file << "    \"path_length\": " << simulator_->metrics.get_path_length() << ",\n";
	file << "    \"scatter_events\": " << simulator_->metrics.get_scatter_events() << ",\n";
	file << "    \"average_step_size\": " << simulator_->metrics.get_average_step_size() << ",\n";
	file << "    \"diffusion_distance\": " << simulator_->metrics.get_diffusion_distance() << ",\n";
	file << "    \"total_absorption\": " << simulator_->metrics.get_total_absorption() << ",\n";
	file << "    \"total_reflection\": " << simulator_->metrics.get_total_reflection() << ",\n";
	file << "    \"total_transmission\": " << simulator_->metrics.get_total_transmission() << ",\n";
	file << "    \"total_diffusion\": " << simulator_->metrics.get_total_diffusion() << ",\n";
	file << "    \"surface_reflection\": " << simulator_->metrics.get_surface_reflection() << ",\n";
	file << "    \"surface_refraction\": " << simulator_->metrics.get_surface_refraction() << "\n";
	file << "  },\n";
	file << "  \"tissue_properties\": [\n";

	auto tissues = simulator_->get_all_tissues();
	for (size_t i = 0; i < tissues.size(); ++i) {
		const auto& tissue = tissues[i];
		file << "    {\n";
		file << "      \"id\": " << tissue.id << ",\n";
		file << "      \"eta\": " << tissue.eta << ",\n";
		file << "      \"mua\": " << tissue.mu_a << ",\n";
		file << "      \"mus\": " << tissue.mu_s << ",\n";
		file << "      \"ani\": " << tissue.g << "\n";
		file << "    }" << (i < tissues.size() - 1 ? "," : "") << "\n";
	}

	file << "  ]\n";
	file << "}\n";

	file.close();
	std::cout << "Results saved successfully to " << filepath << std::endl;
}

void App::save_results_as_text(const std::string& filepath) {
	if (!simulator_) {
		std::cerr << "No simulator available for saving results." << std::endl;
		return;
	}

	std::ofstream file(filepath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing: " << filepath << std::endl;
		return;
	}

	// Save simulation results as formatted text
	file << "=== Photron Simulation Results ===\n\n";
	file << "Configuration File: " << config_file_ << "\n";
	file << "Timestamp: " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << "\n\n";

	file << "Simulation Parameters:\n";
	file << "  Number of Photons: " << simulator_->config.num_photons() << "\n";
	file << "  Voxel Size: " << simulator_->config.vox_size() << "\n";
	file << "  Grid Dimensions: " << simulator_->config.nx() << " x " << simulator_->config.ny() << " x "
		 << simulator_->config.nz() << "\n";
	file << "  Total Voxels: " << simulator_->get_total_voxel_count() << "\n\n";

	if (simulator_->config.num_photons() > 1) {
		file << "Multi-Photon Statistics:\n";
		double total_weight = simulator_->metrics.get_total_absorption() + simulator_->metrics.get_total_reflection()
							  + simulator_->metrics.get_total_transmission();
		if (total_weight > 0) {
			file << "  Absorption Probability: " << (simulator_->metrics.get_total_absorption() / total_weight) * 100.0
				 << "%\n";
			file << "  Reflection Probability: " << (simulator_->metrics.get_total_reflection() / total_weight) * 100.0
				 << "%\n";
			file << "  Transmission Probability: " << (simulator_->metrics.get_total_transmission() / total_weight) * 100.0
				 << "%\n";
		}
		file << "\nAverages Per Photon:\n";
		file << "  Path Length: " << simulator_->metrics.get_path_length() / simulator_->config.num_photons() << "\n";
		file << "  Scatter Events: " << simulator_->metrics.get_scatter_events() / simulator_->config.num_photons() << "\n";
		file << "  Step Size: " << simulator_->metrics.get_average_step_size() << "\n";
	}
	else {
		file << "Single Photon Results:\n";
		file << "  Path Length: " << simulator_->metrics.get_path_length() << "\n";
		file << "  Scatter Events: " << simulator_->metrics.get_scatter_events() << "\n";
		file << "  Average Step Size: " << simulator_->metrics.get_average_step_size() << "\n";
		file << "  Diffusion Distance: " << simulator_->metrics.get_diffusion_distance() << "\n";
	}

	file << "\nDetailed Results:\n";
	file << "  Total Absorption: " << simulator_->metrics.get_total_absorption() << "\n";
	file << "  Total Reflection: " << simulator_->metrics.get_total_reflection() << "\n";
	file << "  Total Transmission: " << simulator_->metrics.get_total_transmission() << "\n";
	file << "  Total Diffusion: " << simulator_->metrics.get_total_diffusion() << "\n";
	file << "  Surface Reflection: " << simulator_->metrics.get_surface_reflection() << "\n";
	file << "  Surface Refraction: " << simulator_->metrics.get_surface_refraction() << "\n";

	file << "\nTissue Properties:\n";
	auto tissues = simulator_->get_all_tissues();
	for (const auto& tissue : tissues) {
		file << "  Tissue " << tissue.id << ":\n";
		file << "    Refractive Index (eta): " << tissue.eta << "\n";
		file << "    Absorption Coefficient (mua): " << tissue.mu_a << " cm^-1\n";
		file << "    Scattering Coefficient (mus): " << tissue.mu_s << " cm^-1\n";
		file << "    Anisotropy Factor (g): " << tissue.g << "\n";
	}

	file.close();
	std::cout << "Results saved successfully to " << filepath << std::endl;
}
