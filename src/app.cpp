#include "app.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <cxxopts.hpp>

#include "renderer/overlay.hpp"
#include "renderer/renderer.hpp"
#include "simulator/simulator.hpp"
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

std::string App::get_output_path(const std::string& filename) {
	std::filesystem::path out_dir = std::filesystem::path(executable_directory_) / "out";
	
	// Create the out directory if it doesn't exist
	std::filesystem::create_directories(out_dir);
	
	return (out_dir / filename).string();
}

App::App() : window_(nullptr), window_width_(1200), window_height_(800), should_close_(false) {
	app_instance = this;
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
	cxxopts::Options options("Photron", "Monte Carlo Photon Transport Renderer");
	
	options.add_options()
		("c,config", "Configuration file path (optional for GUI mode)", cxxopts::value<std::string>())
		("headless", "Run simulation without GUI")
		("v,verbose", "Enable verbose initialization messages")
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
		bool verbose_mode = result.count("verbose") > 0;
		bool headless_mode = result.count("headless") > 0;

		// If headless mode is requested without a config file, that's an error
		if (headless_mode && !has_config_file) {
			std::cerr << "Error: Headless mode requires a configuration file" << std::endl;
			std::cerr << std::endl << options.help() << std::endl;
			return false;
		}

		// Always initialize Config (with defaults or from file)
		if (has_config_file) {
			Config::initialize(config_file_);
		} else {
			Config::initialize(); // Initialize with defaults
		}

		// Set verbose mode if requested
		if (verbose_mode) {
			Config::get().set_verbose(true);
		}

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
			std::cerr << "Failed to initialize GUI" << std::endl;
			return false;
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
		// Disable UI while loading
		overlay_->set_ui_enabled(false);

		// Update config file path
		config_file_ = filepath;
		
		// Reinitialize Config with the new file
		Config::shutdown();
		Config::initialize(config_file_);

		// Create simulator if it doesn't exist
		if (!simulator_) {
			if (!initialize_simulator()) {
				std::cerr << "Failed to initialize simulator" << std::endl;
				overlay_->set_ui_enabled(true);
				return;
			}
		}

		// Load new configuration into simulator
		simulator_->initialize(filepath.c_str());

		// Run simulation immediately
		std::cout << "Running simulation with config: " << filepath << std::endl;
		simulator_->simulate();
		simulator_->report();

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

// Helper methods to aggregate metrics from all mediums
double App::aggregate_path_length() const {
	if (!simulator_) return 0.0;
	
	double total_path_length = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		const auto& metrics = medium.get_metrics();
		total_path_length += metrics.compute_path_length();
	}
	return total_path_length;
}

double App::aggregate_scatter_events() const {
	if (!simulator_) return 0.0;
	
	double total_scatter_events = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		const auto& metrics = medium.get_metrics();
		total_scatter_events += metrics.get_scatter_events();
	}
	return total_scatter_events;
}

double App::aggregate_average_step_size() const {
	if (!simulator_) return 0.0;
	
	double total_step_size = 0.0;
	double total_steps = 0.0;
	const auto& mediums = simulator_->mediums;
	
	for (const auto& medium : mediums) {
		const auto& metrics = medium.get_metrics();
		double step_size = metrics.compute_average_step_size();
		if (step_size > 0.0) {
			// Weight by number of steps (path vertices - 1)
			size_t num_vertices = metrics.get_path_vertices().size();
			if (num_vertices > 1) {
				double num_steps = static_cast<double>(num_vertices - 1);
				total_step_size += step_size * num_steps;
				total_steps += num_steps;
			}
		}
	}
	
	return (total_steps > 0.0) ? (total_step_size / total_steps) : 0.0;
}

double App::aggregate_diffusion_distance() const {
	if (!simulator_) return 0.0;
	
	// For diffusion distance, we'll take the maximum across all mediums
	// since it represents the overall extent of photon travel
	double max_diffusion_distance = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		const auto& metrics = medium.get_metrics();
		double diffusion_distance = metrics.compute_diffusion_distance();
		max_diffusion_distance = std::max(max_diffusion_distance, diffusion_distance);
	}
	return max_diffusion_distance;
}

// Helper methods to aggregate energy conservation values from all mediums
double App::aggregate_total_absorption() const {
	if (!simulator_) return 0.0;
	
	double total_absorption = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		total_absorption += medium.get_metrics().get_total_absorption();
	}
	return total_absorption;
}

double App::aggregate_total_reflection() const {
	if (!simulator_) return 0.0;
	
	double total_reflection = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		total_reflection += medium.get_metrics().get_diffuse_reflection();
	}
	return total_reflection;
}

double App::aggregate_total_transmission() const {
	if (!simulator_) return 0.0;
	
	double total_transmission = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		total_transmission += medium.get_metrics().get_diffuse_transmission();
	}
	return total_transmission;
}

double App::aggregate_total_diffusion() const {
	return aggregate_total_reflection() + aggregate_total_transmission();
}

double App::aggregate_surface_reflection() const {
	if (!simulator_) return 0.0;
	
	double surface_reflection = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		surface_reflection += medium.get_metrics().get_surface_reflection();
	}
	return surface_reflection;
}

double App::aggregate_surface_refraction() const {
	if (!simulator_) return 0.0;
	
	double surface_refraction = 0.0;
	const auto& mediums = simulator_->mediums;
	for (const auto& medium : mediums) {
		surface_refraction += medium.get_metrics().get_surface_refraction();
	}
	return surface_refraction;
}

void App::save_results_as_json(const std::string& filepath) {
	if (!simulator_) {
		std::cerr << "No simulator available for saving results." << std::endl;
		return;
	}

	// Create results directory if it doesn't exist
	std::filesystem::path results_dir = std::filesystem::path(get_executable_directory()) / "out" / "results";
	if (!std::filesystem::exists(results_dir)) {
		std::filesystem::create_directories(results_dir);
	}
	
	// Construct full path in results subfolder
	std::filesystem::path full_path = results_dir / filepath;

	std::ofstream file(full_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing: " << full_path << std::endl;
		return;
	}

	// Save simulation results as JSON
	file << "{\n";
	
	// Sanitize config file path for JSON (escape backslashes)
	std::string sanitized_config = config_file_;
	std::replace(sanitized_config.begin(), sanitized_config.end(), '\\', '/');
	
	file << "  \"config_file\": \"" << sanitized_config << "\",\n";
	file << "  \"medium_statistics\": [\n";

	// Export per-medium statistics array
	auto& mediums = simulator_->mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		auto& medium = mediums[i];
		const auto& metrics = medium.get_metrics();
		const auto& volume = medium.get_volume();
		const auto& dimensions = volume.dimensions();
		
		// Count surface voxels
		size_t surface_count = 0;
		for (uint64_t idx = 0; idx < volume.size(); ++idx) {
			auto voxel = volume.at(static_cast<uint32_t>(idx));
			if (voxel && voxel->is_surface_voxel) {
				surface_count++;
			}
		}
		
		file << "    {\n";
		file << "      \"medium_id\": " << (i + 1) << ",\n";
		file << "      \"volume_statistics\": {\n";
		file << "        \"grid_size\": [" << dimensions.x << ", " << dimensions.y << ", " << dimensions.z << "],\n";
		file << "        \"voxel_size\": " << volume.voxel_size() << ",\n";
		file << "        \"total_voxels\": " << volume.size() << ",\n";
		file << "        \"surface_voxels\": " << surface_count << "\n";
		file << "      },\n";
		file << "      \"transport_statistics\": {\n";
		file << "        \"total_photons\": " << simulator_->get_paths().size() << ",\n";
		file << "        \"photons_entered\": " << medium.get_metrics().get_photons_entered() << ",\n";
		file << "        \"scatter_events\": " << metrics.get_scatter_events() << ",\n";
		file << "        \"path_length\": " << metrics.compute_path_length() << ",\n";
		file << "        \"average_step_size\": " << metrics.compute_average_step_size() << ",\n";
		file << "        \"diffusion_distance\": " << metrics.compute_diffusion_distance() << "\n";
		file << "      },\n";
		file << "      \"energy_conservation\": {\n";
		file << "        \"total_absorption\": " << medium.get_metrics().get_total_absorption() << ",\n";
		file << "        \"diffuse_reflection\": " << medium.get_metrics().get_diffuse_reflection() << ",\n";
		file << "        \"specular_reflection\": " << medium.get_metrics().get_surface_reflection() << ",\n";
		file << "        \"diffuse_transmission\": " << medium.get_metrics().get_diffuse_transmission() << ",\n";
		file << "        \"specular_transmission\": " << medium.get_metrics().get_specular_transmission() << ",\n";
		file << "        \"surface_refraction\": " << medium.get_metrics().get_surface_refraction() << "\n";
		file << "      },\n";
		file << "      \"tissue_properties\": [\n";
		
		// Export tissues for this medium
		auto& tissues = medium.get_tissues();
		for (size_t j = 0; j < tissues.size(); ++j) {
			const auto& material = tissues[j];
			file << "        {\n";
			file << "          \"id\": " << material.id() << ",\n";
			file << "          \"eta\": " << material.eta() << ",\n";
			file << "          \"mua\": " << material.mu_a() << ",\n";
			file << "          \"mus\": " << material.mu_s() << ",\n";
			file << "          \"ani\": " << material.g() << "\n";
			file << "        }" << (j < tissues.size() - 1 ? "," : "") << "\n";
		}
		file << "      ]\n";
		file << "    }" << (i < mediums.size() - 1 ? "," : "") << "\n";
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

	// Create results directory if it doesn't exist
	std::filesystem::path results_dir = std::filesystem::path(get_executable_directory()) / "out" / "results";
	if (!std::filesystem::exists(results_dir)) {
		std::filesystem::create_directories(results_dir);
	}
	
	// Construct full path in results subfolder
	std::filesystem::path full_path = results_dir / filepath;

	std::ofstream file(full_path);
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing: " << full_path << std::endl;
		return;
	}

	// Save simulation results as formatted text
	file << "=== Photron Simulation Results ===\n\n";
	file << "Configuration File: " << config_file_ << "\n";
	file << "Timestamp: " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << "\n\n";

	// Per-medium statistics section (matching overlay format)
	auto& mediums = simulator_->mediums;
	for (size_t i = 0; i < mediums.size(); ++i) {
		auto& medium = mediums[i];
		const auto& metrics = medium.get_metrics();
		const auto& volume = medium.get_volume();
		const auto& dimensions = volume.dimensions();
		
		file << "Medium " << (i + 1) << "\n";
		file << "Volume Statistics\n";
		file << "  Volume Grid:         " << dimensions.x << "x" << dimensions.y << "x" << dimensions.z << "\n";
		file << "  Total Voxels:        " << volume.size() << "\n";
		
		// Count surface voxels
		size_t surface_count = 0;
		for (uint64_t idx = 0; idx < volume.size(); ++idx) {
			auto voxel = volume.at(static_cast<uint32_t>(idx));
			if (voxel && voxel->is_surface_voxel) {
				surface_count++;
			}
		}
		file << "  Surface Voxels:      " << surface_count << "\n";
		
		file << "\nTransport Statistics\n";
		file << "  Total photons:       " << simulator_->get_paths().size() << "\n";
		file << "  Photons entered:     " << medium.get_metrics().get_photons_entered() << "\n";
		file << "  Scatter events:      " << static_cast<int>(metrics.get_scatter_events()) << "\n";
		file << "  Total path length:   " << std::fixed << std::setprecision(6) << metrics.compute_path_length() << "\n";
		file << "  Average step size:   " << std::fixed << std::setprecision(6) << metrics.compute_average_step_size() << "\n";
		file << "  Diffusion distance:  " << std::fixed << std::setprecision(6) << metrics.compute_diffusion_distance() << "\n";
		
		// Use unified energy conservation calculation (same as console output)
		auto energy = simulator_->calculate_energy_conservation();
		
		file << "\nRadiance Properties\n";
		file << "  Total absorption:    " << std::fixed << std::setprecision(6) << energy.total_absorption << "\n";
		file << "  Total diffusion:     " << std::fixed << std::setprecision(6) << energy.total_diffusion << "\n";
		file << "    Reflection:        " << std::fixed << std::setprecision(6) << energy.total_reflection << "\n";
		file << "    Transmission:      " << std::fixed << std::setprecision(6) << energy.total_transmission << "\n";
		
		file << "\nEnergy Conservation\n";
		if (energy.surface_refraction > 0) {
			// Use TOTAL initial energy as baseline (specular reflection + refracted energy)
			double baseline_energy = energy.surface_reflection + energy.surface_refraction;
			
			// Calculate percentages relative to total initial energy
			double surface_reflection_percent = (energy.surface_reflection / baseline_energy) * 100.0;
			double absorption_percent = (energy.total_absorption / baseline_energy) * 100.0;
			double reflection_percent = (energy.total_reflection / baseline_energy) * 100.0;
			double transmission_percent = (energy.total_transmission / baseline_energy) * 100.0;
			
			// Total should equal baseline_energy for perfect conservation
			double total_accounted = energy.surface_reflection + energy.total_absorption + 
									energy.total_reflection + energy.total_transmission;
			double total_percent = (total_accounted / baseline_energy) * 100.0;
			
			file << "  Surface reflection:  " << std::fixed << std::setprecision(1) << surface_reflection_percent << "%\n";
			file << "  Absorption:          " << std::fixed << std::setprecision(1) << absorption_percent << "%\n";
			file << "  Reflection:  		" << std::fixed << std::setprecision(1) << reflection_percent << "%\n";
			file << "  Transmission:        " << std::fixed << std::setprecision(1) << transmission_percent << "%\n";
			file << "  Total:               " << std::fixed << std::setprecision(1) << total_percent << "%\n";
		}
		else {
			file << "  Surface reflection:  0.0%\n";
			file << "  Absorption:          0.0%\n";
			file << "  Reflection:  		0.0%\n";
			file << "  Transmission:        0.0%\n";
			file << "  Total:               0.0%\n";
		}
		
		// Add material properties for this medium
		file << "\nTissue Properties\n";
		auto& tissues = medium.get_tissues();
		for (const auto& material : tissues) {
			file << "  material " << material.id() << ":\n";
			file << "    Refractive Index (eta): " << material.eta() << "\n";
			file << "    Absorption Coefficient (mua): " << material.mu_a() << " cm^-1\n";
			file << "    Scattering Coefficient (mus): " << material.mu_s() << " cm^-1\n";
			file << "    Anisotropy Factor (g): " << material.g() << "\n";
		}
		
		if (i < mediums.size() - 1) {
			file << "\n======================================\n\n";
		}
	}

	file.close();
	std::cout << "Results saved successfully to " << filepath << std::endl;
}

bool App::initialize_simulator() {
	// Initialize simulator (Config should already be initialized by this point)
	simulator_ = std::make_unique<Simulator>();
	
	if (!simulator_->initialize(config_file_.c_str())) {
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
		std::cerr << "Error: run_simulation_with_progress called without config file!" << std::endl;
		return;
	}
	
	if (!simulator_) {
		std::cerr << "Error: run_simulation_with_progress called without simulator!" << std::endl;
		return;
	}

	std::cout << "=== Starting Monte Carlo Simulation ===" << std::endl;
	std::cout << "Configuration: " << config_file_ << std::endl;
	std::cout << "Number of photons: " << Config::get().num_photons() << std::endl;
	std::cout << std::endl;
	
	// Run simulation
	simulator_->simulate();
	simulator_->report();
	
	std::cout << std::endl;
	std::cout << "=== Simulation Complete ===" << std::endl;
	std::cout << "Initializing GUI..." << std::endl;
}
