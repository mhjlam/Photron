/**
 * @file overlay.cpp
 * @brief Implementation of ImGui-based user interface overlay
 *
 * Implements the Overlay class providing comprehensive GUI controls for
 * Monte Carlo simulation parameters, real-time visualization settings,
 * and results analysis. Features optimized caching for performance and
 * intuitive file dialog management.
 */

#include "overlay.hpp"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "common/config.hpp"
#include "simulator/medium.hpp"
#include "simulator/photon.hpp"
#include "simulator/simulator.hpp"
#include "simulator/volume.hpp"
#include "simulator/voxel.hpp"

// Optimized numerical formatting with multi-slot caching
struct CachedFormat
{
	double last_value = -999999.0; // Use impossible value for invalid state
	char cached_string[16] = {0};
	bool is_valid = false;
};

// Multiple cache slots for concurrent value formatting
static CachedFormat format_cache[16];
static int cache_slot = 0;

static const char* format_8char(double value) {
	// Find existing cache slot or allocate new one
	CachedFormat* slot = nullptr;

	// Check if value is already cached for reuse
	for (int i = 0; i < 16; i++) {
		if (format_cache[i].is_valid && std::abs(format_cache[i].last_value - value) < 1e-9) {
			slot = &format_cache[i];
			break;
		}
	}

	// Use round-robin cache allocation for new values
	if (!slot) {
		slot = &format_cache[cache_slot];
		cache_slot = (cache_slot + 1) % 16;
	}

	// Recalculate format only if value changed
	if (!slot->is_valid || std::abs(slot->last_value - value) >= 1e-9) {
		slot->last_value = value;
		slot->is_valid = true;

		// Try precision from 6 down to 0 until it fits in 8 characters
		bool formatted = false;
		for (int precision = 6; precision >= 0 && !formatted; precision--) {
			snprintf(slot->cached_string, sizeof(slot->cached_string), "%.*f", precision, value);
			if (strlen(slot->cached_string) <= 8) {
				formatted = true;
			}
		}

		// Fallback formatting for very large numbers
		if (!formatted) {
			snprintf(slot->cached_string, sizeof(slot->cached_string), "%.0f", value);

			// Truncate if still too long
			if (strlen(slot->cached_string) > 8) {
				slot->cached_string[8] = '\0';
			}
		}
	}

	return slot->cached_string;
}

// ImGui helper function for consistent tooltips
static void HelpMarker(const char* desc) {
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered()) {
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

Overlay::Overlay() :
	ui_enabled_(true), show_file_dialog_(false), deferred_open_config_(false), current_directory_("config/"),
	show_save_feedback_(false), save_feedback_timer_(0.0f) {
	// Initialize file path buffer and dialog state
	std::fill(std::begin(file_path_buffer_), std::end(file_path_buffer_), '\0');
	file_path_.clear();
}

Overlay::~Overlay() {
	// Clean shutdown of ImGui systems
	shutdown();
}

void Overlay::open_config_dialog() {
	show_file_dialog_ = true;
}

void Overlay::open_config_dialog_deferred() {
	deferred_open_config_ = true;
}

bool Overlay::initialize(GLFWwindow* window) {
	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

	// Disable ImGui's cursor management so we can control it ourselves
	io.ConfigFlags |= ImGuiConfigFlags_NoMouseCursorChange;

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	// Setup Platform/Renderer backends
	if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
		std::cerr << "Failed to initialize ImGui GLFW backend" << std::endl;
		return false;
	}

	if (!ImGui_ImplOpenGL3_Init("#version 330")) {
		std::cerr << "Failed to initialize ImGui OpenGL3 backend" << std::endl;
		return false;
	}

	if (Config::get().log()) {
		std::cout << "ImGui overlay initialized successfully" << std::endl;
	}
	return true;
}

void Overlay::shutdown() {
	// Ensure we have valid ImGui context before shutdown
	if (ImGui::GetCurrentContext() != nullptr) {
		ImGui_ImplOpenGL3_Shutdown();
		ImGui_ImplGlfw_Shutdown();
		ImGui::DestroyContext();
	}
}

void Overlay::handle_keyboard_shortcuts(bool has_simulator) {
	ImGuiIO& io = ImGui::GetIO();

	// Only process shortcuts if UI is enabled and no text input is active
	if (!ui_enabled_ || io.WantTextInput) {
		return;
	}

	bool ctrl_pressed = io.KeyCtrl;

	// Ctrl+O: Load Config (always available)
	if (ctrl_pressed && ImGui::IsKeyPressed(ImGuiKey_O)) {
		show_file_dialog_ = true;
	}

	// R: Run Simulation (only when simulator available)
	if (has_simulator && ImGui::IsKeyPressed(ImGuiKey_R)) {
		if (!ctrl_pressed && run_simulation_callback_) {
			run_simulation_callback_();
		}
		// Ctrl+R: Rerun Simulation
		else if (ctrl_pressed && rerun_simulation_callback_) {
			rerun_simulation_callback_();
		}
	}

	// Ctrl+S: Save Results (only when simulator available)
	if (has_simulator && ctrl_pressed && ImGui::IsKeyPressed(ImGuiKey_S)) {
		if (direct_save_results_callback_) {
			direct_save_results_callback_();
		}
	}
}

void Overlay::render() {
	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// Check for deferred dialog opening
	if (deferred_open_config_) {
		show_file_dialog_ = true;
		deferred_open_config_ = false;
	}

	// Handle keyboard shortcuts
	handle_keyboard_shortcuts(false); // No simulator available

	render_main_menu_bar(false);      // No simulator available
	render_control_panel();

	// Render file dialog if needed
	if (show_file_dialog_) {
		render_file_dialog();
	}

	// Render save feedback if needed
	render_save_feedback();

	// Rendering
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Overlay::render_with_simulator(Simulator& simulator) {
	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// Check for deferred dialog opening
	if (deferred_open_config_) {
		show_file_dialog_ = true;
		deferred_open_config_ = false;
	}

	// Render world text overlays FIRST (so they appear behind other UI elements)
	render_world_text_overlays();

	// Handle keyboard shortcuts
	handle_keyboard_shortcuts(true); // Simulator is available

	render_main_menu_bar(true);      // Simulator is available
	render_control_panel(&simulator);

	// Render file dialog if needed
	if (show_file_dialog_) {
		render_file_dialog();
	}

	// Render save feedback if needed
	render_save_feedback();

	// Rendering
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Overlay::render_main_menu_bar(bool has_simulator) {
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("File")) {
			if (ImGui::MenuItem("Load...", "Ctrl+O", false, ui_enabled_)) {
				show_file_dialog_ = true;
			}
			ImGui::Separator();
			if (ImGui::MenuItem("Exit", "ESC")) {
				glfwSetWindowShouldClose(glfwGetCurrentContext(), GLFW_TRUE);
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("View")) {
			if (ImGui::MenuItem("Reset", nullptr, false, ui_enabled_ && has_simulator)) {
				if (reset_view_callback_) {
					reset_view_callback_();
				}
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Simulation")) {
			if (ImGui::MenuItem("Launch Photon", "R", false, ui_enabled_ && has_simulator)) {
				if (run_simulation_callback_) {
					run_simulation_callback_();
				}
			}
			if (ImGui::MenuItem("Rerun Simulation", "Ctrl+R", false, ui_enabled_ && has_simulator)) {
				if (rerun_simulation_callback_) {
					rerun_simulation_callback_();
				}
			}
			ImGui::Separator();
			if (ImGui::MenuItem("Save Statistics", "Ctrl+S", false, ui_enabled_ && has_simulator)) {
				if (direct_save_results_callback_) {
					direct_save_results_callback_();
				}
			}
			ImGui::EndMenu();
		}

		// Add config filename in center of menu bar (only if simulator is loaded)
		float menu_bar_width = ImGui::GetWindowWidth();
		std::string config_name = "No Configuration Loaded";
		if (has_simulator) {
			const std::string& config_path = Config::get().config_filename();
			std::filesystem::path path(config_path);
			config_name = path.filename().string();
		}

		float text_width = ImGui::CalcTextSize(config_name.c_str()).x;
		ImGui::SetCursorPosX((menu_bar_width - text_width) * 0.5f);
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.6f, 0.6f, 0.6f, 1.0f)); // Gray color
		ImGui::Text("%s", config_name.c_str());
		ImGui::PopStyleColor();

		// Add controls help in the upper-right corner
		ImGui::SetCursorPosX(menu_bar_width - 30);
		ImGui::Text("?");
		if (ImGui::IsItemHovered()) {
			ImGui::BeginTooltip();
			ImGui::Text("Controls");
			ImGui::Separator();
			ImGui::Text("Orbit Mode:");
			ImGui::BulletText("Left Mouse: Hold to orbit");
			ImGui::BulletText("Right Mouse: Hold to zoom");
			ImGui::BulletText("Scroll Wheel: Move along the Y axis");
			ImGui::Text("Free Mode:");
			ImGui::BulletText("WASD: Move forward/back/left/right");
			ImGui::BulletText("Space/Shift: Move up/down");
			ImGui::BulletText("Right Mouse: Hold to look around");
			ImGui::EndTooltip();
		}

		ImGui::EndMainMenuBar();
	}
}

void Overlay::render_control_panel(Simulator* simulator) {
	if (!simulator) {
		return;
	}

	// Set window position and size to be larger to accommodate more information
	ImGui::SetNextWindowPos(ImVec2(10, 30), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(240, 720), ImGuiCond_Always); // Set size every time

	// Create non-movable, non-resizable window like VolRec
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;

	if (ImGui::Begin("##MainWindow", nullptr, window_flags)) {
		// Rendering Options Section (only when simulator is available)
		ImGui::SeparatorText("Rendering");

		// Boolean toggles
		ImGui::Checkbox("Voxels", &settings_.draw_voxels);
		ImGui::SameLine();
		HelpMarker("Enable or disable voxel rendering");

		ImGui::Checkbox("Photons", &settings_.draw_paths);
		ImGui::SameLine();
		HelpMarker("Show photon paths and scatter points");

		// Energy Labels - auto-disabled after 10 photons (handled by renderer)
		bool many_photons = (simulator && simulator->photons.size() > 10);

		ImGui::Checkbox("Energy Labels", &settings_.draw_labels);
		ImGui::SameLine();
		if (many_photons) {
			HelpMarker(
				"Energy labels are auto-disabled after 10 photons to prevent clutter, but you can manually enable them");
		}
		else {
			HelpMarker("Show energy percentage labels at interaction points");
		}

		ImGui::Checkbox("Volume Geometry", &settings_.draw_volume);
		ImGui::SameLine();
		HelpMarker("Show actual geometry boundary with triangles/quads");

		ImGui::Spacing();

		// Camera mode
		{
			const char* camera_modes[] = {"Orbit", "Free"};
			int current_camera = static_cast<int>(settings_.camera_mode);
			ImGui::SetNextItemWidth(100.0f); // Make combo box shorter
			if (ImGui::Combo("Camera Mode", &current_camera, camera_modes, IM_ARRAYSIZE(camera_modes))) {
				settings_.camera_mode = static_cast<CameraMode>(current_camera);
				if (camera_mode_changed_callback_) {
					camera_mode_changed_callback_(settings_.camera_mode == CameraMode::Orbit);
				}
			}
			ImGui::SameLine();
			HelpMarker(
				"Orbit: Traditional orbit camera that rotates around the scene center.\nFree: FPS-style camera with "
				"WASD movement, Space/Shift for up/down, and mouse look.");
		}

		// Voxel mode (disabled when voxels are not drawn)
		{
			const char* voxel_modes[] = {"Layers", "Absorption", "Emittance"};
			int current_voxel = static_cast<int>(settings_.voxel_mode);
			ImGui::SetNextItemWidth(100.0f); // Make combo box shorter

			if (!settings_.draw_voxels) {
				ImGui::BeginDisabled();
			}

			if (ImGui::Combo("Voxel Mode", &current_voxel, voxel_modes, IM_ARRAYSIZE(voxel_modes))) {
				settings_.voxel_mode = static_cast<VoxelMode>(current_voxel);
			}

			if (!settings_.draw_voxels) {
				ImGui::EndDisabled();
			}

			ImGui::SameLine();
			HelpMarker("Visualization mode for voxels:\n"
					   "- Layers: Show layer colors with light transparency\n"
					   "- Absorption: Color by absorbed energy\n"
					   "- Emittance: Color by emitted energy");
		}

		ImGui::Spacing();

		// Simulation Info Section
		if (simulator) {
			ImGui::SeparatorText("Simulation Statistics");

			// General Simulation Statistics
			const auto& config = Config::get();
			ImGui::Text("Total Photons:       %llu", simulator->photons.size());
			ImGui::Text("Volume Grid:         %dx%dx%d", config.nx(), config.ny(), config.nz());
			ImGui::Text("Voxel Size:          %s", format_8char(config.vox_size()));

			// Calculate overall voxel statistics (using first medium for now)
			auto& mediums = simulator->mediums;
			if (!mediums.empty()) {
				const auto& volume = mediums[0].get_volume();
				size_t total_voxels = volume.size();
				size_t material_voxels = 0;

				for (uint64_t idx = 0; idx < volume.size(); ++idx) {
					auto voxel = volume.at(static_cast<uint32_t>(idx));
					if (voxel && voxel->material) {
						material_voxels++;
					}
				}

				ImGui::Text("Total Voxels:        %zu", total_voxels);
				ImGui::Text("Empty Voxels:        %zu", total_voxels - material_voxels);
			}

			// Per-medium statistics (collapsible)
			ImGui::Spacing();
			for (size_t i = 0; i < mediums.size(); ++i) {
				auto& medium = mediums[i];
				auto& metrics = medium.get_metrics();

				// Medium header with ID
				std::string medium_name = "Medium " + std::to_string(i + 1);
				if (ImGui::CollapsingHeader(medium_name.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
					// Medium Voxel Statistics
					const auto& volume = medium.get_volume();
					size_t material_voxels = 0;
					size_t surface_voxels = 0;

					for (uint64_t idx = 0; idx < volume.size(); ++idx) {
						auto voxel = volume.at(static_cast<uint32_t>(idx));
						if (voxel && voxel->material) {
							material_voxels++;
							if (voxel->is_surface_voxel) {
								surface_voxels++;
							}
						}
					}
					// Volume Statistics
					ImGui::Text("Volume Statistics");
					ImGui::Text("  Medium Voxels:       %zu", material_voxels);
					ImGui::Text("  Surface Voxels:      %zu", surface_voxels);
					ImGui::Spacing();

					// Transport Statistics (8-character adaptive formatting)
					ImGui::Text("Transport Statistics");
					ImGui::Text("  Photons entered:     %d", metrics.get_photons_entered());
					ImGui::Text("  Scatter events:      %.0f", metrics.get_scatter_events());
					ImGui::Text("  Total path length:   %s", format_8char(metrics.compute_path_length()));
					ImGui::Text("  Average step size:   %s", format_8char(metrics.compute_average_step_size()));
					ImGui::Text("  Diffusion distance:  %s", format_8char(metrics.compute_diffusion_distance()));

					ImGui::Spacing();

					// Use cached energy display data (single call, event-driven updates)
					const auto& energy_data = get_cached_energy_data(simulator);
					const auto& energy = energy_data.conservation;
					const auto& percentages = energy_data.percentages;

					ImGui::Text("Radiance Properties");
					ImGui::Text("  Total absorption:    %s", format_8char(energy.total_absorption));
					ImGui::Text("  Total diffusion:     %s", format_8char(energy.total_diffusion));
					ImGui::Text("    Reflection:        %s", format_8char(energy.total_reflection));
					ImGui::Text("    Transmission:      %s", format_8char(energy.total_transmission));

					ImGui::Spacing();

					// Energy Conservation
					ImGui::Text("Energy Conservation");

					if (energy_data.is_valid) {
						// Energy conservation total as percentage with color coding (red/green)
						ImVec4 total_color = (!percentages.is_conserved) ? ImVec4(1.0f, 0.3f, 0.3f, 1.0f)
																		 : ImVec4(0.3f, 1.0f, 0.3f, 1.0f);

						ImGui::Text("  Specularity:         %7.1f%%", percentages.surface_reflection_percent);
						ImGui::Text("  Absorption:          %7.1f%%", percentages.absorption_percent);
						ImGui::Text("  Reflection:          %7.1f%%", percentages.reflection_percent);
						ImGui::Text("  Transmission:        %7.1f%%", percentages.transmission_percent);
						ImGui::TextColored(total_color, "  Total:               %7.1f%%", percentages.total_percent);
					}
					else {
						ImGui::Text("  Specularity:         0.0%%");
						ImGui::Text("  Absorption:          0.0%%");
						ImGui::Text("  Reflection:          0.0%%");
						ImGui::Text("  Transmission:        0.0%%");
						ImGui::Text("  Total:               0.0%%");
					}

					ImGui::Spacing();
				}
			}
		}
	}
	ImGui::End();
}

void Overlay::render_file_dialog() {
	const char* title = "Load Configuration File";

	// Center the dialog on the screen
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 center = ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f);
	ImGui::SetNextWindowPos(center, ImGuiCond_Always, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_FirstUseEver);

	// Prevent the dialog from being minimized/collapsed
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoCollapse;
	if (ImGui::Begin(title, &show_file_dialog_, window_flags)) {
		// File path input
		ImGui::Text("File path:");
		ImGui::SetNextItemWidth(-1);
		if (ImGui::InputText("##filepath", file_path_buffer_, sizeof(file_path_buffer_))) {
			file_path_ = std::string(file_path_buffer_);
		}

		ImGui::Separator();

		// Quick access buttons for common directories and files
		ImGui::Text("Available Config Files:");
		ImGui::BeginChild("ConfigList", ImVec2(0, 200), true);

		// Dynamically list all .in files in the config directory
		try {
			for (const auto& entry : std::filesystem::directory_iterator("config")) {
				if (entry.is_regular_file()) {
					const auto& path = entry.path();
					if (path.extension() == ".toml") {
						std::string file_str = path.string();
						if (ImGui::Selectable(file_str.c_str())) {
							file_path_ = file_str;
							std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
							file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
						}
					}
				}
			}
		}
		catch (const std::exception& e) {
			ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error reading config directory!");
			std::cerr << "Error accessing config directory: " << e.what() << std::endl;
		}
		ImGui::EndChild();

		ImGui::Separator();

		// Dialog buttons
		if (ImGui::Button("OK") && !file_path_.empty()) {
			if (open_config_callback_) {
				open_config_callback_(file_path_);
			}
			show_file_dialog_ = false;
			std::fill(std::begin(file_path_buffer_), std::end(file_path_buffer_), '\0');
			file_path_.clear();
		}

		ImGui::SameLine();
		if (ImGui::Button("Cancel")) {
			show_file_dialog_ = false;
			std::fill(std::begin(file_path_buffer_), std::end(file_path_buffer_), '\0');
			file_path_.clear();
		}
	}
	ImGui::End();
}

void Overlay::add_world_text(const std::string& text, float screen_x, float screen_y, const glm::vec4& color) {
	world_text_queue_.push_back({text, screen_x, screen_y, color});
}

void Overlay::clear_world_text() {
	world_text_queue_.clear();
}

void Overlay::render_world_text_overlays() {
	// Only render if we have text to show
	if (world_text_queue_.empty()) {
		return;
	}

	// Create an invisible fullscreen window for text rendering behind UI
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove
									| ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse
									| ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoBackground
									| ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing
									| ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoDecoration;

	// Set window to cover entire viewport, positioned behind other windows
	ImGuiViewport* viewport = ImGui::GetMainViewport();
	ImGui::SetNextWindowPos(viewport->Pos);
	ImGui::SetNextWindowSize(viewport->Size);

	if (ImGui::Begin("##WorldTextOverlay", nullptr, window_flags)) {
		ImDrawList* draw_list = ImGui::GetWindowDrawList();

		for (const auto& world_text : world_text_queue_) {
			// Skip invalid coordinates
			if (world_text.screen_x < 0 || world_text.screen_y < 0) {
				continue;
			}

			// Convert glm::vec4 color to ImU32
			ImU32 color = IM_COL32(static_cast<int>(world_text.color.r * 255),
								   static_cast<int>(world_text.color.g * 255),
								   static_cast<int>(world_text.color.b * 255),
								   static_cast<int>(world_text.color.a * 255));

			// Position text using direct screen coordinates (already calculated relative to viewport)
			ImVec2 text_pos(world_text.screen_x, world_text.screen_y);

			// Black outline for readability
			for (int dx = -1; dx <= 1; dx++) {
				for (int dy = -1; dy <= 1; dy++) {
					if (dx != 0 || dy != 0) {
						draw_list->AddText(
							ImVec2(text_pos.x + dx, text_pos.y + dy), IM_COL32_BLACK, world_text.text.c_str());
					}
				}
			}

			// Main text
			draw_list->AddText(text_pos, color, world_text.text.c_str());
		}
	}
	ImGui::End();
}

void Overlay::show_save_feedback(const std::string& message) {
	save_feedback_message_ = message;
	show_save_feedback_ = true;
	save_feedback_timer_ = 3.0f; // Show for 3 seconds
}

void Overlay::render_save_feedback() {
	if (!show_save_feedback_)
		return;

	// Advance timer
	ImGuiIO& io = ImGui::GetIO();
	save_feedback_timer_ -= io.DeltaTime;

	if (save_feedback_timer_ <= 0.0f) {
		show_save_feedback_ = false;
		return;
	}

	// Render feedback popup in center of screen
	ImGuiViewport* viewport = ImGui::GetMainViewport();
	ImVec2 center = viewport->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Always, ImVec2(0.5f, 0.5f));

	ImGui::SetNextWindowBgAlpha(0.9f);
	auto flags = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove;
	if (ImGui::Begin("Save Status", nullptr, flags)) {
		// White text
		ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), "%s", save_feedback_message_.c_str());

		// Calculate text width to make progress bar match
		float text_width = ImGui::CalcTextSize(save_feedback_message_.c_str()).x;

		// Progress bar showing completion (0% to 100%)
		float progress = 1.0f - (save_feedback_timer_ / 3.0f);

		// Set nice pastel blue color for progress bar
		ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.6f, 0.8f, 1.0f, 0.8f)); // Pastel blue
		ImGui::ProgressBar(progress, ImVec2(text_width, 0.0f));
		ImGui::PopStyleColor();
	}
	ImGui::End();
}

/***********************************************************
 * ENERGY DATA CACHING
 * Event-driven caching to eliminate per-frame energy calculations
 ***********************************************************/
const Metrics::EnergyDisplayData& Overlay::get_cached_energy_data(Simulator* simulator) const {
	if (!simulator) {
		// Return invalid data if no simulator
		static const Metrics::EnergyDisplayData invalid_data {};
		return invalid_data;
	}

	const uint64_t current_version = simulator->get_simulation_version();

	// Check if we need to refresh the cache (version changed = new photons or rerun)
	if (!cached_energy_data_ || cached_simulation_version_ != current_version) {
		cached_energy_data_ = simulator->get_metrics().get_energy_display_data(*simulator);
		cached_simulation_version_ = current_version;
	}

	return *cached_energy_data_;
}
