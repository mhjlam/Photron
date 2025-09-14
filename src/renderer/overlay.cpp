#include "overlay.hpp"

#include <cstring>
#include <filesystem>
#include <iostream>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "simulator/simulator.hpp"

// Helper function for ImGui tooltips
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
	show_options_window_(true), show_info_window_(false), show_demo_window_(false), ui_enabled_(true),
	show_file_dialog_(false), file_dialog_mode_(FileDialogMode::LoadConfig), current_directory_("config/") {
	std::fill(std::begin(file_path_buffer_), std::end(file_path_buffer_), '\0');
	file_path_.clear();
}

Overlay::~Overlay() {
	shutdown();
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

	std::cout << "ImGui overlay initialized successfully" << std::endl;
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

void Overlay::handle_keyboard_shortcuts() {
	ImGuiIO& io = ImGui::GetIO();

	// Only process shortcuts if UI is enabled and no text input is active
	if (!ui_enabled_ || io.WantTextInput) {
		return;
	}

	bool ctrl_pressed = io.KeyCtrl;

	// Ctrl+O: Load Config
	if (ctrl_pressed && ImGui::IsKeyPressed(ImGuiKey_O)) {
		show_file_dialog_ = true;
		file_dialog_mode_ = FileDialogMode::LoadConfig;
	}

	// R: Run Simulation
	if (ImGui::IsKeyPressed(ImGuiKey_R)) {
		if (!ctrl_pressed && run_simulation_callback_) {
			run_simulation_callback_();
		}
		// Shift+R: Rerun Simulation
		else if (ctrl_pressed && rerun_simulation_callback_) {
			rerun_simulation_callback_();
		}
	}

	// Ctrl+S: Save Results
	if (ctrl_pressed && ImGui::IsKeyPressed(ImGuiKey_S)) {
		show_file_dialog_ = true;
		file_dialog_mode_ = FileDialogMode::SaveResults;
	}
}

void Overlay::render() {
	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	render_main_menu_bar();
	render_control_panel();

	// Rendering
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Overlay::render_with_simulator(Simulator& simulator) {
	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// Render world text overlays FIRST (so they appear behind other UI elements)
	render_world_text_overlays();

	// Handle keyboard shortcuts
	handle_keyboard_shortcuts();

	render_main_menu_bar();
	render_control_panel(&simulator);

	// Render file dialog if needed
	if (show_file_dialog_) {
		render_file_dialog();
	}

	// Rendering
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Overlay::render_main_menu_bar() {
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("File")) {
			if (ImGui::MenuItem("Load Config...", "Ctrl+O", false, ui_enabled_)) {
				show_file_dialog_ = true;
				file_dialog_mode_ = FileDialogMode::LoadConfig;
			}
			ImGui::Separator();
			if (ImGui::MenuItem("Exit")) {
				glfwSetWindowShouldClose(glfwGetCurrentContext(), GLFW_TRUE);
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Simulation")) {
			if (ImGui::MenuItem("Add Run", "R", false, ui_enabled_)) {
				if (run_simulation_callback_) {
					run_simulation_callback_();
				}
			}
			if (ImGui::MenuItem("Re-Run", "Ctrl+R", false, ui_enabled_)) {
				if (rerun_simulation_callback_) {
					rerun_simulation_callback_();
				}
			}
			ImGui::Separator();
			if (ImGui::MenuItem("Save Results...", "Ctrl+S", false, ui_enabled_)) {
				show_file_dialog_ = true;
				file_dialog_mode_ = FileDialogMode::SaveResults;
			}
			ImGui::EndMenu();
		}

		// Add controls help in the upper-right corner like VolRec
		float menu_bar_width = ImGui::GetWindowWidth();
		ImGui::SetCursorPosX(menu_bar_width - 30);
		ImGui::Text("?");
		if (ImGui::IsItemHovered()) {
			ImGui::BeginTooltip();
			ImGui::Text("Controls");
			ImGui::Separator();
			ImGui::Text("ESC: Exit application");
			ImGui::Text("Orbit Mode:");
			ImGui::BulletText("Left Mouse: Hold to orbit");
			ImGui::BulletText("Right Mouse: Hold to zoom");
			ImGui::BulletText("Scroll Wheel: Move along the medium layers");
			ImGui::Text("Free Mode:");
			ImGui::BulletText("WASD: Move forward/back/left/right");
			ImGui::BulletText("QE: Move up/down");
			ImGui::BulletText("Right Mouse: Hold to look around");
			ImGui::EndTooltip();
		}

		ImGui::EndMainMenuBar();
	}
}

void Overlay::render_control_panel(Simulator* simulator) {
	// Set window position and size to be larger to accommodate more information
	ImGui::SetNextWindowPos(ImVec2(10, 30), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(240, 480), ImGuiCond_Always); // Force size update every time

	// Create non-movable, non-resizable window like VolRec
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;

	if (ImGui::Begin("##MainWindow", nullptr, window_flags)) {
		// Rendering Options Section
		ImGui::SeparatorText("Rendering");

		// Boolean toggles
		ImGui::Checkbox("Voxels", &settings_.draw_voxels);
		ImGui::SameLine();
		HelpMarker("Enable or disable voxel rendering");

		ImGui::Checkbox("Photons", &settings_.draw_paths);
		ImGui::SameLine();
		HelpMarker("Show photon paths and scatter points");

		ImGui::Checkbox("Geometry Bounds", &settings_.draw_volume);
		ImGui::SameLine();
		HelpMarker("Show actual geometry boundary with triangles/quads");

		// Energy Labels - auto-disabled after 10 photons (handled by renderer)
		bool many_photons = (simulator && simulator->paths.size() > 10);
		
		ImGui::Checkbox("Energy Labels", &settings_.draw_labels);
		ImGui::SameLine();
		if (many_photons) {
			HelpMarker("Energy labels are auto-disabled after 10 photons to prevent clutter, but you can manually enable them");
		}
		else {
			HelpMarker("Show energy percentage labels at interaction points");
		}

		ImGui::Spacing();

		// Camera mode
		{
			const char* camera_modes[] = {"Orbit", "Free"};
			int current_camera = static_cast<int>(settings_.camera_mode);
			ImGui::SetNextItemWidth(120.0f); // Make combo box shorter
			if (ImGui::Combo("Camera Mode", &current_camera, camera_modes, IM_ARRAYSIZE(camera_modes))) {
				settings_.camera_mode = static_cast<CameraMode>(current_camera);
				if (camera_mode_changed_callback_) {
					camera_mode_changed_callback_(settings_.camera_mode == CameraMode::Orbit);
				}
			}
			ImGui::SameLine();
			HelpMarker(
				"Orbit: Traditional orbit camera that rotates around the scene center.\nFree: FPS-style camera with "
				"WASD movement and mouse look.");
		}

		// Voxel mode (disabled when voxels are not drawn)
		{
			const char* voxel_modes[] = {"Layers", "Absorption", "Emittance"};
			int current_voxel = static_cast<int>(settings_.voxel_mode);
			ImGui::SetNextItemWidth(120.0f); // Make combo box shorter

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
			HelpMarker(
				"Visualization mode for voxels:\n• Layers: Show layer colors with light transparency\n• Absorption: Color by absorbed energy\n• "
				"Emittance: Color by emitted energy");
		}

		ImGui::Spacing();

		// Reset View button
		if (ImGui::Button("Reset View") && ui_enabled_) {
			if (reset_view_callback_) {
				reset_view_callback_();
			}
		}

		ImGui::Spacing();

		// Simulation Info Section
		if (simulator) {
			ImGui::SeparatorText("Simulation Statistics");

			ImGui::Text("Grid: %dx%dx%d", simulator->config.nx(), simulator->config.ny(), simulator->config.nz());
			ImGui::Text("Voxels: %zu (size: %.4f)", simulator->get_total_voxel_count(), simulator->config.vox_size());

			ImGui::Spacing();

			// Basic simulation info
			ImGui::Text("Photons Traced: %zu", simulator->paths.size());
			
			// Per-medium statistics (collapsible)
			ImGui::Spacing();
			auto& mediums = simulator->mediums;
			for (size_t i = 0; i < mediums.size(); ++i) {
				auto& medium = mediums[i];
				const auto& record = medium.get_record();
				auto& metrics = medium.get_metrics();
				
				// Medium header with ID
				std::string medium_name = "Medium " + std::to_string(i + 1);
				if (ImGui::CollapsingHeader(medium_name.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
					// Grid info
					ImGui::Text("Total Voxels: %llu", medium.get_volume().size());
					
					// Count surface voxels manually
					size_t surface_count = 0;
					const auto& volume = medium.get_volume();
					for (uint64_t idx = 0; idx < volume.size(); ++idx) {
						auto voxel = volume.at(static_cast<uint32_t>(idx));
						if (voxel && voxel->is_surface_voxel) {
							surface_count++;
						}
					}
					ImGui::Text("Surface Voxels: %zu", surface_count);
					
					ImGui::Spacing();
					
					// Transport Statistics (matching console format)
					ImGui::Text("Transport Statistics");
					ImGui::Text("  Path length:        %7.6f", metrics.compute_path_length());
					ImGui::Text("  Scatter events:     %.0f", metrics.get_scatter_events());
					ImGui::Text("  Average step size:  %.6f", metrics.compute_average_step_size());
					ImGui::Text("  Diffusion distance: %.6f", metrics.compute_diffusion_distance());
					
					ImGui::Spacing();
					
					// Surface Interactions (matching console format)
					ImGui::Text("Surface Interaction");
					ImGui::Text("  Surface reflection: %.6f", record.specular_reflection);
					ImGui::Text("  Surface refraction: %.6f", record.surface_refraction);

					ImGui::Spacing();
					
					// CORRECTED VOXEL-BASED ENERGY CALCULATIONS (matching console format)
					double total_voxel_absorption = 0.0;
					double total_voxel_emittance = 0.0;
					
					// Calculate voxel-based totals using existing volume variable
					for (uint64_t idx = 0; idx < volume.size(); ++idx) {
						auto voxel = volume.at(static_cast<uint32_t>(idx));
						if (voxel) {
							total_voxel_absorption += voxel->absorption;
							total_voxel_emittance += voxel->emittance;
						}
					}
					
					// Radiance Properties (using direction-classified medium records)
					double total_reflection = record.diffuse_reflection; // Direction-classified reflection
					double total_transmission = record.diffuse_transmission; // Direction-classified transmission  
					double total_diffusion = record.diffuse_reflection + record.diffuse_transmission;

					ImGui::Text("Radiance Properties");
					ImGui::Text("  Total absorption:   %.6f", total_voxel_absorption);
					ImGui::Text("  Total diffusion:    %.6f", total_diffusion);
					ImGui::Text("    Reflection:       %.6f", total_reflection);
					ImGui::Text("    Transmission:     %.6f", total_transmission);
					
					ImGui::Spacing();
					
					// Energy Conservation (using corrected voxel-based totals)
					double corrected_total_energy = total_voxel_absorption + total_reflection + total_transmission;
					ImGui::Text("Energy Conservation");
					if (record.surface_refraction > 0) {
						double conservation_ratio = corrected_total_energy / record.surface_refraction;
						
						if (corrected_total_energy > 0) {
							ImGui::Text("  Absorption:         %.1f%%", (total_voxel_absorption / record.surface_refraction) * 100.0);
							ImGui::Text("  Reflection:         %.1f%%", (total_reflection / record.surface_refraction) * 100.0);
							ImGui::Text("  Transmission:       %.1f%%", (total_transmission / record.surface_refraction) * 100.0);
							
							// Energy conservation status
							if (std::abs(conservation_ratio - 1.0) > 0.05) {
								ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "  WARNING: Energy not conserved!");
								ImGui::Text("  Ratio = %.3f", conservation_ratio);
							} else {
								ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "  Energy conserved!");
								ImGui::Text("  Ratio = %.3f", conservation_ratio);
							}
						}
					}
				}
			}
		}
	}
	ImGui::End();
}

void Overlay::render_file_dialog() {
	const char* title =
		(file_dialog_mode_ == FileDialogMode::LoadConfig) ? "Load Configuration File" : "Save Results File";

	ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowPos(ImVec2(200, 100), ImGuiCond_FirstUseEver);

	if (ImGui::Begin(title, &show_file_dialog_, ImGuiWindowFlags_Modal)) {
		// File path input
		ImGui::Text("File path:");
		ImGui::SetNextItemWidth(-1);
		if (ImGui::InputText("##filepath", file_path_buffer_, sizeof(file_path_buffer_))) {
			file_path_ = std::string(file_path_buffer_);
		}

		ImGui::Separator();

		// Quick access buttons for common directories and files
		if (file_dialog_mode_ == FileDialogMode::LoadConfig) {
			ImGui::Text("Quick Select:");
			if (ImGui::Button("config1.in")) {
				file_path_ = "config/config1.in";
				std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
				file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
			}
			ImGui::SameLine();
			if (ImGui::Button("config2.in")) {
				file_path_ = "config/config2.in";
				std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
				file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
			}
			ImGui::SameLine();
			if (ImGui::Button("multi-photon-100.in")) {
				file_path_ = "config/multi-photon-100.in";
				std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
				file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
			}

			// List available config files
			ImGui::Separator();
			ImGui::Text("Available Config Files:");
			ImGui::BeginChild("ConfigList", ImVec2(0, 200), true);

			// Dynamically list all .in files in the config directory
			try {
				for (const auto& entry : std::filesystem::directory_iterator("config")) {
					if (entry.is_regular_file()) {
						const auto& path = entry.path();
						if (path.extension() == ".in") {
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
		}
		else { // Save Results
			ImGui::Text("Save as:");
			if (ImGui::Button("results.json")) {
				file_path_ = "results.json";
				std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
				file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
			}
			ImGui::SameLine();
			if (ImGui::Button("simulation_data.txt")) {
				file_path_ = "simulation_data.txt";
				std::strncpy(file_path_buffer_, file_path_.c_str(), sizeof(file_path_buffer_) - 1);
				file_path_buffer_[sizeof(file_path_buffer_) - 1] = '\0';
			}
		}

		ImGui::Separator();

		// Dialog buttons
		if (ImGui::Button("OK") && !file_path_.empty()) {
			if (file_dialog_mode_ == FileDialogMode::LoadConfig) {
				if (open_config_callback_) {
					open_config_callback_(file_path_);
				}
			}
			else {
				if (save_results_callback_) {
					save_results_callback_(file_path_);
				}
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
			ImU32 color =
				IM_COL32(static_cast<int>(world_text.color.r * 255), static_cast<int>(world_text.color.g * 255),
						 static_cast<int>(world_text.color.b * 255), static_cast<int>(world_text.color.a * 255));

			// Position text using direct screen coordinates (already calculated relative to viewport)
			ImVec2 text_pos(world_text.screen_x, world_text.screen_y);

			// Black outline for readability
			for (int dx = -1; dx <= 1; dx++) {
				for (int dy = -1; dy <= 1; dy++) {
					if (dx != 0 || dy != 0) {
						draw_list->AddText(ImVec2(text_pos.x + dx, text_pos.y + dy), IM_COL32_BLACK,
										   world_text.text.c_str());
					}
				}
			}

			// Main text
			draw_list->AddText(text_pos, color, world_text.text.c_str());
		}
	}
	ImGui::End();
}
