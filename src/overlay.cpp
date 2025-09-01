#include "overlay.hpp"
#include "simulator/simulator.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <iostream>

Overlay::Overlay()
    : show_options_window_(true)
    , show_info_window_(false)
    , show_demo_window_(false) {
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

void Overlay::render_with_simulator(Simulator* simulator) {
    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    render_main_menu_bar();
    render_control_panel(simulator);

    // Rendering
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Overlay::render_main_menu_bar() {
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Open Configuration...")) {
                // TODO: Implement file dialog for configuration
            }
            if (ImGui::MenuItem("Save Screenshot...")) {
                // TODO: Implement screenshot saving
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) {
                // TODO: Implement exit
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("View")) {
            ImGui::MenuItem("Draw Bounds", nullptr, &settings_.draw_bounds);
            ImGui::MenuItem("Draw Frame", nullptr, &settings_.draw_frame);
            ImGui::MenuItem("Draw Paths", nullptr, &settings_.draw_paths);
            
            ImGui::Separator();
            
            // Voxel mode submenu
            if (ImGui::BeginMenu("Voxel Mode")) {
                bool none_selected = (settings_.voxel_mode == VoxelMode::None);
                bool absorption_selected = (settings_.voxel_mode == VoxelMode::Absorption);
                bool emittance_selected = (settings_.voxel_mode == VoxelMode::Emittance);
                
                if (ImGui::MenuItem("None", nullptr, none_selected)) {
                    settings_.voxel_mode = VoxelMode::None;
                }
                if (ImGui::MenuItem("Absorption", nullptr, absorption_selected)) {
                    settings_.voxel_mode = VoxelMode::Absorption;
                }
                if (ImGui::MenuItem("Emittance", nullptr, emittance_selected)) {
                    settings_.voxel_mode = VoxelMode::Emittance;
                }
                ImGui::EndMenu();
            }
            
            // Geometry mode submenu
            if (ImGui::BeginMenu("Geometry Mode")) {
                bool none_selected = (settings_.geometry_mode == GeometryMode::None);
                bool white_selected = (settings_.geometry_mode == GeometryMode::White);
                bool color_selected = (settings_.geometry_mode == GeometryMode::Color);
                
                if (ImGui::MenuItem("None", nullptr, none_selected)) {
                    settings_.geometry_mode = GeometryMode::None;
                }
                if (ImGui::MenuItem("White", nullptr, white_selected)) {
                    settings_.geometry_mode = GeometryMode::White;
                }
                if (ImGui::MenuItem("Color", nullptr, color_selected)) {
                    settings_.geometry_mode = GeometryMode::Color;
                }
                ImGui::EndMenu();
            }
            
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("Simulation")) {
            if (ImGui::MenuItem("Run Simulation...")) {
                // TODO: Implement simulation dialog
            }
            if (ImGui::MenuItem("Load Results...")) {
                // TODO: Implement results loading
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Reset View")) {
                // TODO: Reset camera to default position
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("Help")) {
            if (ImGui::MenuItem("About")) {
                // TODO: Show about dialog
            }
            if (ImGui::MenuItem("Controls")) {
                // TODO: Show controls help
            }
            ImGui::EndMenu();
        }
        
        ImGui::EndMainMenuBar();
    }
}

void Overlay::render_options_window() {
    ImGui::Begin("Render Options", &show_options_window_);
    
    // Camera mode
    {
        const char* camera_modes[] = { "Arc", "Free" };
        int current_camera = static_cast<int>(settings_.camera_mode);
        if (ImGui::Combo("Camera Mode", &current_camera, camera_modes, IM_ARRAYSIZE(camera_modes))) {
            settings_.camera_mode = static_cast<CameraMode>(current_camera);
        }
    }
    
    // Voxel mode
    {
        const char* voxel_modes[] = { "None", "Absorption", "Emittance" };
        int current_voxel = static_cast<int>(settings_.voxel_mode);
        if (ImGui::Combo("Voxel Mode", &current_voxel, voxel_modes, IM_ARRAYSIZE(voxel_modes))) {
            settings_.voxel_mode = static_cast<VoxelMode>(current_voxel);
        }
    }
    
    // Text mode
    {
        const char* text_modes[] = { "None", "Heads Up Display", "All" };
        int current_text = static_cast<int>(settings_.text_mode);
        if (ImGui::Combo("Text Mode", &current_text, text_modes, IM_ARRAYSIZE(text_modes))) {
            settings_.text_mode = static_cast<TextMode>(current_text);
        }
    }
    
    // Grid mode
    {
        const char* grid_modes[] = { "None", "Partial", "All" };
        int current_grid = static_cast<int>(settings_.grid_mode);
        if (ImGui::Combo("Grid Mode", &current_grid, grid_modes, IM_ARRAYSIZE(grid_modes))) {
            settings_.grid_mode = static_cast<GridMode>(current_grid);
        }
    }
    
    // Geometry mode
    {
        const char* geometry_modes[] = { "None", "White", "Color" };
        int current_geometry = static_cast<int>(settings_.geometry_mode);
        if (ImGui::Combo("Geometry Mode", &current_geometry, geometry_modes, IM_ARRAYSIZE(geometry_modes))) {
            settings_.geometry_mode = static_cast<GeometryMode>(current_geometry);
        }
    }
    
    ImGui::Separator();
    
    // Boolean toggles
    ImGui::Checkbox("Draw Paths", &settings_.draw_paths);
    ImGui::Checkbox("Draw Frame", &settings_.draw_frame);
    ImGui::Checkbox("Draw Bounds", &settings_.draw_bounds);
    
    ImGui::End();
}

void Overlay::render_info_window(Simulator* simulator) {
    ImGui::Begin("Simulation Info", &show_info_window_);
    
    if (simulator) {
        // Display simulation information
        ImGui::Text("Simulation Status: Complete");
        ImGui::Separator();
        
        // Display configuration info
        ImGui::Text("Grid Size: %dx%dx%d", simulator->config.nx, simulator->config.ny, simulator->config.nz);
        ImGui::Text("Voxel Size: %.4f", simulator->config.vox_size);
        ImGui::Text("Total Photons: %lld", simulator->config.num_photons);
        
        ImGui::Separator();
        
        // Display bounds
        ImGui::Text("Bounds:");
        ImGui::Text("  X: [%.3f, %.3f]", simulator->bounds.x_min, simulator->bounds.x_max);
        ImGui::Text("  Y: [%.3f, %.3f]", simulator->bounds.y_min, simulator->bounds.y_max);
        ImGui::Text("  Z: [%.3f, %.3f]", simulator->bounds.z_min, simulator->bounds.z_max);
        
        ImGui::Separator();
        
        // Display voxel count
        ImGui::Text("Total Voxels: %zu", simulator->voxels.size());
    } else {
        ImGui::Text("No simulation data available");
    }
    
    ImGui::End();
}

void Overlay::render_control_panel(Simulator* simulator) {
    // Set window position and size to be static at top-left
    ImGui::SetNextWindowPos(ImVec2(10, 30), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);
    
    // Create non-movable, non-resizable window like VolRec
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;
    
    if (ImGui::Begin("Control Panel", nullptr, window_flags)) {
        // Rendering Options Section
        ImGui::SeparatorText("Rendering Options");
        
        // Boolean toggles
        ImGui::Checkbox("Draw Bounds", &settings_.draw_bounds);
        ImGui::Checkbox("Draw Frame", &settings_.draw_frame);
        ImGui::Checkbox("Draw Paths", &settings_.draw_paths);
        
        ImGui::Spacing();
        
        // Camera mode
        {
            const char* camera_modes[] = { "Arc", "Free" };
            int current_camera = static_cast<int>(settings_.camera_mode);
            if (ImGui::Combo("Camera Mode", &current_camera, camera_modes, IM_ARRAYSIZE(camera_modes))) {
                settings_.camera_mode = static_cast<CameraMode>(current_camera);
            }
        }
        
        // Voxel mode
        {
            const char* voxel_modes[] = { "None", "Absorption", "Emittance" };
            int current_voxel = static_cast<int>(settings_.voxel_mode);
            if (ImGui::Combo("Voxel Mode", &current_voxel, voxel_modes, IM_ARRAYSIZE(voxel_modes))) {
                settings_.voxel_mode = static_cast<VoxelMode>(current_voxel);
            }
        }
        
        // Grid mode
        {
            const char* grid_modes[] = { "None", "Partial", "All" };
            int current_grid = static_cast<int>(settings_.grid_mode);
            if (ImGui::Combo("Grid Mode", &current_grid, grid_modes, IM_ARRAYSIZE(grid_modes))) {
                settings_.grid_mode = static_cast<GridMode>(current_grid);
            }
        }
        
        // Geometry mode
        {
            const char* geometry_modes[] = { "None", "White", "Color" };
            int current_geometry = static_cast<int>(settings_.geometry_mode);
            if (ImGui::Combo("Geometry Mode", &current_geometry, geometry_modes, IM_ARRAYSIZE(geometry_modes))) {
                settings_.geometry_mode = static_cast<GeometryMode>(current_geometry);
            }
        }
        
        // Simulation Info Section
        if (simulator) {
            ImGui::SeparatorText("Simulation Info");
            
            ImGui::Text("Status: Complete");
            ImGui::Text("Grid: %dx%dx%d", simulator->config.nx, simulator->config.ny, simulator->config.nz);
            ImGui::Text("Voxel Size: %.4f", simulator->config.vox_size);
            ImGui::Text("Photons: %lld", simulator->config.num_photons);
            ImGui::Text("Voxels: %zu", simulator->voxels.size());
            
            ImGui::Spacing();
            ImGui::Text("Bounds:");
            ImGui::Text("  X: [%.3f, %.3f]", simulator->bounds.x_min, simulator->bounds.x_max);
            ImGui::Text("  Y: [%.3f, %.3f]", simulator->bounds.y_min, simulator->bounds.y_max);
            ImGui::Text("  Z: [%.3f, %.3f]", simulator->bounds.z_min, simulator->bounds.z_max);
        }
    }
    ImGui::End();
}
