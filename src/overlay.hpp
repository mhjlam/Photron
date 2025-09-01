#pragma once

#include <GLFW/glfw3.h>
#include "structs/settings.hpp"

class Simulator;

class Overlay {
public:
    Overlay();
    ~Overlay();

    bool initialize(GLFWwindow* window);
    void shutdown();
    void render();
    void render_with_simulator(Simulator* simulator);

    // Getters for settings
    const Settings& get_settings() const { return settings_; }
    Settings& get_settings() { return settings_; }

private:
    void render_main_menu_bar();
    void render_options_window();
    void render_info_window(Simulator* simulator);
    void render_control_panel(Simulator* simulator = nullptr);

    Settings settings_;
    bool show_options_window_;
    bool show_info_window_;
    bool show_demo_window_;
};
