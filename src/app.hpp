#pragma once

#include <GLFW/glfw3.h>
#include <memory>
#include <string>

class Overlay;
class Renderer;
class Simulator;

class App {
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

private:
    void setup_callbacks();
    void update();
    void render();

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
};
