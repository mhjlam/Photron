#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "app.hpp"
#include "renderer/renderer.hpp"
#include "overlay.hpp"
#include "simulator/simulator.hpp"
#include "utilities/utilities.hpp"

#include <iostream>
#include <stdexcept>

// Static pointer for callbacks
static App* app_instance = nullptr;

App::App() 
    : window_(nullptr)
    , window_width_(1200)
    , window_height_(800)
    , should_close_(false) {
    app_instance = this;
}

App::~App() {
    shutdown();
    app_instance = nullptr;
}

bool App::initialize(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc < 2) {
        std::cerr << "Error: program requires an input file." << std::endl;
        return false;
    }

    config_file_ = argv[1];

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
    if (!simulator_->initialize(config_file_.c_str())) {
        std::cerr << "Failed to initialize simulator" << std::endl;
        return false;
    }

    // Run simulation
    std::cout << "Running photon transport simulation..." << std::endl;
    simulator_->simulate();
    simulator_->report();
    std::cout << "Simulation completed!" << std::endl;

    // Stop if rendering is disabled
    if (argc > 2 && equals(argv[2], "norender")) {
        should_close_ = true;
        return true;
    }

    // Initialize renderer
    renderer_ = std::make_unique<Renderer>();
    if (!renderer_->initialize()) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return false;
    }

    // Initialize overlay
    overlay_ = std::make_unique<Overlay>();
    if (!overlay_->initialize(window_)) {
        std::cerr << "Failed to initialize overlay" << std::endl;
        return false;
    }

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

void App::update() {
    // Update camera and other systems
    if (renderer_) {
        renderer_->update();
    }
}

void App::render() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Update renderer settings from overlay before rendering
    if (renderer_ && overlay_) {
        renderer_->set_settings(overlay_->get_settings());
    }
    
    // Render 3D scene FIRST
    if (renderer_ && simulator_) {
        renderer_->render(simulator_.get());
    }
    
    // Render ImGui overlay AFTER 3D scene (so it appears on top)
    if (overlay_) {
        // Reset OpenGL state for ImGui
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        overlay_->render_with_simulator(simulator_.get());
        
        // Restore depth testing for next frame
        glEnable(GL_DEPTH_TEST);
    }
}

// GLFW callbacks
void App::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
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
            case GLFW_KEY_ESCAPE:
                glfwSetWindowShouldClose(window, GLFW_TRUE);
                break;
        }
    }
    
    if (app_instance && app_instance->renderer_) {
        app_instance->renderer_->handle_key_input(key, scancode, action, mods);
    }
}

void App::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos) {
    if (app_instance && app_instance->renderer_) {
        app_instance->renderer_->handle_mouse_move(static_cast<float>(xpos), static_cast<float>(ypos));
    }
}

void App::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (app_instance && app_instance->renderer_) {
        app_instance->renderer_->handle_mouse_button(button, action, mods);
    }
}

void App::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    if (app_instance && app_instance->renderer_) {
        app_instance->renderer_->handle_mouse_scroll(static_cast<float>(xoffset), static_cast<float>(yoffset));
    }
}
