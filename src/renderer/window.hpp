#pragma once

#include <Windows.h>
#include <GL/gl.h>
#include <iostream>
#include <string>

// Mouse state for camera control
struct MouseState {
    bool left_button_down = false;
    bool right_button_down = false;
    bool middle_button_down = false;
    int last_x = 0;
    int last_y = 0;
    int delta_x = 0;
    int delta_y = 0;
    float scroll_delta = 0.0f;
    
    void update(int x, int y) {
        delta_x = x - last_x;
        delta_y = y - last_y;
        last_x = x;
        last_y = y;
    }
    
    void reset_deltas() {
        delta_x = 0;
        delta_y = 0;
        scroll_delta = 0.0f;
    }
};

// Simple Windows OpenGL context manager (without GLFW dependency)
class GLWindow {
public:
    GLWindow() : hwnd_(nullptr), hdc_(nullptr), hglrc_(nullptr) {}
    ~GLWindow() { cleanup(); }
    
    bool create(const std::string& title, int width = 1200, int height = 800);
    void swap_buffers();
    bool should_close();
    void poll_events();
    void cleanup();
    
    // Mouse input access
    const MouseState& get_mouse_state() const { return mouse_state_; }
    
private:
    HWND hwnd_;
    HDC hdc_;
    HGLRC hglrc_;
    bool should_close_ = false;
    MouseState mouse_state_;
    
    static LRESULT CALLBACK window_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);
};
