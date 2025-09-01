#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "structs/camera.hpp"
#include "structs/settings.hpp"

#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

class Simulator;

class Renderer {
public:
    Renderer();
    ~Renderer();

    bool initialize();
    void update();
    void render(Simulator* simulator);
    void set_viewport(int width, int height);

    // Input handling
    void handle_key_input(int key, int scancode, int action, int mods);
    void handle_mouse_move(float xpos, float ypos);
    void handle_mouse_button(int button, int action, int mods);
    void handle_mouse_scroll(float xoffset, float yoffset);

    // Settings
    void set_settings(const Settings& settings) { settings_ = settings; }
    const Settings& get_settings() const { return settings_; }

    // Consolidated shader-based rendering methods
    void begin_lines();
    void add_line(const glm::vec3& start, const glm::vec3& end, const glm::vec4& color = glm::vec4(1.0f));
    void end_lines();
    void draw_lines();

    void begin_points();
    void add_point(const glm::vec3& position, const glm::vec4& color = glm::vec4(1.0f));
    void end_points();
    void draw_points();

    void begin_triangles();
    void add_triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec4& color = glm::vec4(1.0f));
    void end_triangles();
    void draw_triangles();

    // Text billboard rendering for energy labels
    void draw_energy_labels(const Settings& settings);
    
    // Adaptive color mapping helper
    glm::vec4 get_adaptive_energy_color(float energy, float min_energy, float max_energy);

private:
    void setup_opengl();
    void update_camera();
    void update_camera_target(Simulator* simulator);
    
    // Modern shader-based drawing functions
    void draw_bounds_modern(Simulator* simulator);
    void draw_coordinate_axes();
    void draw_test_geometry();
    void draw_voxels_modern(const Settings& settings);
    void draw_paths_modern(const Settings& settings);  
    void draw_photons_modern(const Settings& settings);
    void draw_layers_modern(Simulator* simulator);

    // Shader management methods
    bool setup_line_rendering();
    bool setup_point_rendering();
    bool setup_triangle_rendering();
    std::string load_shader_source(const std::string& file_path);
    GLuint compile_shader(const std::string& source, GLenum shader_type);
    GLuint create_shader_program(const std::string& vertex_source, const std::string& fragment_source);

    // Vertex structures for consolidated rendering
    struct LineVertex {
        glm::vec3 position;
        glm::vec4 color;
    };

    struct PointVertex {
        glm::vec3 position;
        glm::vec4 color;
    };

    struct TriangleVertex {
        glm::vec3 position;
        glm::vec4 color;
    };

    // Energy label structure for billboard text rendering
    struct EnergyLabel {
        glm::vec3 world_position;
        std::string text;
        glm::vec4 color;
        float scale = 1.0f;
    };

    // OpenGL resources for shader-based rendering
    GLuint line_vao_, line_vbo_;
    GLuint point_vao_, point_vbo_;
    GLuint triangle_vao_, triangle_vbo_;
    GLuint line_shader_program_, point_shader_program_, triangle_shader_program_;
    
    // Vertex data containers
    std::vector<LineVertex> line_vertices_;
    std::vector<PointVertex> point_vertices_;
    std::vector<TriangleVertex> triangle_vertices_;

    // Matrix system
    glm::mat4 projection_matrix_;
    glm::mat4 view_matrix_;
    glm::mat4 model_matrix_;
    glm::mat4 mvp_matrix_;

    // Camera controls
    float camera_distance_;
    float camera_rotation_x_;
    float camera_rotation_y_;
    float camera_target_x_;
    float camera_target_y_;
    float camera_target_z_;

    // Mouse state
    struct MouseState {
        bool left_pressed = false;
        bool right_pressed = false;
        bool middle_pressed = false;
        float last_x = 0.0f;
        float last_y = 0.0f;
        bool first_mouse = true;
    } mouse_state_;

    Settings settings_;
    int viewport_width_;
    int viewport_height_;
    
    // Keep reference to current simulator
    Simulator* simulator_;
};
