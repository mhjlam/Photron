#pragma once

#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <concepts>
#include <ranges>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "renderer/camera.hpp"
#include "renderer/settings.hpp"

class Simulator;

class Renderer
{
public:
	Renderer();
	~Renderer();

	bool initialize();
	void update();
	void render(Simulator& simulator);
	void set_viewport(int width, int height);

	// Input handling
	void handle_key_input(int key, int scancode, int action, int mods);
	void handle_mouse_move(float xpos, float ypos);
	void handle_mouse_button(int button, int action, int mods);
	void handle_mouse_scroll(float xoffset, float yoffset);

	// Camera control
	void reset_camera();
	void set_camera_mode(bool is_arc_mode);
	bool should_capture_mouse() const;
	Camera& get_camera() { return camera_; }

	// Callback for camera mode changes - Modern C++20 with perfect forwarding
	template<typename Callable>
		requires std::invocable<Callable, bool>
	void set_camera_mode_change_callback(Callable&& callback) {
		camera_mode_change_callback_ = std::forward<Callable>(callback);
	}

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
	void add_triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3,
					  const glm::vec4& color = glm::vec4(1.0f));
	void end_triangles();
	void draw_triangles();
	void draw_triangles_with_clipping(const std::vector<glm::vec4>& clipping_planes);

	// Instanced voxel rendering (high performance)
	void begin_voxel_instances();
	void add_voxel_instance(const glm::vec3& position, const glm::vec4& color, float scale = 1.0f, float depth = 0.0f);
	void end_voxel_instances();
	void draw_voxel_instances();

	// Text billboard rendering for energy labels
	void draw_labels(const Settings& settings);
	void cache_energy_labels();                  // Cache energy labels from current simulation data
	void invalidate_energy_label_cache();        // Force recalculation of energy labels
	void update_energy_label_screen_positions(); // Update screen positions once per frame
	void auto_manage_energy_labels(Settings& settings); // Auto-disable energy labels for performance

	// World-to-screen coordinate conversion for text rendering
	glm::vec2 world_to_screen(const glm::vec3& world_pos, int screen_width, int screen_height) const;

	// Text rendering callback - Modern C++20 with concepts for type safety
	template<typename Callable>
		requires std::invocable<Callable, const std::string&, float, float, const glm::vec4&>
	void set_text_render_callback(Callable&& callback) {
		text_render_callback_ = std::forward<Callable>(callback);
	}

	// Adaptive color mapping helper
	glm::vec4 get_adaptive_energy_color(float energy, float min_energy, float max_energy);
	glm::vec4 get_layer_specific_energy_color(float energy, float min_energy, float max_energy, uint8_t tissue_id);

private:
	void setup_opengl();
	void update_camera();
	void update_camera_target(const Simulator& simulator);

	// Shader-based drawing functions
	void draw_volume(const Simulator& simulator);
	void draw_coordinate_axes();
	void draw_test_geometry();
	void draw_voxels(const Settings& settings);
	void draw_voxels_instanced(const Settings& settings); // High-performance instanced version
	void draw_paths(const Settings& settings);
	void draw_paths_instanced(const Settings& settings); // High-performance instanced version
	void draw_emitters(const Settings& settings);  // Draw true exit points
	
	// Utility methods
	bool is_point_inside_mesh(const glm::vec3& point, const Simulator& simulator) const;

	// Shader management methods
	bool setup_line_rendering();
	bool setup_point_rendering();
	bool setup_triangle_rendering();
	bool setup_voxel_instanced_rendering();
	bool setup_line_instanced_rendering();
	bool setup_point_instanced_rendering();
	
	// Performance optimization: Cache energy range calculation
	void update_cached_energy_range(const Settings& settings) const;
	void invalidate_energy_cache() const;
	
	std::string load_shader_source(const std::string& file_path);
	GLuint compile_shader(const std::string& source, GLenum shader_type);
	GLuint create_shader_program(const std::string& vertex_source, const std::string& fragment_source);

	// Vertex structures for consolidated rendering
	struct LineVertex
	{
		glm::vec3 position{};
		glm::vec4 color{1.0f};
	};

	struct PointVertex
	{
		glm::vec3 position{};
		glm::vec4 color{1.0f};
	};

	struct TriangleVertex
	{
		glm::vec3 position{};
		glm::vec4 color{1.0f};
	};

	// Instance data structure for voxel rendering
	struct VoxelInstance
	{
		glm::vec3 position{};
		glm::vec4 color{1.0f};
		float scale{1.0f};
		float depth{0.0f}; // For depth sorting
	};

	// Instance data structure for line rendering
	struct LineInstance
	{
		glm::vec3 start{};
		glm::vec3 end{};
		glm::vec4 color{1.0f};
	};

	// Instance data structure for point rendering
	struct PointInstance
	{
		glm::vec3 position{};
		glm::vec4 color{1.0f};
		float size{1.0f};
	};

	// Energy label structure for billboard text rendering
	struct EnergyLabel
	{
		glm::vec3 world_position;
		std::string text;
		glm::vec4 color;
		float scale {1.0f};
		glm::vec2 screen_position;          // Cached screen position
		bool screen_position_valid {false}; // Whether screen position is current
	};

	// OpenGL resources for shader-based rendering
	GLuint line_vao_ {0}, line_vbo_ {0};
	GLuint point_vao_ {0}, point_vbo_ {0};
	GLuint triangle_vao_ {0}, triangle_vbo_ {0};
	GLuint line_shader_program_ {0}, point_shader_program_ {0}, triangle_shader_program_ {0};

	// Instanced voxel rendering resources
	GLuint voxel_cube_vao_ {0}, voxel_cube_vbo_ {0}, voxel_instance_vbo_ {0};
	GLuint voxel_shader_program_ {0};
	std::vector<VoxelInstance> voxel_instances_;

	// Instanced line rendering resources
	GLuint line_instanced_vao_ {0}, line_instanced_vbo_ {0}, line_instance_vbo_ {0};
	GLuint line_instanced_shader_program_ {0};
	std::vector<LineInstance> line_instances_;

	// Instanced point rendering resources
	GLuint point_instanced_vao_ {0}, point_instanced_vbo_ {0}, point_instance_vbo_ {0};
	GLuint point_instanced_shader_program_ {0};
	std::vector<PointInstance> point_instances_;

	// Vertex data containers
	std::vector<LineVertex> line_vertices_;
	std::vector<PointVertex> point_vertices_;
	std::vector<TriangleVertex> triangle_vertices_;

	// Camera system
	Camera camera_;
	bool is_arc_camera_mode_ {true}; // true for Orbit mode, false for Free mode

	Settings settings_;
	int viewport_width_ {800};
	int viewport_height_ {600};

	// Text rendering callback for ImGui-based text display
	std::function<void(const std::string&, float, float, const glm::vec4&)> text_render_callback_;

	// Cached energy labels (computed once when simulation updates)
	std::vector<EnergyLabel> cached_energy_labels_;
	bool energy_labels_cached_ {false};

	// Performance optimization: Cache energy analysis results
	mutable float cached_min_energy_ {1.0f};
	mutable float cached_max_energy_ {0.0f};
	mutable bool energy_range_cached_ {false};
	mutable size_t last_path_count_ {0};
	mutable size_t last_voxel_data_version_ {0};

	// Camera state tracking for label position updates
	glm::vec3 last_camera_position_;
	float last_camera_distance_;
	float last_camera_azimuth_;
	float last_camera_elevation_;
	bool camera_state_changed_ {true};

	// Key state tracking for smooth movement
	struct KeyState
	{
		bool w_pressed {false};
		bool a_pressed {false};
		bool s_pressed {false};
		bool d_pressed {false};
		bool q_pressed {false};
		bool e_pressed {false};
	} key_state_;

	// Callback for camera mode changes
	std::function<void(bool)> camera_mode_change_callback_;

	// Keep reference to current simulator
	Simulator* simulator_ {nullptr};
};
