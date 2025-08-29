#pragma once

#include <windows.h>

#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "simulator/simulator.hpp"
#include "structs/camera.hpp"
#include "structs/mouse.hpp"
#include "structs/settings.hpp"
#include "structs/window_size.hpp"

#include <memory>
#include <string>
#include <vector>

// Forward declarations
class GLWindow;

class Renderer
{
private:
	// Simple OpenGL rendering
	static std::unique_ptr<GLWindow> gl_window_;

	// Modern matrix system
	static glm::mat4 projection_matrix_;
	static glm::mat4 view_matrix_;
	static glm::mat4 model_matrix_;
	static glm::mat4 mvp_matrix_;

	// Interactive camera controls
	static float camera_distance_;
	static float camera_rotation_x_;
	static float camera_rotation_y_;

	static Mouse mouse_;
	static Camera camera_;
	static WindowSize window_size_;
	static Settings settings_;
	static Simulator* simulator_;

private:
	// Draw functions
	static void draw_bounds();
	static void draw_layers();
	static void draw_voxels();
	static void draw_photons();

	// Helper functions
	static void update_camera();
	static void handle_mouse_input();

public:
	// Simple initialization
	static void initialize();
	static void initialize(int argc, char* argv[]);
	static void render(Simulator* simulator);

	// Matrix functions
	static void setup_perspective(int width, int height, float fov, float near_plane, float far_plane);
	static void set_view_matrix(const glm::vec3& eye, const glm::vec3& center, const glm::vec3& up);
	static void update_mvp();
};
