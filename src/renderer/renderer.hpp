#pragma once

#include "../simulator/simulator.hpp"
#include "../structs/camera.hpp"
#include "../structs/keyboard.hpp"
#include "../structs/mouse.hpp"
#include "../structs/settings.hpp"
#include "../structs/window.hpp"

// Simple OpenGL headers
#include <Windows.h>
#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <string>
#include <vector>
#include <memory>

// Forward declarations
struct GLFWwindow;
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

	// Legacy counters (keep for compatibility)
	static unsigned long numvoxall;             // number of voxels (all)
	static unsigned long numvoxrec;             // number of voxels (recorded)

private:
	static GLFWwindow* glfwWindow;              // Legacy GLFW handle (unused)
	static Mouse mouse;
	static Camera camera;
	static WindowDimensions window;
	static Keyboard keyboard;
	static Settings settings;
	static Simulator* simulator;

private:
	// Simple render functions
	static void SetupViewport(int width, int height);
	
	// Simple draw functions
	static void DrawBounds();
	static void DrawLayers();
	static void DrawVoxels();
	static void DrawPhotons();

	// Helper functions
	static void UpdateCamera();
	static void HandleMouseInput();

public:
	// Simple initialization
	static void Initialize();
	static void Initialize(int argc, char* argv[]);
	static void Render(Simulator* sim);
	
	// Matrix functions
	static void SetupPerspective(int width, int height, float fov, float near_plane, float far_plane);
	static void SetViewMatrix(const glm::vec3& eye, const glm::vec3& center, const glm::vec3& up);
	static void UpdateMVP();
};
