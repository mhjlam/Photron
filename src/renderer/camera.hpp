#pragma once

#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Camera
{
private:
	glm::vec3 position_;
	glm::vec3 target_;
	glm::vec3 up_;

	// Camera mode
	bool is_fps_mode_; // true for FPS, false for orbital

	// FPS camera vectors
	glm::vec3 front_;
	glm::vec3 right_;
	glm::vec3 world_up_;

	// FPS camera angles (different from orbital angles)
	float yaw_;   // Horizontal look angle
	float pitch_; // Vertical look angle

	// Spherical coordinates for orbiting
	float distance_;
	float azimuth_;   // Horizontal rotation (degrees)
	float elevation_; // Vertical rotation (degrees)

	// Matrices
	glm::mat4 view_matrix_;
	glm::mat4 projection_matrix_;

	// Camera parameters
	float fov_; // Field of view in degrees
	float near_plane_;
	float far_plane_;
	float aspect_ratio_;

	// Tissue elevation bounds for constraining Y movement
	float min_elevation_;
	float max_elevation_;
	bool elevation_bounds_set_;

	// Camera Y-offset for world space vertical movement
	float camera_y_offset_;

	// Initial camera state for reset functionality
	struct InitialCameraState
	{
		float distance = 6.0f; // Closer default distance
		float azimuth = 45.0f;
		float elevation = 30.0f;
		glm::vec3 target{0.0f};
	} initial_state_;

	// Mouse state for camera control
	struct MouseState
	{
		bool left_pressed = false;
		bool right_pressed = false;
		float last_x = 0.0f;
		float last_y = 0.0f;
		bool first_mouse = true;
		// Track center position for FPS mode mouse capture
		float center_x = 0.0f;
		float center_y = 0.0f;
		bool center_set = false;
	} mouse_state_;

public:
	Camera();

	// Matrix getters
	const glm::mat4& get_view_matrix() const { return view_matrix_; }
	const glm::mat4& get_projection_matrix() const { return projection_matrix_; }
	glm::mat4 get_mvp_matrix(const glm::mat4& model = glm::mat4(1.0f)) const {
		return projection_matrix_ * view_matrix_ * model;
	}

	// Camera state getters
	const glm::vec3& get_position() const { return position_; }
	const glm::vec3& get_target() const { return target_; }
	const glm::vec3& get_up() const { return up_; }
	float get_distance() const { return distance_; }
	float get_azimuth() const { return azimuth_; }
	float get_elevation() const { return elevation_; }

	// Camera controls
	void orbit(float delta_azimuth, float delta_elevation);
	void zoom(float delta);
	void set_target(const glm::vec3& target);
	void adjust_camera_elevation(float delta);
	void set_elevation_bounds(float min_y, float max_y);
	void set_distance(float distance);
	void set_angles(float azimuth, float elevation);

	// FPS camera controls
	void move_forward(float distance);
	void move_backward(float distance);
	void move_left(float distance);
	void move_right(float distance);
	void move_up(float distance);
	void move_down(float distance);
	void set_fps_mode(bool fps_mode);
	bool is_fps_mode() const { return is_fps_mode_; }

	// Camera reset
	void reset_to_initial();

	// Mouse input handling
	void handle_mouse_move(float xpos, float ypos);
	void handle_mouse_button(int button, int action);
	void handle_mouse_scroll(float yoffset);

	// Mouse capture for FPS mode
	bool should_capture_mouse() const { return is_fps_mode_ && mouse_state_.right_pressed; }

	// Set center position for mouse capture
	void set_mouse_center(float center_x, float center_y);

	// Projection settings
	void set_perspective(float fov, float aspect, float near_plane, float far_plane);
	void set_aspect_ratio(float aspect);

private:
	void update_position_from_spherical();
	void update_view_matrix();
	void update_projection_matrix();
	void update_fps_vectors();
};
