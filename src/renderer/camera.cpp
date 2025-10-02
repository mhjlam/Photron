/**
 * @file camera.cpp
 * @brief Implementation of 3D camera system with orbital and FPS modes
 *
 * Implements a flexible 3D camera system supporting both orbital (arc-ball)
 * and first-person shooter (FPS) camera modes. Features smooth movement,
 * configurable bounds, and efficient matrix calculations for real-time
 * 3D visualization of Monte Carlo simulation results.
 */

#include "camera.hpp"

namespace
{
// Camera system default parameters
constexpr float DEFAULT_DISTANCE = 6.0f;
constexpr float DEFAULT_AZIMUTH = 45.0f;
constexpr float DEFAULT_ELEVATION = 30.0f;
constexpr float DEFAULT_YAW = -90.0f;
constexpr float DEFAULT_PITCH = 0.0f;
constexpr float DEFAULT_FOV = 45.0f;
constexpr float DEFAULT_NEAR_PLANE = 0.1f;
constexpr float DEFAULT_FAR_PLANE = 100.0f;
constexpr float DEFAULT_ASPECT_RATIO = 1.0f;
constexpr float MIN_ELEVATION_DEFAULT = -10.0f;
constexpr float MAX_ELEVATION_DEFAULT = 10.0f;
} // namespace

Camera::Camera() {
	// Initialize camera coordinate system
	target_ = glm::vec3(0.0f);
	up_ = glm::vec3(0.0f, 1.0f, 0.0f);
	world_up_ = glm::vec3(0.0f, 1.0f, 0.0f);

	// Configure orbital camera defaults for isometric viewing
	distance_ = DEFAULT_DISTANCE;   // Optimal distance for simulation visualization
	azimuth_ = DEFAULT_AZIMUTH;     // 45 degrees from front for good perspective
	elevation_ = DEFAULT_ELEVATION; // 30 degrees elevation for isometric view

	// Configure FPS camera defaults
	is_fps_mode_ = false;
	yaw_ = DEFAULT_YAW; // Point towards negative Z axis
	pitch_ = DEFAULT_PITCH;
	position_ = glm::vec3(0.0f, 0.0f, 3.0f);

	// Configure projection parameters
	fov_ = DEFAULT_FOV;
	near_plane_ = DEFAULT_NEAR_PLANE;
	far_plane_ = DEFAULT_FAR_PLANE;
	aspect_ratio_ = DEFAULT_ASPECT_RATIO;

	// Initialize elevation bounds and vertical offset
	min_elevation_ = MIN_ELEVATION_DEFAULT;
	max_elevation_ = MAX_ELEVATION_DEFAULT;
	elevation_bounds_set_ = false;
	camera_y_offset_ = 0.0f;

	// Store initial state for reset functionality
	initial_state_.distance = distance_;
	initial_state_.azimuth = azimuth_;
	initial_state_.elevation = elevation_;
	initial_state_.target = target_;

	// Initialize camera matrices and vectors
	update_fps_vectors();
	update_position_from_spherical();
	update_view_matrix();
	update_projection_matrix();
}

void Camera::orbit(float delta_azimuth, float delta_elevation) {
	// Update orbital camera angles with elevation clamping
	azimuth_ += delta_azimuth;
	elevation_ = std::clamp(elevation_ + delta_elevation, -89.0f, 89.0f);

	// Recalculate camera position and view matrix
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::zoom(float delta) {
	// Adjust camera distance with minimum distance constraint
	distance_ = std::max(0.1f, distance_ - delta);

	// Update camera position for new distance
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::set_target(const glm::vec3& target) {
	// Update camera focus target and refresh view matrix
	target_ = target;
	update_view_matrix();
}

void Camera::adjust_camera_elevation(float delta) {
	// Adjust vertical camera position within configured bounds
	float new_offset = camera_y_offset_ + delta;
	if (elevation_bounds_set_) {
		new_offset = std::clamp(new_offset, min_elevation_, max_elevation_);
	}
	camera_y_offset_ = new_offset;

	// Update camera position and view matrix
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::set_elevation_bounds(float min_y, float max_y) {
	// Configure vertical movement bounds based on scene geometry
	float range = max_y - min_y;
	min_elevation_ = -range; // Allow moving down by material height
	max_elevation_ = range;  // Allow moving up by material height
	elevation_bounds_set_ = true;

	// Clamp current position to new bounds
	camera_y_offset_ = std::clamp(camera_y_offset_, min_elevation_, max_elevation_);
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::set_distance(float distance) {
	distance_ = std::max(0.1f, distance);
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::set_angles(float azimuth, float elevation) {
	azimuth_ = azimuth;
	elevation_ = std::clamp(elevation, -89.0f, 89.0f);
	update_position_from_spherical();
	update_view_matrix();
}

void Camera::set_perspective(float fov, float aspect, float near_plane, float far_plane) {
	fov_ = fov;
	aspect_ratio_ = aspect;
	near_plane_ = near_plane;
	far_plane_ = far_plane;
	update_projection_matrix();
}

void Camera::set_aspect_ratio(float aspect) {
	aspect_ratio_ = aspect;
	update_projection_matrix();
}

void Camera::update_position_from_spherical() {
	float azimuth_rad = glm::radians(azimuth_);
	float elevation_rad = glm::radians(elevation_);

	float x = distance_ * cos(elevation_rad) * cos(azimuth_rad);
	float y = distance_ * sin(elevation_rad);
	float z = distance_ * cos(elevation_rad) * sin(azimuth_rad);

	// Apply Y-offset to move entire camera up/down in world space
	position_ = target_ + glm::vec3(x, y + camera_y_offset_, z);
}

void Camera::update_view_matrix() {
	if (is_fps_mode_) {
		// FPS camera: look in the direction of the front vector
		view_matrix_ = glm::lookAt(position_, position_ + front_, up_);
	}
	else {
		// Orbital camera: look at target with Y-offset
		glm::vec3 adjusted_target = target_ + glm::vec3(0.0f, camera_y_offset_, 0.0f);
		view_matrix_ = glm::lookAt(position_, adjusted_target, up_);
	}
}

void Camera::update_projection_matrix() {
	projection_matrix_ = glm::perspective(glm::radians(fov_), aspect_ratio_, near_plane_, far_plane_);
}

void Camera::handle_mouse_move(float xpos, float ypos) {
	if (mouse_state_.first_mouse) {
		mouse_state_.last_x = xpos;
		mouse_state_.last_y = ypos;
		mouse_state_.first_mouse = false;
	}

	float x_offset = xpos - mouse_state_.last_x;
	float y_offset = mouse_state_.last_y - ypos; // Reversed since y-coordinates go from bottom to top

	mouse_state_.last_x = xpos;
	mouse_state_.last_y = ypos;

	if (is_fps_mode_) {
		// FPS mode: mouse look only when right button is held
		if (mouse_state_.right_pressed) {
			// If we have a center position set and mouse capture is active,
			// calculate offset from center instead of last position
			if (mouse_state_.center_set) {
				x_offset = xpos - mouse_state_.center_x;
				y_offset = mouse_state_.center_y - ypos; // Reversed for proper look
			}

			const float sensitivity = 0.1f;
			x_offset *= sensitivity;
			y_offset *= sensitivity;

			yaw_ += x_offset;
			pitch_ += y_offset;

			// Constrain pitch
			pitch_ = std::clamp(pitch_, -89.0f, 89.0f);

			update_fps_vectors();
			update_view_matrix();
		}
	}
	else {
		// Orbital mode - require mouse button press
		if (mouse_state_.left_pressed) {
			// Convert to camera's orbit system (azimuth and elevation) with high sensitivity
			orbit(x_offset * 0.1f, -y_offset * 0.1f);
		}

		// Right mouse button for zoom
		if (mouse_state_.right_pressed) {
			zoom(y_offset * 0.02f);
		}
	}
}

void Camera::handle_mouse_button(int button, int action) {
	// Assuming GLFW constants: GLFW_MOUSE_BUTTON_LEFT = 0, GLFW_MOUSE_BUTTON_RIGHT = 1, GLFW_PRESS = 1
	if (button == 0) {                              // GLFW_MOUSE_BUTTON_LEFT
		mouse_state_.left_pressed = (action == 1);  // GLFW_PRESS
	}
	if (button == 1) {                              // GLFW_MOUSE_BUTTON_RIGHT
		mouse_state_.right_pressed = (action == 1); // GLFW_PRESS
	}
}

void Camera::handle_mouse_scroll(float y_offset) {
	if (is_fps_mode_) {
		// In FPS mode, scroll adjusts movement speed or does nothing
		// For now, do nothing in FPS mode
		return;
	}
	else {
		// Scroll wheel moves camera physically up/down along world Y-axis in orbital mode
		adjust_camera_elevation(y_offset * 0.05f); // Reduced sensitivity from 0.2f to 0.05f
	}
}

// FPS camera methods
void Camera::move_forward(float distance) {
	if (is_fps_mode_) {
		position_ += front_ * distance;
		update_view_matrix();
	}
}

void Camera::move_backward(float distance) {
	if (is_fps_mode_) {
		position_ -= front_ * distance;
		update_view_matrix();
	}
}

void Camera::move_left(float distance) {
	if (is_fps_mode_) {
		position_ -= right_ * distance;
		update_view_matrix();
	}
}

void Camera::move_right(float distance) {
	if (is_fps_mode_) {
		position_ += right_ * distance;
		update_view_matrix();
	}
}

void Camera::move_up(float distance) {
	if (is_fps_mode_) {
		position_ += world_up_ * distance;
		update_view_matrix();
	}
}

void Camera::move_down(float distance) {
	if (is_fps_mode_) {
		position_ -= world_up_ * distance;
		update_view_matrix();
	}
}

void Camera::set_fps_mode(bool fps_mode) {
	if (fps_mode != is_fps_mode_) {
		is_fps_mode_ = fps_mode;
		if (is_fps_mode_) {
			// Switch to FPS mode: preserve both position AND look direction exactly
			// Calculate the current look direction from orbital camera
			glm::vec3 look_direction = glm::normalize(target_ - position_);

			// Keep the exact same position - no movement!
			// Just convert the look direction to yaw/pitch for FPS camera
			yaw_ = glm::degrees(atan2(look_direction.z, look_direction.x));
			pitch_ = glm::degrees(asin(look_direction.y));

			// Constrain pitch to valid range
			pitch_ = std::clamp(pitch_, -89.0f, 89.0f);

			world_up_ = glm::vec3(0.0f, 1.0f, 0.0f);
			update_fps_vectors();
		}
		else {
			// Switch to orbit mode: preserve current position but make it look at target
			// Calculate orbital parameters based on current FPS position relative to target
			glm::vec3 offset = position_ - target_;

			// Calculate distance (radius of orbit)
			distance_ = glm::length(offset);

			// Calculate azimuth (horizontal angle around Y-axis)
			azimuth_ = glm::degrees(atan2(offset.z, offset.x));

			// Calculate elevation (vertical angle from horizontal plane)
			float horizontal_distance = sqrt(offset.x * offset.x + offset.z * offset.z);
			elevation_ = glm::degrees(atan2(offset.y, horizontal_distance));

			// Reset Y-offset and up vector for clean orbital mode
			camera_y_offset_ = 0.0f;
			up_ = world_up_;
		}

		// Update view matrix to look at target
		update_view_matrix();
	}
}

void Camera::reset_to_initial() {
	// Always reset to initial arc camera position
	distance_ = initial_state_.distance;
	azimuth_ = initial_state_.azimuth;
	elevation_ = initial_state_.elevation;
	target_ = initial_state_.target;
	camera_y_offset_ = 0.0f;

	// Reset camera orientation completely
	up_ = world_up_; // Reset up vector to world up

	// If we're in FPS mode, also reset FPS vectors
	if (is_fps_mode_) {
		yaw_ = -90.0f;
		pitch_ = 0.0f;
		update_fps_vectors();
	}

	update_position_from_spherical();
	update_view_matrix();
}

void Camera::update_fps_vectors() {
	// Calculate front vector from yaw and pitch
	glm::vec3 direction;
	direction.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
	direction.y = sin(glm::radians(pitch_));
	direction.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
	front_ = glm::normalize(direction);

	// Calculate right and up vectors
	right_ = glm::normalize(glm::cross(front_, world_up_));
	up_ = glm::normalize(glm::cross(right_, front_));
}

void Camera::set_mouse_center(float center_x, float center_y) {
	mouse_state_.center_x = center_x;
	mouse_state_.center_y = center_y;
	mouse_state_.center_set = true;
}
