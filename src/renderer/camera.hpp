/**
 * @file camera.hpp
 * @brief Dual-mode 3D camera system for interactive visualization
 *
 * Provides both orbital (arc-ball) and free-flight camera modes for flexible
 * 3D scene navigation. The camera system is optimized for scientific visualization
 * with features like automatic framing, elevation constraints, and smooth transitions.
 */

#pragma once

#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

/**
 * @class Camera
 * @brief Dual-mode 3D camera with orbital and free-flight navigation
 *
 * The Camera class provides two distinct interaction paradigms:
 *
 * **Orbital Mode (Arc-ball)**:
 * - Camera orbits around a fixed target point
 * - Uses spherical coordinates (distance, azimuth, elevation)
 * - Ideal for examining objects from all angles
 * - Prevents camera from "getting lost" in 3D space
 *
 * **Free-flight Mode (FPS)**:
 * - First-person style camera movement
 * - Uses Euler angles (yaw, pitch) for orientation
 * - Allows arbitrary positioning and orientation
 * - Better for navigating through complex scenes
 *
 * Key features:
 * - Smooth interpolation between modes
 * - Automatic scene framing capabilities
 * - Elevation constraints for scientific data visualization
 * - Reset functionality to return to initial viewpoint
 * - Optimized matrix calculations with caching
 */
class Camera
{
private:
	glm::vec3 position_; ///< Current camera position in world space
	glm::vec3 target_;   ///< Current camera target point
	glm::vec3 up_;       ///< Current camera up vector

	bool is_fps_mode_;   ///< Camera mode flag: true for FPS, false for orbital

	// FPS camera vectors
	glm::vec3 front_;    ///< FPS mode forward direction vector
	glm::vec3 right_;    ///< FPS mode right direction vector
	glm::vec3 world_up_; ///< World up vector for FPS calculations

	// FPS camera angles
	float yaw_;   ///< Horizontal look angle for FPS mode
	float pitch_; ///< Vertical look angle for FPS mode

	// Spherical coordinates for orbital mode
	float distance_;  ///< Distance from camera to target in orbital mode
	float azimuth_;   ///< Horizontal rotation angle in degrees
	float elevation_; ///< Vertical rotation angle in degrees

	// Transformation matrices
	glm::mat4 view_matrix_;       ///< Current view transformation matrix
	glm::mat4 projection_matrix_; ///< Current projection transformation matrix

	// Camera parameters
	float fov_;          ///< Field of view in degrees
	float near_plane_;   ///< Near clipping plane distance
	float far_plane_;    ///< Far clipping plane distance
	float aspect_ratio_; ///< Viewport aspect ratio (width/height)

	// Elevation constraints
	float min_elevation_;           ///< Minimum allowed elevation angle
	float max_elevation_;           ///< Maximum allowed elevation angle
	bool elevation_bounds_set_;     ///< Flag indicating if elevation bounds are active

	float camera_y_offset_;         ///< Vertical offset for world space movement

	struct InitialCameraState
	{
		float distance {6.0f};      ///< Default camera distance from target
		float azimuth {45.0f};      ///< Default horizontal rotation angle
		float elevation {30.0f};    ///< Default vertical rotation angle
		glm::vec3 target {0.0f};    ///< Default camera target point
	} initial_state_;               ///< Stored initial camera state for reset functionality

	struct MouseState
	{
		bool left_pressed {false};  ///< Left mouse button state
		bool right_pressed {false}; ///< Right mouse button state
		float last_x {0.0f};        ///< Last recorded mouse X position
		float last_y {0.0f};        ///< Last recorded mouse Y position
		bool first_mouse {true};    ///< Flag for first mouse movement detection
		float center_x {0.0f};      ///< Center X position for FPS mode mouse capture
		float center_y {0.0f};      ///< Center Y position for FPS mode mouse capture
		bool center_set {false};    ///< Flag indicating if center position is set
	} mouse_state_;                 ///< Current mouse interaction state

public:
	Camera();

	// Matrix getters
	const glm::mat4& get_view_matrix() const { return view_matrix_; }
	const glm::mat4& get_projection_matrix() const { return projection_matrix_; }
	glm::mat4 get_mvp_matrix(const glm::mat4& model = glm::mat4(1.0f)) const {
		return projection_matrix_ * view_matrix_ * model;
	}

	/**
	 * @brief Get current camera position in world coordinates
	 * @return const glm::vec3& Camera position vector
	 */
	const glm::vec3& get_position() const { return position_; }

	/**
	 * @brief Get current camera target point
	 * @return const glm::vec3& Target position vector
	 */
	const glm::vec3& get_target() const { return target_; }

	/**
	 * @brief Get current camera up direction vector
	 * @return const glm::vec3& Up direction vector
	 */
	const glm::vec3& get_up() const { return up_; }

	/**
	 * @brief Get orbital distance from camera to target
	 * @return float Distance in world units
	 */
	float get_distance() const { return distance_; }

	/**
	 * @brief Get horizontal rotation angle in orbital mode
	 * @return float Azimuth angle in degrees
	 */
	float get_azimuth() const { return azimuth_; }

	/**
	 * @brief Get vertical rotation angle in orbital mode
	 * @return float Elevation angle in degrees
	 */
	float get_elevation() const { return elevation_; }

	/**
	 * @brief Rotate camera around target in orbital mode
	 * @param delta_azimuth Horizontal rotation change in degrees
	 * @param delta_elevation Vertical rotation change in degrees
	 */
	void orbit(float delta_azimuth, float delta_elevation);

	/**
	 * @brief Adjust camera distance from target (zoom in/out)
	 * @param delta Distance change (negative = zoom in, positive = zoom out)
	 */
	void zoom(float delta);

	/**
	 * @brief Set new camera target point for orbital mode
	 * @param target New target position in world coordinates
	 */
	void set_target(const glm::vec3& target);

	/**
	 * @brief Adjust vertical position constraints for camera movement
	 * @param delta Elevation adjustment amount
	 */
	void adjust_camera_elevation(float delta);

	/**
	 * @brief Set vertical movement bounds for camera positioning
	 * @param min_y Minimum allowed Y coordinate
	 * @param max_y Maximum allowed Y coordinate
	 */
	void set_elevation_bounds(float min_y, float max_y);

	/**
	 * @brief Set orbital distance from camera to target
	 * @param distance New distance in world units
	 */
	void set_distance(float distance);

	/**
	 * @brief Set orbital angles directly
	 * @param azimuth Horizontal angle in degrees
	 * @param elevation Vertical angle in degrees
	 */
	void set_angles(float azimuth, float elevation);

	/**
	 * @brief Move camera forward in current facing direction (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_forward(float distance);

	/**
	 * @brief Move camera backward from current facing direction (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_backward(float distance);

	/**
	 * @brief Move camera left relative to current facing direction (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_left(float distance);

	/**
	 * @brief Move camera right relative to current facing direction (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_right(float distance);

	/**
	 * @brief Move camera vertically upward (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_up(float distance);

	/**
	 * @brief Move camera vertically downward (FPS mode)
	 * @param distance Movement distance in world units
	 */
	void move_down(float distance);

	/**
	 * @brief Switch between FPS and orbital camera modes
	 * @param fps_mode True for FPS mode, false for orbital mode
	 */
	void set_fps_mode(bool fps_mode);

	/**
	 * @brief Check if camera is in FPS mode
	 * @return bool True if in FPS mode, false if in orbital mode
	 */
	bool is_fps_mode() const { return is_fps_mode_; }

	/**
	 * @brief Reset camera to initial position and orientation
	 */
	void reset_to_initial();

	/**
	 * @brief Process mouse movement for camera control
	 * @param xpos Current mouse X coordinate
	 * @param ypos Current mouse Y coordinate
	 */
	void handle_mouse_move(float xpos, float ypos);

	/**
	 * @brief Process mouse button press/release events
	 * @param button Mouse button identifier (GLFW constants)
	 * @param action Button action (press/release)
	 */
	void handle_mouse_button(int button, int action);

	/**
	 * @brief Process mouse scroll wheel events for zoom control
	 * @param yoffset Scroll wheel offset (positive = scroll up)
	 */
	void handle_mouse_scroll(float yoffset);

	/**
	 * @brief Check if mouse cursor should be captured and hidden
	 * @return bool True if mouse should be captured for FPS camera control
	 */
	bool should_capture_mouse() const { return is_fps_mode_ && mouse_state_.right_pressed; }

	/**
	 * @brief Set center position for mouse capture in FPS mode
	 * @param center_x Screen center X coordinate
	 * @param center_y Screen center Y coordinate
	 */
	void set_mouse_center(float center_x, float center_y);

	/**
	 * @brief Configure camera projection parameters
	 * @param fov Field of view angle in degrees
	 * @param aspect Aspect ratio (width/height)
	 * @param near_plane Near clipping plane distance
	 * @param far_plane Far clipping plane distance
	 */
	void set_perspective(float fov, float aspect, float near_plane, float far_plane);

	/**
	 * @brief Update camera aspect ratio for viewport changes
	 * @param aspect New aspect ratio (width/height)
	 */
	void set_aspect_ratio(float aspect);

private:
	void update_position_from_spherical();
	void update_view_matrix();
	void update_projection_matrix();
	void update_fps_vectors();
};
