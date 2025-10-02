/**
 * @file input_handler.hpp
 * @brief Input handling system for camera control and mode management
 *
 * Provides centralized input processing for keyboard and mouse events,
 * camera movement, and seamless mode switching between orbital and FPS modes.
 */

#pragma once

#include <concepts>
#include <functional>

// Forward declarations
class Camera;

/**
 * @class InputHandler
 * @brief Centralized input processing for camera control
 *
 * Handles all keyboard and mouse input processing for camera movement
 * and mode switching. Encapsulates key state tracking, automatic mode
 * switching, and provides clean callbacks for UI integration.
 *
 * Features:
 * - WASD key state tracking for smooth FPS movement
 * - Automatic orbit→FPS mode switching on WASD press
 * - Mouse input forwarding to camera
 * - Camera mode change callbacks
 * - Clean separation from rendering logic
 */
class InputHandler
{
public:
	/**
	 * @brief Construct InputHandler with initial camera mode
	 * @param initial_orbit_mode Initial camera mode (true = orbit, false = FPS)
	 */
	explicit InputHandler(bool initial_orbit_mode = true);

	/**
	 * @brief Process keyboard input events
	 *
	 * Handles WASD movement keys, automatic mode switching, and key state tracking.
	 *
	 * @param key GLFW key code
	 * @param scancode Platform-specific scan code
	 * @param action Key action (GLFW_PRESS, GLFW_RELEASE, GLFW_REPEAT)
	 * @param mods Modifier key flags
	 */
	void handle_key_input(int key, int scancode, int action, int mods);

	/**
	 * @brief Process mouse movement events
	 * @param camera Reference to camera for immediate processing
	 * @param xpos Mouse X coordinate
	 * @param ypos Mouse Y coordinate
	 * @return true if camera state changed
	 */
	bool handle_mouse_move(Camera& camera, float xpos, float ypos);

	/**
	 * @brief Process mouse button events
	 * @param camera Reference to camera for immediate processing
	 * @param button Mouse button code
	 * @param action Button action (GLFW_PRESS, GLFW_RELEASE)
	 * @param mods Modifier key flags
	 */
	void handle_mouse_button(Camera& camera, int button, int action, int mods);

	/**
	 * @brief Process mouse scroll events
	 * @param camera Reference to camera for immediate processing
	 * @param xoffset Horizontal scroll offset
	 * @param yoffset Vertical scroll offset
	 */
	void handle_mouse_scroll(Camera& camera, float xoffset, float yoffset);

	/**
	 * @brief Update camera based on current input state
	 *
	 * Processes accumulated key states and applies smooth FPS movement.
	 * Should be called every frame in the update loop.
	 *
	 * @param camera Reference to camera to update
	 * @return true if camera state changed (for cache invalidation)
	 */
	bool update(Camera& camera);

	/**
	 * @brief Get current camera mode
	 * @return true if in orbit mode, false if in FPS mode
	 */
	bool is_orbit_mode() const { return orbit_mode_; }

	/**
	 * @brief Set camera mode explicitly
	 * @param orbit_mode true for orbit mode, false for FPS mode
	 */
	void set_camera_mode(bool orbit_mode);

	/**
	 * @brief Check if mouse should be captured for FPS mode
	 * @param camera Reference to camera for capture check
	 * @return true if mouse should be captured and hidden
	 */
	bool should_capture_mouse(const Camera& camera) const;

	/**
	 * @brief Reset camera to initial state
	 * @param camera Reference to camera to reset
	 */
	void reset_camera(Camera& camera) const;

	/**
	 * @brief Set callback for camera mode change notifications
	 *
	 * Registers a callback function that will be invoked whenever the camera
	 * mode changes (either through automatic WASD switching or explicit setting).
	 *
	 * @tparam Callable Function object type that accepts a boolean parameter
	 * @param callback Function to call when mode changes (true = orbit, false = FPS)
	 */
	template<typename Callable>
		requires std::invocable<Callable, bool>
	void set_camera_mode_change_callback(Callable&& callback) {
		mode_change_callback_ = std::forward<Callable>(callback);
	}

private:
	/**
	 * @struct KeyState
	 * @brief Tracks pressed state of movement keys for smooth FPS control
	 */
	struct KeyState
	{
		bool w_pressed {false};                      ///< W key (forward) pressed state
		bool a_pressed {false};                      ///< A key (left) pressed state
		bool s_pressed {false};                      ///< S key (backward) pressed state
		bool d_pressed {false};                      ///< D key (right) pressed state
		bool space_pressed {false};                  ///< Space key (up) pressed state
		bool shift_pressed {false};                  ///< Shift key (down) pressed state
	};

	bool orbit_mode_ {true};                         ///< Current camera mode
	KeyState key_state_;                             ///< Current key press states
	std::function<void(bool)> mode_change_callback_; ///< Callback for mode changes

	static constexpr float MOVEMENT_SPEED = 0.02f;   ///< FPS movement speed per frame

	/**
	 * @brief Notify callback of mode change
	 * @param new_mode New camera mode
	 */
	void notify_mode_change(bool new_mode);
};
