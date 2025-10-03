/**
 * @file input_handler.cpp
 * @brief Implementation of centralized input processing system
 */

#include "input_handler.hpp"

#include <iostream>

#include <GLFW/glfw3.h>

#include "renderer/camera.hpp"

InputHandler::InputHandler(bool initial_orbit_mode) : orbit_mode_(initial_orbit_mode) {
	// Initialize all key states to unpressed
	key_state_ = {};
}

void InputHandler::handle_key_input(int key, int scancode, int action, int mods) {
	(void)scancode; // Unused parameter
	(void)mods;     // Unused parameter

	// Check if WASD is pressed while in Orbit mode - switch to FPS mode automatically
	if (orbit_mode_ && action == GLFW_PRESS) {
		if (key == GLFW_KEY_W || key == GLFW_KEY_A || key == GLFW_KEY_S || key == GLFW_KEY_D) {
			set_camera_mode(false); // Switch to FPS mode
		}
	}

	// Track key states for smooth movement (works in press, release, and repeat)
	bool key_pressed = (action != GLFW_RELEASE);

	switch (key) {
		case GLFW_KEY_W: key_state_.w_pressed = key_pressed; break;
		case GLFW_KEY_A: key_state_.a_pressed = key_pressed; break;
		case GLFW_KEY_S: key_state_.s_pressed = key_pressed; break;
		case GLFW_KEY_D: key_state_.d_pressed = key_pressed; break;
		case GLFW_KEY_SPACE: key_state_.space_pressed = key_pressed; break;
		case GLFW_KEY_LEFT_SHIFT:
		case GLFW_KEY_RIGHT_SHIFT: key_state_.shift_pressed = key_pressed; break;
		default:
			// Ignore other keys
			break;
	}
}

bool InputHandler::handle_mouse_move(Camera& camera, float xpos, float ypos) {
	// Forward mouse movement directly to camera
	camera.handle_mouse_move(xpos, ypos);
	// Mouse movement always indicates camera state change for label updates
	return true;
}

void InputHandler::handle_mouse_button(Camera& camera, int button, int action, int mods) {
	(void)mods; // Unused parameter
	// Forward mouse button events directly to camera
	camera.handle_mouse_button(button, action);
}

void InputHandler::handle_mouse_scroll(Camera& camera, float xoffset, float yoffset) {
	(void)xoffset; // Unused - only vertical scroll matters
	// Forward scroll events directly to camera
	camera.handle_mouse_scroll(yoffset);
}

bool InputHandler::update(Camera& camera) {
	bool camera_state_changed = false;

	// Apply smooth movement for FPS mode only
	if (!orbit_mode_) {
		// Process accumulated key states for smooth FPS movement
		if (key_state_.w_pressed) {
			camera.move_forward(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
		if (key_state_.s_pressed) {
			camera.move_backward(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
		if (key_state_.a_pressed) {
			camera.move_left(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
		if (key_state_.d_pressed) {
			camera.move_right(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
		if (key_state_.space_pressed) {
			camera.move_up(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
		if (key_state_.shift_pressed) {
			camera.move_down(MOVEMENT_SPEED);
			camera_state_changed = true;
		}
	}

	return camera_state_changed;
}

void InputHandler::set_camera_mode(bool orbit_mode) {
	if (orbit_mode != orbit_mode_) {
		orbit_mode_ = orbit_mode;
		notify_mode_change(orbit_mode);
	}
}

bool InputHandler::should_capture_mouse(const Camera& camera) const {
	return camera.should_capture_mouse();
}

void InputHandler::reset_camera(Camera& camera) const {
	camera.reset_to_initial();
}

void InputHandler::notify_mode_change(bool new_mode) {
	if (mode_change_callback_) {
		mode_change_callback_(new_mode);
	}
}
