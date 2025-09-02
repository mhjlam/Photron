#pragma once

#include "glm_types.hpp"
#include <cstdint>

struct WindowSize {
	glm::uvec2 size;

	// Convenience accessors for backwards compatibility
	uint32_t& width = size.x;
	uint32_t& height = size.y;

	WindowSize() {
		size = glm::uvec2(1280, 720);
	}

	WindowSize(uint32_t w, uint32_t h) {
		size.x = w > 0 ? w : 1280;
		size.y = h > 0 ? h : 720;
	}

	void update(uint32_t w, uint32_t h) {
		size = glm::uvec2(w, h);
	}
};
