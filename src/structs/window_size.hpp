#pragma once

#include <cstdint>

struct WindowSize {
	uint32_t width;
	uint32_t height;

	WindowSize() {
		width = 1280;
		height = 720;
	}

	WindowSize(uint32_t w, uint32_t h) {
		width = static_cast<uint32_t>(w > 0 ? w : 1280);
		height = static_cast<uint32_t>(h > 0 ? h : 720);
	}

	void update(uint32_t w, uint32_t h) {
		width = static_cast<uint32_t>(w);
		height = static_cast<uint32_t>(h);
	}
};
