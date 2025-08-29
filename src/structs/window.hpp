#pragma once

struct WindowDimensions {
	unsigned short width;
	unsigned short height;

	WindowDimensions() {
		width = 800;
		height = 600;
	}

	WindowDimensions(short w, short h) {
		width = static_cast<unsigned short>(w > 0 ? w : 800);
		height = static_cast<unsigned short>(h > 0 ? h : 600);
	}

	void Update(int w, int h) {
		width = static_cast<unsigned short>(w);
		height = static_cast<unsigned short>(h);
	}
};
