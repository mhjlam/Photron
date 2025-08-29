#pragma once

#include <cstdint>

struct Mouse {
	int x;
	int y;

	uint32_t lmb;
	uint32_t rmb;

	Mouse() {
		x = y = 0;
		lmb = rmb = 1;
	}

	void update(int nx, int ny) {
		x = nx;
		y = ny;
	}
};
