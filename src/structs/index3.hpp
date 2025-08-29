#pragma once

#include <cstdint>

struct Index3 {
	uint32_t a, b, c;

	Index3() { a = b = c = 0; }

	Index3(uint32_t aa, uint32_t bb, uint32_t cc) {
		a = aa;
		b = bb;
		c = cc;
	}
};
