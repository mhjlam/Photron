#pragma once

struct Index3 {
	unsigned int a, b, c;

	Index3() { a = b = c = 0; }

	Index3(unsigned int aa, unsigned int bb, unsigned int cc) {
		a = aa;
		b = bb;
		c = cc;
	}
};
