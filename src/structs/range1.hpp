#pragma once

#include <cmath>

struct Range1 {
	double x_min, x_max;
	double width;

	Range1() {
		x_min = x_max = 0;
		width = 0;
	}

	Range1(double xx0, double xx1) {
		x_min = xx0;
		x_max = xx1;
		width = fabs(x_max - x_min);
	}

	bool includes(double x) { return (x >= x_min && x <= x_max); }
};
