#pragma once

#include <cmath>

struct Range1 {
	double xmin, xmax;
	double width;

	Range1() {
		xmin = xmax = 0;
		width = 0;
	}

	Range1(double xx0, double xx1) {
		xmin = xx0;
		xmax = xx1;
		width = fabs(xmax - xmin);
	}

	bool includes(double x) { return (x >= xmin && x <= xmax); }
};
