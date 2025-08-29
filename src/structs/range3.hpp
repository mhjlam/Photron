#pragma once

#include <cmath>

struct Range3 {
	/*
	 * Right-handed coordinate system
	 *
	 *      | Y
	 * 		|
	 * 		|
	 * 		O ----- X
	 * 	   /
	 * 	  /
	 * 	 Z
	 */

	double x_min, x_max; // left and right boundaries
	double y_min, y_max; // bottom and top boundaries
	double z_min, z_max; // rear and front boundaries

	double width;
	double height;
	double depth;

	Range3() {
		x_min = x_max = 0;
		y_min = y_max = 0;
		z_min = z_max = 0;

		width = height = depth = 0;
	}

	Range3(double xx0, double xx1, double yy0, double yy1, double zz0, double zz1) {
		x_min = xx0;
		x_max = xx1;
		y_min = yy0;
		y_max = yy1;
		z_min = zz0;
		z_max = zz1;

		width = std::fabs(x_min - x_max);
		height = std::fabs(y_min - y_max);
		depth = std::fabs(z_min - z_max);
	}

	bool includes(double x, double y, double z) {
		return (x >= x_min && y >= y_min && z >= z_min && x <= x_max && y <= y_max && z <= z_max);
	}
};
