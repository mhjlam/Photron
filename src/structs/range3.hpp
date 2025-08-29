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

	double xmin, xmax; // left and right boundaries
	double ymin, ymax; // bottom and top boundaries
	double zmin, zmax; // rear and front boundaries

	double width;
	double height;
	double depth;

	Range3() {
		xmin = xmax = 0;
		ymin = ymax = 0;
		zmin = zmax = 0;

		width = height = depth = 0;
	}

	Range3(double xx0, double xx1, double yy0, double yy1, double zz0, double zz1) {
		xmin = xx0;
		xmax = xx1;
		ymin = yy0;
		ymax = yy1;
		zmin = zz0;
		zmax = zz1;

		width = fabs(xmin - xmax);
		height = fabs(ymin - ymax);
		depth = fabs(zmin - zmax);
	}

	bool includes(double x, double y, double z) {
		if (x >= xmin && y >= ymin && z >= zmin && x <= xmax && y <= ymax && z <= zmax)
			return true;

		return false;
	}
};
