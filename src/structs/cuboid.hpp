#pragma once

struct Cuboid {
	double min_x, min_y, min_z;
	double max_x, max_y, max_z;

	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) {
		min_x = x1;
		min_y = y1;
		min_z = z1;
		max_x = x2;
		max_y = y2;
		max_z = z2;
	}
};
