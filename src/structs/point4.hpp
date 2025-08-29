#pragma once

struct Point4 {
	double x, y, z, w;

	Point4() { x = y = z = w = 0; }

	Point4(double xyzw) { x = y = z = w = xyzw; }

	Point4(double px, double py, double pz, double pw) {
		x = px;
		y = py;
		z = pz;
		w = pw;
	}

	bool operator==(Point4 p) { return (x == p.x && y == p.y && z == p.z && w == p.w); }
};
