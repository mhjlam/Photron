#pragma once

struct Point3 {
	double x, y, z;

	Point3() { x = y = z = 0; }

	Point3(double xyz) { x = y = z = xyz; }

	Point3(double px, double py, double pz) {
		x = px;
		y = py;
		z = pz;
	}

	bool operator==(Point3& p) { return (x == p.x && y == p.y && z == p.z); }
};
