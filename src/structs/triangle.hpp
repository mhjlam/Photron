#pragma once

#include "point3.hpp"
#include "vector3.hpp"

struct Triangle {
	Point3 v0;
	Point3 v1;
	Point3 v2;
	Vector3 normal;

	Triangle() {
		v0 = 0;
		v1 = 0;
		v2 = 0;
	}

	Triangle(Point3 v0, Point3 v1, Point3 v2) {
		v0 = v0;
		v1 = v1;
		v2 = v2;
	}

	bool is_invalid() { return (v0 == v1 || v0 == v2 || v1 == v2); }
};
