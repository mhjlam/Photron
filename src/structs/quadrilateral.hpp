#pragma once

#include "point3.hpp"
#include "vector3.hpp"

struct Quadrilateral {
	Point3 v0;
	Point3 v1;
	Point3 v2;
	Point3 v3;
	Vector3 normal;

	Quadrilateral(Point3 v0, Point3 v1, Point3 v2, Point3 v3) {
		v0 = v0;
		v1 = v1;
		v2 = v2;
		v3 = v3;
	}

	bool is_axis_aligned() {
		if ((v0.x == v1.x) && (v0.x == v2.x) && (v0.x == v3.x)) {
			return true;
		}
		else if ((v0.y == v1.y) && (v0.y == v2.y) && (v0.y == v3.y)) {
			return true;
		}
		else if ((v0.z == v1.z) && (v0.z == v2.z) && (v0.z == v3.z)) {
			return true;
		}
		return false;
	}
};
