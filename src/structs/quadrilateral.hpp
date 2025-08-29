#pragma once

#include "point3.hpp"
#include "vector3.hpp"

struct Quadrilateral {
	Point3 vertex0;
	Point3 vertex1;
	Point3 vertex2;
	Point3 vertex3;
	Vector3 normal;

	Quadrilateral(Point3 v0, Point3 v1, Point3 v2, Point3 v3) {
		vertex0 = v0;
		vertex1 = v1;
		vertex2 = v2;
		vertex3 = v3;
	}

	bool isaxisaligned() {
		if ((vertex0.x == vertex1.x) && (vertex0.x == vertex2.x) && (vertex0.x == vertex3.x))
			return true;

		else if ((vertex0.y == vertex1.y) && (vertex0.y == vertex2.y) && (vertex0.y == vertex3.y))
			return true;

		else if ((vertex0.z == vertex1.z) && (vertex0.z == vertex2.z) && (vertex0.z == vertex3.z))
			return true;

		return false;
	}
};
