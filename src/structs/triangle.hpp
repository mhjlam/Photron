#pragma once

#include "point3.hpp"
#include "vector3.hpp"

struct Triangle {
	Point3 vertex0;
	Point3 vertex1;
	Point3 vertex2;
	Vector3 normal;

	Triangle() {
		vertex0 = 0;
		vertex1 = 0;
		vertex2 = 0;
	}

	Triangle(Point3 v0, Point3 v1, Point3 v2) {
		vertex0 = v0;
		vertex1 = v1;
		vertex2 = v2;
	}

	bool isinvalid() { return (vertex0 == vertex1 || vertex0 == vertex2 || vertex1 == vertex2); }
};
