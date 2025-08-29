#pragma once

#include "point3.hpp"
#include "vector3.hpp"

struct Ray {
	Point3 origin;
	Vector3 direction;

	Ray(Point3 p, Vector3 v, bool normalize = false) : origin(p), direction(v) {
		if (normalize) {
			direction.normalize();
		}
	}
};
