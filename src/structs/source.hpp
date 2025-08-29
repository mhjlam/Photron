#pragma once

#include "point3.hpp"
#include "triangle.hpp"
#include "vector3.hpp"

#include <cstdint>

struct Source {
	uint64_t id;                // identifier
	Point3 origin;              // origin
	Vector3 direction;          // direction of incidence
	Vector3 specular_direction; // direction of specular reflectance

	Point3 intersect;           // intersection point
	Triangle triangle;          // triangle at intersection point

	Source() {
		id = 0;
		origin = Point3();
		direction = Vector3();
		intersect = Point3();
		triangle = Triangle();
		specular_direction = Vector3();
	}

	Source(uint64_t i, Point3 p, Vector3 v) {
		id = i;
		origin = p;
		direction = v;
		intersect = Point3();
		triangle = Triangle();
		specular_direction = Vector3();
	}
};
