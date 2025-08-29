#pragma once

#include "point3.hpp"
#include "triangle.hpp"
#include "vector3.hpp"

struct Source {
	unsigned long id;  // identifier
	Point3 orig;       // origin
	Vector3 dir;       // direction of incidence
	Vector3 dirspec;   // direction of specular reflectance

	Point3 inter;      // intersection point
	Triangle intertri; // triangle at intersection point

	Source() {
		id = 0;
		orig = Point3();
		dir = Vector3();

		inter = Point3();
		intertri = Triangle();

		dirspec = Vector3();
	}

	Source(unsigned long i, Point3 p, Vector3 v) {
		id = i;
		orig = p;
		dir = v;

		inter = Point3();
		intertri = Triangle();

		dirspec = Vector3();
	}
};
