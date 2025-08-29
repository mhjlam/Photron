#pragma once

#include "point3.hpp"

struct Segment {
	Point3 origin;
	Point3 target;

	Segment() {
		origin = Point3(0.0f);
		target = Point3(1.0f);
	}

	Segment(Point3 p0, Point3 p1) {
		origin = p0;
		target = p1;
	}
};
