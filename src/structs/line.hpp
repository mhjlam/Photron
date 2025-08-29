#pragma once

#include "point3.hpp"

#include <cstdint>

struct Line {
	Point3 p1;
	Point3 p2;

	double slope;
	double intercept;

	Line(Point3 pt1, Point3 pt2) {
		p1 = pt1;
		p2 = pt2;

		slope = (p2.y - p1.y) / (p2.x - p1.x);
		intercept = p1.y - (slope * p1.x);
	}
};
