#pragma once

#include "point3.hpp"

struct Vertex {
	double x, y, z;
	double value;

	Vertex* prev; // previous internal vertex
	Vertex* next; // next internal vertex
	Vertex* emit; // external vertex

	Vertex(double xx, double yy, double zz, double v) {
		x = xx;
		y = yy;
		z = zz;
		value = v;
		prev = 0;
		next = 0;
		emit = 0;
	}

	Vertex(Point3 position, double v) {
		x = position.x;
		y = position.y;
		z = position.z;
		value = v;
		prev = 0;
		next = 0;
		emit = 0;
	}
};
