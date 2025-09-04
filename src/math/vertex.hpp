#pragma once

#include "glm_types.hpp"

struct Vertex
{
	glm::dvec3 position;
	double value;

	Vertex* prev; // previous internal vertex
	Vertex* next; // next internal vertex
	Vertex* emit; // external vertex

	Vertex(double xx, double yy, double zz, double v) 
		: position(xx, yy, zz), value(v), prev(nullptr), next(nullptr), emit(nullptr) {
	}

	Vertex(const glm::dvec3& pos, double v) 
		: position(pos), value(v), prev(nullptr), next(nullptr), emit(nullptr) {
	}
};
