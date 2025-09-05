#pragma once

#include "glm_types.hpp"

struct Vertex
{
	glm::dvec3 position;
	double value;

	Vertex* prev = nullptr; // previous internal vertex
	Vertex* next = nullptr; // next internal vertex  
	Vertex* emit = nullptr; // external vertex

	Vertex(double xx, double yy, double zz, double v) noexcept
		: position(xx, yy, zz), value(v) {
	}

	Vertex(const glm::dvec3& pos, double v) noexcept
		: position(pos), value(v) {
	}
};
