#pragma once

#include <memory>
#include "glm_types.hpp"

struct Vertex
{
	glm::dvec3 position;
	double value;

	std::shared_ptr<Vertex> prev = nullptr; // previous internal vertex
	std::shared_ptr<Vertex> next = nullptr; // next internal vertex  
	std::shared_ptr<Vertex> emit = nullptr; // external vertex

	Vertex(double xx, double yy, double zz, double v) noexcept
		: position(xx, yy, zz), value(v) {
	}

	Vertex(const glm::dvec3& pos, double v) noexcept
		: position(pos), value(v) {
	}
};
