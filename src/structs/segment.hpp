#pragma once

#include "glm_types.hpp"

struct Segment
{
	glm::dvec3 origin;
	glm::dvec3 target;

	Segment() {
		origin = glm::dvec3(0.0);
		target = glm::dvec3(1.0);
	}

	Segment(glm::dvec3 p0, glm::dvec3 p1) {
		origin = p0;
		target = p1;
	}
};
