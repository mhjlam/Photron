#pragma once

#include "glm_types.hpp"

struct Ray
{
	glm::dvec3 origin;
	glm::dvec3 direction;

	Ray(glm::dvec3 p, glm::dvec3 v, bool normalize = false) : origin(p), direction(v) {
		if (normalize) {
			direction = glm::normalize(direction);
		}
	}
};
