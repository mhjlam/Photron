#pragma once

#include <cstdint>

#include "math/glm_types.hpp"
#include "math/triangle.hpp"

struct Source
{
	uint64_t id;                   // identifier
	glm::dvec3 origin;             // origin
	glm::dvec3 direction;          // direction of incidence
	glm::dvec3 specular_direction; // direction of specular reflectance

	glm::dvec3 intersect;          // intersection point
	Triangle triangle;             // triangle at intersection point

	Source() {
		id = 0;
		origin = glm::dvec3(0);
		direction = glm::dvec3(0);
		intersect = glm::dvec3(0);
		triangle = Triangle();
		specular_direction = glm::dvec3(0);
	}

	Source(uint64_t i, glm::dvec3 p, glm::dvec3 v) {
		id = i;
		origin = p;
		direction = v;
		intersect = glm::dvec3(0);
		triangle = Triangle();
		specular_direction = glm::dvec3(0);
	}
};
