#pragma once

#include "glm_types.hpp"

struct Triangle {
	glm::dvec3 v0;
	glm::dvec3 v1;
	glm::dvec3 v2;
	glm::dvec3 normal;

	Triangle() {
		v0 = glm::dvec3(0);
		v1 = glm::dvec3(0);
		v2 = glm::dvec3(0);
	}

	Triangle(glm::dvec3 v0, glm::dvec3 v1, glm::dvec3 v2) {
		this->v0 = v0;
		this->v1 = v1;
		this->v2 = v2;
	}

	bool is_invalid() { return (v0 == v1 || v0 == v2 || v1 == v2); }
};
