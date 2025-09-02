#pragma once

#include "glm_types.hpp"

struct Cuboid {
	glm::dvec3 min_point;
	glm::dvec3 max_point;

	// Convenience accessors for backwards compatibility
	double& min_x = min_point.x;
	double& min_y = min_point.y;
	double& min_z = min_point.z;
	double& max_x = max_point.x;
	double& max_y = max_point.y;
	double& max_z = max_point.z;

	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2) {
		min_point = glm::dvec3(x1, y1, z1);
		max_point = glm::dvec3(x2, y2, z2);
	}
};
