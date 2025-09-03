#pragma once

#include <cmath>

#include "glm_types.hpp"

struct Range3
{
	/*
	 * Right-handed coordinate system
	 *
	 *      | Y
	 * 		|
	 * 		|
	 * 		O ----- X
	 * 	   /
	 * 	  /
	 * 	 Z
	 */

	// Keep original interface but use GLM internally
	double x_min, x_max; // left and right boundaries
	double y_min, y_max; // bottom and top boundaries
	double z_min, z_max; // rear and front boundaries

	double width;
	double height;
	double depth;

	Range3() {
		x_min = x_max = 0;
		y_min = y_max = 0;
		z_min = z_max = 0;
		width = height = depth = 0;
	}

	Range3(double x0, double x1, double y0, double y1, double z0, double z1) {
		x_min = x0;
		x_max = x1;
		y_min = y0;
		y_max = y1;
		z_min = z0;
		z_max = z1;

		width = std::fabs(x_max - x_min);
		height = std::fabs(y_max - y_min);
		depth = std::fabs(z_max - z_min);
	}

	Range3(const glm::dvec3& min_pt, const glm::dvec3& max_pt) {
		x_min = min_pt.x;
		x_max = max_pt.x;
		y_min = min_pt.y;
		y_max = max_pt.y;
		z_min = min_pt.z;
		z_max = max_pt.z;

		width = std::fabs(x_max - x_min);
		height = std::fabs(y_max - y_min);
		depth = std::fabs(z_max - z_min);
	}

	bool includes(double x, double y, double z) {
		return (x >= x_min && y >= y_min && z >= z_min && x <= x_max && y <= y_max && z <= z_max);
	}

	bool includes(const glm::dvec3& point) { return includes(point.x, point.y, point.z); }

	// GLM convenience methods
	glm::dvec3 min_bound() const { return glm::dvec3(x_min, y_min, z_min); }
	glm::dvec3 max_bound() const { return glm::dvec3(x_max, y_max, z_max); }
	glm::dvec3 center() const { return (min_bound() + max_bound()) * 0.5; }
	glm::dvec3 size() const { return glm::dvec3(width, height, depth); }
};
