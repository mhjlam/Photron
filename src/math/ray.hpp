#pragma once

#include <vector>

#include "glm_types.hpp"

// Forward declarations
class Triangle;
class Cuboid;
struct Layer;

class Ray
{
private:
	glm::dvec3 origin_;
	glm::dvec3 direction_;

public:
	Ray(const glm::dvec3& origin, const glm::dvec3& direction, bool normalize = false);

	// Getters
	const glm::dvec3& origin() const { return origin_; }
	const glm::dvec3& direction() const { return direction_; }

	// Setters
	void set_origin(const glm::dvec3& origin) { origin_ = origin; }
	void set_direction(const glm::dvec3& direction, bool normalize = false);

	// Intersection methods
	bool intersect_triangle(Triangle& triangle, glm::dvec3& intersection) const;
	std::pair<bool, glm::dvec3> intersect_triangle(Triangle& triangle) const;
	bool intersect_plane(const glm::dvec3& normal, const glm::dvec3& point, glm::dvec3& intersection) const;

	void intersect_triangles(std::vector<Triangle>& triangles, std::vector<glm::dvec3>& intersections) const;
	double intersect_first_triangle(std::vector<Triangle>& triangles, glm::dvec3& intersection,
									Triangle& hit_triangle) const;
	double intersect_first_triangle_from_layers(const std::vector<Layer>& layers, glm::dvec3& intersection,
												Triangle& hit_triangle) const;
	double intersect_cuboid_internal(Cuboid& cuboid, glm::dvec3& intersection, glm::dvec3& normal) const;
};
