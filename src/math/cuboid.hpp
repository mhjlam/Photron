#pragma once

#include "glm_types.hpp"

// Forward declaration
class Ray;

class Cuboid
{
private:
	glm::dvec3 min_point_;
	glm::dvec3 max_point_;

public:
	Cuboid(double x1, double y1, double z1, double x2, double y2, double z2);
	Cuboid(const glm::dvec3& min_pt, const glm::dvec3& max_pt);

	// Getters
	const glm::dvec3& min_point() const { return min_point_; }
	const glm::dvec3& max_point() const { return max_point_; }

	// Setters
	void set_bounds(const glm::dvec3& min_pt, const glm::dvec3& max_pt);

	// Utility methods
	glm::dvec3 center() const;
	glm::dvec3 size() const;
	bool contains(const glm::dvec3& point) const;
	double volume() const;

	// Intersection methods
	double intersect_ray_internal(const Ray& ray, glm::dvec3& intersection, glm::dvec3& normal) const;
};
