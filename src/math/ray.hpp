#pragma once

#include <span>
#include <utility>
#include <vector>

#include "math.hpp"

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
	const glm::dvec3& origin() const noexcept { return origin_; }
	const glm::dvec3& direction() const noexcept { return direction_; }

	// Setters
	void set_origin(const glm::dvec3& origin) noexcept { origin_ = origin; }
	void set_direction(const glm::dvec3& direction, bool normalize = false);

	// Intersection methods - modernized with noexcept and span support
	bool intersect_triangle(Triangle& triangle, glm::dvec3& intersection) const;
	std::pair<bool, glm::dvec3> intersect_triangle(Triangle& triangle) const;
	bool intersect_plane(const glm::dvec3& normal, const glm::dvec3& point, glm::dvec3& intersection) const noexcept;

	// C++20: Use span for better array interface without copying
	void intersect_triangles(std::span<Triangle> triangles, std::vector<glm::dvec3>& intersections) const;
	double intersect_first_triangle(std::span<Triangle> triangles, glm::dvec3& intersection,
									Triangle& hit_triangle) const;
	double intersect_first_triangle_from_layers(std::span<const Layer> layers, glm::dvec3& intersection,
												Triangle& hit_triangle) const;
	double intersect_cuboid_internal(Cuboid& cuboid, glm::dvec3& intersection, glm::dvec3& normal) const;
};
