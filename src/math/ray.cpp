#include "ray.hpp"

#include <algorithm>
#include <limits>

#include "cuboid.hpp"
#include "triangle.hpp"
#include "../simulator/layer.hpp"

Ray::Ray(const glm::dvec3& origin, const glm::dvec3& direction, bool normalize) :
	origin_(origin), direction_(direction) {
	if (normalize) {
		direction_ = glm::normalize(direction_);
	}
}

void Ray::set_direction(const glm::dvec3& direction, bool normalize) {
	direction_ = direction;
	if (normalize) {
		direction_ = glm::normalize(direction_);
	}
}

/***********************************************************
 * Ray-Triangle intersection using Möller-Trumbore algorithm
 * Modern implementation with consistent epsilon handling and vectorized operations
 * Reference: http://www.graphics.cornell.edu/pubs/1997/MT97.html
 ***********************************************************/
bool Ray::intersect_triangle(Triangle& triangle, glm::dvec3& intersection) const {
	// Use consistent epsilon throughout the system
	constexpr double EPSILON = 1e-12;
	
	// Get triangle vertices
	const glm::dvec3 vert0 = triangle.v0();
	const glm::dvec3 vert1 = triangle.v1();
	const glm::dvec3 vert2 = triangle.v2();

	// Calculate edge vectors using GLM operations
	const glm::dvec3 edge1 = vert1 - vert0;
	const glm::dvec3 edge2 = vert2 - vert0;

	// Begin calculating determinant - also used to calculate u parameter
	const glm::dvec3 pvec = glm::cross(direction_, edge2);
	const double det = glm::dot(edge1, pvec);

	// If determinant is near zero, ray lies in plane of triangle or ray is parallel to plane
	if (std::abs(det) < EPSILON) {
		return false;
	}

	const double inv_det = 1.0 / det;

	// Calculate distance from vert0 to ray origin
	const glm::dvec3 tvec = origin_ - vert0;

	// Calculate u parameter and test bound
	const double u = glm::dot(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0) {
		return false;
	}

	// Prepare to test v parameter
	const glm::dvec3 qvec = glm::cross(tvec, edge1);

	// Calculate v parameter and test bound
	const double v = glm::dot(direction_, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0) {
		return false;
	}

	// Calculate t, ray intersects triangle
	const double t = glm::dot(edge2, qvec) * inv_det;

	if (t <= EPSILON) { // Use consistent epsilon for self-intersection avoidance
		return false;
	}

	// Calculate intersection point using vectorized operations
	intersection = origin_ + t * direction_;

	return true;
}

std::pair<bool, glm::dvec3> Ray::intersect_triangle(Triangle& triangle) const {
	std::pair<bool, glm::dvec3> result;
	result.first = intersect_triangle(triangle, result.second);
	return result;
}

/***********************************************************
 * Modern Ray-Plane intersection with vectorized operations
 ***********************************************************/
bool Ray::intersect_plane(const glm::dvec3& normal, const glm::dvec3& point, glm::dvec3& intersection) const {
	constexpr double EPSILON = 1e-12;
	
	// Ensure normal is normalized for accurate calculations
	const glm::dvec3 norm = glm::normalize(normal);
	
	// Calculate denominator (ray direction dot normal)
	const double denominator = glm::dot(direction_, norm);

	// Ray is parallel to plane if denominator is near zero
	if (std::abs(denominator) < EPSILON) {
		return false;
	}

	// Calculate distance from ray origin to plane
	const glm::dvec3 origin_to_point = point - origin_;
	const double t = glm::dot(origin_to_point, norm) / denominator;

	// Intersection is behind ray origin
	if (t <= EPSILON) {
		return false;
	}

	// Calculate intersection point
	intersection = origin_ + t * direction_;
	return true;
}

void Ray::intersect_triangles(std::vector<Triangle>& triangles, std::vector<glm::dvec3>& intersections) const {
	intersections.clear();

	for (Triangle& triangle : triangles) {
		auto result = intersect_triangle(triangle);
		if (result.first) {
			intersections.push_back(result.second);
		}
	}
}

double Ray::intersect_first_triangle(std::vector<Triangle>& triangles, glm::dvec3& intersection,
									 Triangle& hit_triangle) const {
	double shortest_dist = std::numeric_limits<double>::max();
	bool found_intersection = false;

	for (Triangle& triangle : triangles) {
		glm::dvec3 temp_intersection;
		if (intersect_triangle(triangle, temp_intersection)) {
			double dist = glm::distance(origin_, temp_intersection);
			if (dist < shortest_dist) {
				shortest_dist = dist;
				intersection = temp_intersection;
				hit_triangle = triangle;
				found_intersection = true;
			}
		}
	}

	return found_intersection ? shortest_dist : -1.0;
}

double Ray::intersect_first_triangle_from_layers(const std::vector<Layer>& layers, glm::dvec3& intersection,
												 Triangle& hit_triangle) const {
	double shortest_dist = std::numeric_limits<double>::max();
	bool found_intersection = false;

	for (const auto& layer : layers) {
		for (const auto& triangle : layer.mesh) {
			glm::dvec3 temp_intersection;
			if (intersect_triangle(const_cast<Triangle&>(triangle), temp_intersection)) {
				double dist = glm::distance(origin_, temp_intersection);
				if (dist < shortest_dist) {
					shortest_dist = dist;
					intersection = temp_intersection;
					hit_triangle = triangle;
					found_intersection = true;
				}
			}
		}
	}

	return found_intersection ? shortest_dist : -1.0;
}

double Ray::intersect_cuboid_internal(Cuboid& cuboid, glm::dvec3& intersection, glm::dvec3& normal) const {
	return cuboid.intersect_ray_internal(*this, intersection, normal);
}
