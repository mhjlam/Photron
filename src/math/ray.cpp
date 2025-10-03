/**
 * @file ray.cpp
 * @brief Ray casting and intersection algorithms for 3D geometry
 *
 * Implements high-performance ray-geometry intersection tests using optimized algorithms:
 * - Möller-Trumbore ray-triangle intersection for mesh geometry
 * - Vectorized ray-plane intersection for planar surfaces
 * - Multi-triangle intersection with distance sorting
 * - Layer-based intersection for multi-material geometry
 *
 * All algorithms use consistent epsilon handling for numerical stability
 * and GLM vectorized operations for performance optimization.
 */

#include "ray.hpp"

#include <algorithm>
#include <limits>

#include "math/cuboid.hpp"
#include "math/triangle.hpp"
#include "simulator/layer.hpp"

/**
 * @brief Construct ray with origin and direction
 * @param origin Starting point of the ray in 3D space
 * @param direction Direction vector of the ray
 * @param normalize Whether to normalize the direction vector
 *
 * Creates a ray for geometric intersection testing. Direction normalization
 * is optional but recommended for consistent distance calculations.
 */
Ray::Ray(const glm::dvec3& origin, const glm::dvec3& direction, bool normalize) :
	origin_(origin), direction_(direction) {
	// Normalize direction for consistent distance calculations
	if (normalize) {
		direction_ = glm::normalize(direction_);
	}
}

/**
 * @brief Update ray direction with optional normalization
 * @param direction New direction vector
 * @param normalize Whether to normalize the direction vector
 *
 * Allows dynamic ray direction updates during ray casting operations.
 */
void Ray::set_direction(const glm::dvec3& direction, bool normalize) {
	direction_ = direction;
	// Apply normalization if requested for consistent behavior
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
	if (std::abs(det) < MathConstants::GEOMETRIC_EPSILON) {
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

	if (t <= MathConstants::GEOMETRIC_EPSILON) { // Use consistent epsilon for self-intersection avoidance
		return false;
	}

	// Calculate intersection point using vectorized operations
	intersection = origin_ + t * direction_;

	return true;
}

/**
 * @brief Convenient wrapper for ray-triangle intersection
 * @param triangle Target triangle for intersection test
 * @return Pair containing intersection success and point
 *
 * Provides a more convenient interface for single triangle intersection
 * tests when the intersection point is always needed.
 */
std::pair<bool, glm::dvec3> Ray::intersect_triangle(Triangle& triangle) const {
	std::pair<bool, glm::dvec3> result;
	// Use the main intersection method and package result
	result.first = intersect_triangle(triangle, result.second);
	return result;
}

/**
 * @brief Ray-plane intersection using vectorized operations
 * @param normal Plane normal vector (will be normalized internally)
 * @param point Any point on the plane
 * @param intersection Output intersection point if found
 * @return True if intersection exists and is in front of ray origin
 *
 * Uses GLM vectorized operations for high performance plane intersection.
 * Handles parallel ray cases and ensures intersection is forward along ray.
 */
bool Ray::intersect_plane(const glm::dvec3& normal, const glm::dvec3& point, glm::dvec3& intersection) const noexcept {
	// Normalize input normal for accurate geometric calculations
	const glm::dvec3 norm = glm::normalize(normal);

	// Calculate ray-plane angle via dot product
	const double denominator = glm::dot(direction_, norm);

	// Check for parallel ray (no intersection possible)
	if (std::abs(denominator) < MathConstants::GEOMETRIC_EPSILON) {
		return false;
	}

	// Calculate parametric distance to intersection point
	const glm::dvec3 origin_to_point = point - origin_;
	const double t = glm::dot(origin_to_point, norm) / denominator;

	// Reject intersections behind ray origin (using epsilon for stability)
	if (t <= MathConstants::GEOMETRIC_EPSILON) {
		return false;
	}

	// Compute final intersection point using parametric ray equation
	intersection = origin_ + t * direction_;
	return true;
}

/**
 * @brief Find all ray-triangle intersections in a collection
 * @param triangles Span of triangles to test against
 * @param intersections Output vector of intersection points
 *
 * Efficiently tests ray against multiple triangles and collects all
 * intersection points. Uses span for zero-overhead iteration and
 * reserves memory based on expected hit rate.
 */
void Ray::intersect_triangles(std::span<Triangle> triangles, std::vector<glm::dvec3>& intersections) const {
	// Clear previous results and reserve memory for efficiency
	intersections.clear();
	intersections.reserve(triangles.size() / 10); // Assume ~10% hit rate for typical scenes

	// Test each triangle and collect intersection points
	for (Triangle& triangle : triangles) {
		auto result = intersect_triangle(triangle);
		if (result.first) {
			intersections.emplace_back(result.second);
		}
	}
}

/**
 * @brief Find closest ray-triangle intersection in a collection
 * @param triangles Span of triangles to test
 * @param intersection Output closest intersection point
 * @param hit_triangle Output reference to the intersected triangle
 * @return Distance to closest intersection, or -1.0 if no intersection
 *
 * Performs distance-sorted intersection finding for ray-mesh collision.
 * Essential for proper visibility determination in ray tracing applications.
 */
double Ray::intersect_first_triangle(std::span<Triangle> triangles,
									 glm::dvec3& intersection,
									 Triangle& hit_triangle) const {
	// Initialize search for closest intersection
	double shortest_dist = std::numeric_limits<double>::max();
	bool found_intersection = false;

	// Test all triangles and track closest hit
	for (Triangle& triangle : triangles) {
		glm::dvec3 temp_intersection;
		if (intersect_triangle(triangle, temp_intersection)) {
			// Calculate distance to intersection point
			double dist = glm::distance(origin_, temp_intersection);

			// Set closest hit if this is nearer
			if (dist < shortest_dist) {
				shortest_dist = dist;
				intersection = temp_intersection;
				hit_triangle = triangle;
				found_intersection = true;
			}
		}
	}

	// Return distance to closest hit, or -1.0 for no intersection
	return found_intersection ? shortest_dist : -1.0;
}

/**
 * @brief Find closest intersection with a material layer
 * @param layer Material layer containing mesh geometry
 * @param intersection Output closest intersection point
 * @return Distance to closest intersection, or -1.0 if no intersection
 *
 * Specialized intersection test for material layers in multi-material
 * simulations. Handles layer-specific geometry efficiently.
 */
double Ray::intersect_layer(const Layer& layer, glm::dvec3& intersection) const {
	// Initialize closest distance search
	double shortest_dist = std::numeric_limits<double>::max();
	bool found_intersection = false;

	// Test all triangles in layer mesh
	for (const auto& triangle : layer.mesh) {
		glm::dvec3 temp_intersection;
		if (intersect_triangle(const_cast<Triangle&>(triangle), temp_intersection)) {
			// Check if this intersection is closer
			const double dist = glm::distance(origin_, temp_intersection);
			if (dist < shortest_dist) {
				shortest_dist = dist;
				intersection = temp_intersection;
				found_intersection = true;
			}
		}
	}

	// Return closest distance or indicate no intersection
	return found_intersection ? shortest_dist : -1.0;
}

/**
 * @brief Find closest intersection across multiple material layers
 * @param layers Span of material layers to test
 * @param intersection Output closest intersection point
 * @param hit_triangle Output reference to intersected triangle
 * @return Distance to closest intersection, or -1.0 if no intersection
 *
 * Performs hierarchical intersection testing across multiple material layers.
 * Critical for multi-material Monte Carlo simulations where ray traversal
 * must respect material boundaries.
 */
double Ray::intersect_first_triangle(std::span<const Layer> layers,
									 glm::dvec3& intersection,
									 Triangle& hit_triangle) const {
	// Initialize search across all layers
	double shortest_dist = std::numeric_limits<double>::max();
	bool found_intersection = false;

	// Iterate through all material layers
	for (const auto& layer : layers) {
		// Test each triangle in current layer
		for (const auto& triangle : layer.mesh) {
			glm::dvec3 temp_intersection;
			if (intersect_triangle(const_cast<Triangle&>(triangle), temp_intersection)) {
				// Calculate distance and set closest hit
				const double dist = glm::distance(origin_, temp_intersection);
				if (dist < shortest_dist) {
					shortest_dist = dist;
					intersection = temp_intersection;
					hit_triangle = triangle;
					found_intersection = true;
				}
			}
		}
	}

	// Return distance to closest hit across all layers
	return found_intersection ? shortest_dist : -1.0;
}

/**
 * @brief Intersect ray with cuboid geometry (internal surface)
 * @param cuboid Target cuboid for intersection test
 * @param intersection Output intersection point
 * @param normal Output surface normal at intersection
 * @return Distance to intersection, or negative if no intersection
 *
 * Delegates to cuboid's internal intersection method for proper
 * geometric handling of ray-box intersection from inside the volume.
 */
double Ray::intersect_cuboid_internal(Cuboid& cuboid, glm::dvec3& intersection, glm::dvec3& normal) const {
	// Use cuboid's optimized internal intersection algorithm
	return cuboid.intersect_ray_internal(*this, intersection, normal);
}
