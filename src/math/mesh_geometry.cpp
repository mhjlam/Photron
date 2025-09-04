#include "mesh_geometry.hpp"

#include <algorithm>
#include <limits>

/**
 * Calculate the bounding box of the mesh
 */
void MeshGeometry::calculate_bounds() {
	if (triangles_.empty()) {
		has_bounds_ = false;
		return;
	}

	double min_x = std::numeric_limits<double>::max();
	double min_y = std::numeric_limits<double>::max();
	double min_z = std::numeric_limits<double>::max();

	double max_x = -std::numeric_limits<double>::max();
	double max_y = -std::numeric_limits<double>::max();
	double max_z = -std::numeric_limits<double>::max();

	// Find the min and max coordinates across all vertices
	for (const auto& triangle : triangles_) {
		const glm::dvec3& v0 = triangle.v0();
		const glm::dvec3& v1 = triangle.v1();
		const glm::dvec3& v2 = triangle.v2();

		// Update min bounds
		min_x = std::min(min_x, std::min(v0.x, std::min(v1.x, v2.x)));
		min_y = std::min(min_y, std::min(v0.y, std::min(v1.y, v2.y)));
		min_z = std::min(min_z, std::min(v0.z, std::min(v1.z, v2.z)));

		// Update max bounds
		max_x = std::max(max_x, std::max(v0.x, std::max(v1.x, v2.x)));
		max_y = std::max(max_y, std::max(v0.y, std::max(v1.y, v2.y)));
		max_z = std::max(max_z, std::max(v0.z, std::max(v1.z, v2.z)));
	}

	// Set the bounds
	bounds_.min_bounds = glm::dvec3(min_x, min_y, min_z);
	bounds_.max_bounds = glm::dvec3(max_x, max_y, max_z);
	has_bounds_ = true;
}

/**
 * Extend the current bounds to include a triangle
 */
void MeshGeometry::extend_bounds(const Triangle& triangle) {
	if (!has_bounds_) {
		calculate_bounds();
		return;
	}

	const glm::dvec3& v0 = triangle.v0();
	const glm::dvec3& v1 = triangle.v1();
	const glm::dvec3& v2 = triangle.v2();

	// Update min bounds
	bounds_.min_bounds.x = std::min(bounds_.min_bounds.x, std::min(v0.x, std::min(v1.x, v2.x)));
	bounds_.min_bounds.y = std::min(bounds_.min_bounds.y, std::min(v0.y, std::min(v1.y, v2.y)));
	bounds_.min_bounds.z = std::min(bounds_.min_bounds.z, std::min(v0.z, std::min(v1.z, v2.z)));

	// Update max bounds
	bounds_.max_bounds.x = std::max(bounds_.max_bounds.x, std::max(v0.x, std::max(v1.x, v2.x)));
	bounds_.max_bounds.y = std::max(bounds_.max_bounds.y, std::max(v0.y, std::max(v1.y, v2.y)));
	bounds_.max_bounds.z = std::max(bounds_.max_bounds.z, std::max(v0.z, std::max(v1.z, v2.z)));
}

/**
 * Check if a point is inside the mesh using BVH-accelerated ray casting
 */
bool MeshGeometry::contains(const glm::dvec3& point) const {
	if (triangles_.empty() || !might_contain(point)) {
		return false;
	}

	// Use BVH for fast point-in-mesh testing
	return bvh_.contains_point(point);
}

/**
 * Find the first intersection of a ray with the mesh using BVH
 * Returns the distance to intersection and the triangle hit
 */
std::pair<double, Triangle> MeshGeometry::find_first_intersection(const Ray& ray) const {
	if (triangles_.empty()) {
		return std::make_pair(std::numeric_limits<double>::max(), Triangle());
	}

	double distance;
	glm::dvec3 intersection;
	Triangle hit_triangle;

	if (bvh_.intersect(ray, distance, intersection, hit_triangle)) {
		return std::make_pair(distance, hit_triangle);
	} else {
		// Return max distance and an invalid triangle to indicate no intersection
		return std::make_pair(std::numeric_limits<double>::max(), Triangle());
	}
}

/**
 * Find all intersections of a ray with the mesh (uses BVH for initial filtering)
 */
std::vector<std::pair<double, Triangle>> MeshGeometry::find_all_intersections(const Ray& ray) const {
	std::vector<std::pair<double, Triangle>> intersections;

	if (triangles_.empty()) {
		return intersections;
	}

	// For now, fall back to brute force for all intersections
	// TODO: Implement BVH traversal that finds all intersections
	for (const auto& triangle : triangles_) {
		Triangle triangle_copy = triangle; // Create mutable copy
		glm::dvec3 intersection;

		if (ray.intersect_triangle(triangle_copy, intersection)) {
			// Calculate distance to intersection
			double dist = glm::distance(ray.origin(), intersection);
			intersections.push_back(std::make_pair(dist, triangle));
		}
	}

	// Sort intersections by distance
	std::sort(intersections.begin(), intersections.end(),
			  [](const auto& a, const auto& b) { return a.first < b.first; });

	return intersections;
}
