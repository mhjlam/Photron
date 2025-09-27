#include "layer.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <ranges>
#include <sstream>

#include "logger.hpp"

/**
 * Calculate the bounding box of the mesh
 */
void Layer::calculate_bounds() {
	if (mesh.empty()) {
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
	for (const auto& triangle : mesh) {
		const glm::dvec3& v0 = triangle.v0();
		const glm::dvec3& v1 = triangle.v1();
		const glm::dvec3& v2 = triangle.v2();

		// Update min bounds using C++20 ranges
		min_x = std::min(min_x, std::ranges::min({v0.x, v1.x, v2.x}));
		min_y = std::min(min_y, std::ranges::min({v0.y, v1.y, v2.y}));
		min_z = std::min(min_z, std::ranges::min({v0.z, v1.z, v2.z}));

		// Update max bounds using C++20 ranges
		max_x = std::max(max_x, std::ranges::max({v0.x, v1.x, v2.x}));
		max_y = std::max(max_y, std::ranges::max({v0.y, v1.y, v2.y}));
		max_z = std::max(max_z, std::ranges::max({v0.z, v1.z, v2.z}));
	}

	// Set the bounds
	bounds_.min_bounds = glm::dvec3(min_x, min_y, min_z);
	bounds_.max_bounds = glm::dvec3(max_x, max_y, max_z);
	has_bounds_ = true;
}

/**
 * Extend the current bounds to include a triangle
 */
void Layer::extend_bounds(const Triangle& triangle) {
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
bool Layer::contains_point(const glm::dvec3& point) const {
	if (mesh.empty() || !might_contain(point)) {
		return false;
	}

	// Use BVH for fast point-in-mesh testing
	return bvh_.contains_point(point);
}

/**
 * Find the first intersection of a ray with the mesh using BVH
 * Returns the distance to intersection and the triangle hit
 */
std::pair<double, Triangle> Layer::find_first_intersection(const Ray& ray) const {
	if (mesh.empty()) {
		return std::make_pair(std::numeric_limits<double>::max(), Triangle());
	}

	double distance;
	glm::dvec3 intersection;
	Triangle hit_triangle;

	if (bvh_.intersect(ray, distance, intersection, hit_triangle)) {
		return std::make_pair(distance, hit_triangle);
	}
	else {
		// Return max distance and an invalid triangle to indicate no intersection
		return std::make_pair(std::numeric_limits<double>::max(), Triangle());
	}
}

/**
 * Find all intersections of a ray with the mesh (uses BVH for initial filtering)
 */
std::vector<std::pair<double, Triangle>> Layer::find_all_intersections(const Ray& ray) const {
	if (mesh.empty()) {
		return {};
	}

	// Modern C++20: Use ranges and transform_reduce for performance
	std::vector<std::pair<double, Triangle>> intersections;
	intersections.reserve(mesh.size() / 10); // Reserve space for ~10% hit rate

	// Use ranges::for_each with lambda capture for better performance
	std::ranges::for_each(mesh, [&](const Triangle& triangle) {
		Triangle triangle_copy = triangle; // Create mutable copy
		glm::dvec3 intersection_point;

		if (ray.intersect_triangle(triangle_copy, intersection_point)) {
			// Calculate distance to intersection - use squared distance for performance
			const double dist = glm::distance(ray.origin(), intersection_point);
			intersections.emplace_back(dist, triangle);
		}
	});

	// Modern C++20: Use ranges::sort for cleaner syntax
	std::ranges::sort(intersections, [](const auto& a, const auto& b) noexcept { return a.first < b.first; });

	return intersections;
}

/**
 * Validate and fix normal orientations to ensure they point outward from the geometry
 */
void Layer::validate_and_fix_normals() {
	if (mesh.empty()) {
		return;
	}

	// Calculate the centroid of the mesh for reference
	glm::dvec3 mesh_centroid(0.0);
	size_t vertex_count = 0;
	
	for (const auto& triangle : mesh) {
		mesh_centroid += triangle.v0() + triangle.v1() + triangle.v2();
		vertex_count += 3;
	}
	
	if (vertex_count > 0) {
		mesh_centroid /= static_cast<double>(vertex_count);
	}

	// Check and fix each triangle's normal orientation
	size_t flipped_count = 0;
	for (auto& triangle : mesh) {
		// Calculate triangle center
		glm::dvec3 triangle_center = triangle.center();
		
		// Vector from mesh centroid to triangle center (should align with outward normal)
		glm::dvec3 centroid_to_triangle = triangle_center - mesh_centroid;
		
		// Get current normal
		glm::dvec3 current_normal = triangle.normal();
		
		// Check if normal points outward (positive dot product with centroid-to-triangle vector)
		double alignment = glm::dot(current_normal, glm::normalize(centroid_to_triangle));
		
		// If normal points inward (negative dot product), flip it
		if (alignment < 0.0) {
			// Flip the normal by swapping two vertices (reverses winding order)
			glm::dvec3 v0 = triangle.v0();
			glm::dvec3 v1 = triangle.v1();
			glm::dvec3 v2 = triangle.v2();
			
			// Swap v1 and v2 to reverse winding order
			triangle.set_vertices(v0, v2, v1);
			flipped_count++;
		}
	}

	if (flipped_count > 0) {
		std::ostringstream debug_msg;
		debug_msg << "Layer " << static_cast<int>(id) << ": Fixed " << flipped_count 
				  << " inward-pointing normals out of " << mesh.size() << " triangles.";
		Logger::instance().log_info(debug_msg.str());
	}

	// Rebuild spatial acceleration structure with corrected normals
	bvh_.build(mesh);
}
