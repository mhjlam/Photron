/**
 * @file layer.cpp
 * @brief Implementation of geometric layer processing for multi-material simulations
 *
 * Implements the Layer class which represents a single material layer within
 * a multi-layered medium. Handles triangle mesh processing, bounding volume
 * calculations, and provides BVH-accelerated geometric queries for efficient
 * point-in-mesh testing and ray-triangle intersections.
 */

#include "layer.hpp"

// Add includes that were removed from the header
#include <algorithm>
#include <iostream>
#include <limits>
#include <ranges>
#include <sstream>

#include "common/logger.hpp"
#include "math/math.hpp"
#include "math/ray.hpp"
#include "simulator/layer.hpp"
#include "simulator/material.hpp"

void Layer::calculate_bounds() {
	// Calculate axis-aligned bounding box from triangle vertices
	if (mesh.empty()) {
		has_bounds_ = false;
		return;
	}

	// Initialize bounds for min/max search
	double min_x = std::numeric_limits<double>::max();
	double min_y = std::numeric_limits<double>::max();
	double min_z = std::numeric_limits<double>::max();

	double max_x = -std::numeric_limits<double>::max();
	double max_y = -std::numeric_limits<double>::max();
	double max_z = -std::numeric_limits<double>::max();

	// Process all triangle vertices to find coordinate extrema
	for (const auto& triangle : mesh) {
		const glm::dvec3& v0 = triangle.v0();
		const glm::dvec3& v1 = triangle.v1();
		const glm::dvec3& v2 = triangle.v2();

		// Update minimum bounds using modern C++20 ranges algorithms
		min_x = std::min(min_x, std::ranges::min({v0.x, v1.x, v2.x}));
		min_y = std::min(min_y, std::ranges::min({v0.y, v1.y, v2.y}));
		min_z = std::min(min_z, std::ranges::min({v0.z, v1.z, v2.z}));

		// Update maximum bounds using modern C++20 ranges algorithms
		max_x = std::max(max_x, std::ranges::max({v0.x, v1.x, v2.x}));
		max_y = std::max(max_y, std::ranges::max({v0.y, v1.y, v2.y}));
		max_z = std::max(max_z, std::ranges::max({v0.z, v1.z, v2.z}));
	}

	// Store calculated bounds and mark as valid
	bounds_.min_bounds = glm::dvec3(min_x, min_y, min_z);
	bounds_.max_bounds = glm::dvec3(max_x, max_y, max_z);
	has_bounds_ = true;
}

void Layer::extend_bounds(const Triangle& triangle) {
	// Initialize bounds if not yet calculated
	if (!has_bounds_) {
		calculate_bounds();
		return;
	}

	// Extract triangle vertices for bounds extension
	const glm::dvec3& v0 = triangle.v0();
	const glm::dvec3& v1 = triangle.v1();
	const glm::dvec3& v2 = triangle.v2();

	// Extend minimum bounds to include new triangle
	bounds_.min_bounds.x = std::min(bounds_.min_bounds.x, std::min(v0.x, std::min(v1.x, v2.x)));
	bounds_.min_bounds.y = std::min(bounds_.min_bounds.y, std::min(v0.y, std::min(v1.y, v2.y)));
	bounds_.min_bounds.z = std::min(bounds_.min_bounds.z, std::min(v0.z, std::min(v1.z, v2.z)));

	// Extend maximum bounds to include new triangle
	bounds_.max_bounds.x = std::max(bounds_.max_bounds.x, std::max(v0.x, std::max(v1.x, v2.x)));
	bounds_.max_bounds.y = std::max(bounds_.max_bounds.y, std::max(v0.y, std::max(v1.y, v2.y)));
	bounds_.max_bounds.z = std::max(bounds_.max_bounds.z, std::max(v0.z, std::max(v1.z, v2.z)));
}

bool Layer::contains_point(const glm::dvec3& point) const {
	// Quick bounds check before expensive point-in-mesh test
	if (mesh.empty() || !might_contain(point)) {
		return false;
	}

	// BVH-accelerated point containment test
	return bvh_.contains_point(point);
}

std::pair<double, Triangle> Layer::find_first_intersection(const Ray& ray) const {
	// Handle empty mesh case
	if (mesh.empty()) {
		return std::make_pair(std::numeric_limits<double>::max(), Triangle());
	}

	// BVH-accelerated ray intersection query
	double distance;
	glm::dvec3 intersection;
	Triangle hit_triangle;

	if (bvh_.intersect(ray, distance, intersection, hit_triangle)) {
		return std::make_pair(distance, hit_triangle);
	}
	else {
		// No intersection found
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
 * Ensure all triangle normals point outward from mesh centroid
 */
void Layer::validate_and_fix_normals() {
	if (mesh.empty()) {
		return;
	}

	// Compute mesh centroid as reference point
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

		// Test normal orientation relative to mesh center
		double alignment = glm::dot(current_normal, glm::normalize(centroid_to_triangle));

		// Fix inward-pointing normals by reversing winding order
		if (alignment < 0.0) {
			// Swap vertices to flip normal direction
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

	// Rebuild BVH with corrected triangle orientations
	bvh_.build(mesh);
}

void Layer::update_geometry() noexcept {
	calculate_bounds();
	bvh_.build(mesh);
}

void Layer::add_triangle(const Triangle& triangle) {
	mesh.push_back(triangle);
	if (has_bounds_) {
		// Update bounds
		extend_bounds(triangle);
	}
	// Rebuild BVH when triangles are added
	bvh_.build(mesh);
}
