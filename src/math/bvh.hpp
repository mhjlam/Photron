#pragma once

#include <memory>
#include <vector>

#include "concepts.hpp"
#include "cuboid.hpp"
#include "math.hpp"
#include "ray.hpp"
#include "triangle.hpp"

/**
 * BVH Node - can be either internal or leaf
 * Uses unified Cuboid class for bounding boxes
 */
struct BVHNode
{
	Cuboid bounds;
	std::unique_ptr<BVHNode> left;
	std::unique_ptr<BVHNode> right;
	std::vector<int> triangle_indices; // Only used for leaf nodes

	[[nodiscard]] bool is_leaf() const noexcept { return left == nullptr && right == nullptr; }
};

/**
 * Bounding Volume Hierarchy for fast ray-triangle intersection
 * and point-in-mesh testing
 */
class BVH
{
private:
	std::unique_ptr<BVHNode> root_;
	std::vector<Triangle> triangles_;

	// Build the BVH recursively
	std::unique_ptr<BVHNode> build_recursive(std::vector<int>& triangle_indices, int depth = 0);

	// Calculate Cuboid bounds for a set of triangles
	[[nodiscard]] Cuboid calculate_bounds(const std::vector<int>& triangle_indices) const;

	// Find best split using Surface Area Heuristic (SAH)
	[[nodiscard]] int find_best_split(std::vector<int>& triangle_indices, int& best_axis) const;

	// Ray intersection traversal
	[[nodiscard]] bool intersect_recursive(const BVHNode* node, const Ray& ray, double& closest_t,
										   glm::dvec3& intersection, Triangle& hit_triangle) const;

	// Point-in-mesh testing traversal
	[[nodiscard]] int count_intersections_recursive(const BVHNode* node, const Ray& ray, const glm::dvec3& point) const;

public:
	// Maximum triangles per leaf node
	static constexpr int MAX_TRIANGLES_PER_LEAF = 4;

	// Maximum tree depth
	static constexpr int MAX_DEPTH = 20;

	BVH() = default;

	explicit BVH(const std::vector<Triangle>& triangles);

	// Delete copy constructor and copy assignment (non-copyable due to unique_ptr)
	BVH(const BVH&) = delete;
	BVH& operator=(const BVH&) = delete;

	// Enable move constructor and move assignment
	BVH(BVH&&) = default;
	BVH& operator=(BVH&&) = default;

	// Build the BVH from triangles
	void build(const std::vector<Triangle>& triangles);

	// Find closest ray-triangle intersection
	[[nodiscard]] bool intersect(const Ray& ray, double& distance, glm::dvec3& intersection,
								 Triangle& hit_triangle) const;

	// Test if point is inside the mesh using ray casting
	[[nodiscard]] bool contains_point(const glm::dvec3& point) const;

	// Get all triangles (for debugging)
	[[nodiscard]] const std::vector<Triangle>& triangles() const noexcept { return triangles_; }

	// Check if BVH is empty
	[[nodiscard]] constexpr bool empty() const noexcept { return triangles_.empty() || !root_; }
};
